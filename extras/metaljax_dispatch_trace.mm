// Optional diagnostics for the checksum-pinned private MLX runtime.
// Uses Objective-C method hooks; adds no GPU commands or synchronization.
#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#import <objc/runtime.h>
#include <CommonCrypto/CommonDigest.h>
#include <mlx/stream.h>
#include <mlx/backend/metal/metal.h>
#include <array>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <dlfcn.h>
#include <map>
#include <mutex>
#include <set>
#include <sstream>
#include <string>
#include <time.h>
#include <thread>
#include <chrono>

namespace mlx::core::metal {
class CommandEncoder;
CommandEncoder& get_command_encoder(mlx::core::Stream);
}
namespace {
namespace mx = mlx::core;
uint64_t now() { return clock_gettime_nsec_np(CLOCK_UPTIME_RAW); }
std::string quote(const std::string& s) {
  std::string out = "\"";
  for (unsigned char c : s) {
    if (c == '\\' || c == '"') out += '\\';
    if (c >= 32) out += c;
  }
  return out + '"';
}
std::string phase() {
  const char* s = std::getenv("FBC_PROFILE_ACTIVE");
  return s && s[0] != '0' ? s : "";
}
bool copies() { static bool on = std::getenv("FBC_MLX_TRACE_COPIES") != nullptr; return on; }
bool reflectionEnabled() { static bool on = std::getenv("FBC_MLX_TRACE_REFLECTION") != nullptr; return on; }
struct Binding { uintptr_t pointer = 0; uint64_t bytes = 0, offset = 0; };
struct Encoder {
  uintptr_t buffer = 0;
  std::string pipeline;
  std::array<Binding, 32> bindings{};
  std::map<unsigned, std::string> constants;
  bool collect_sample = false;
  bool capture_copy = false;
};
struct Kernel {
  uint64_t count = 0, cpu_dispatch_ns = 0, binding_span_bytes = 0;
  std::set<std::string> samples;
  std::string arguments = "[]";
};
struct Buffer {
  uint64_t dispatches = 0, encode_begin_ns = 0;
  std::map<std::string, uint64_t> kernels;
};
struct State {
  std::mutex mutex;
  std::map<uintptr_t, Encoder> encoders;
  std::map<uintptr_t, Buffer> buffers;
  std::map<std::pair<std::string, std::string>, Kernel> kernels;
  std::map<std::pair<Class, SEL>, IMP> originals;
  std::map<uintptr_t, std::string> pipeline_names;
  std::map<uintptr_t, std::string> pipeline_arguments;
  FILE* output = nullptr;
  uint64_t serial = 0;
  uint64_t pending = 0;
};
State& state() { static State* value = new State; return *value; }
template<class F> F original(id object, SEL selector) {
  return reinterpret_cast<F>(state().originals.at({object_getClass(object), selector}));
}
void hook(Class cls, SEL selector, IMP replacement) {
  auto key = std::make_pair(cls, selector);
  if (state().originals.count(key)) return;
  Method method = class_getInstanceMethod(cls, selector);
  if (!method) return;
  state().originals[key] = method_getImplementation(method);
  if (!class_addMethod(cls, selector, replacement, method_getTypeEncoding(method)))
    method_setImplementation(class_getInstanceMethod(cls, selector), replacement);
}

// Read-only prefix of the frozen CommandEncoder. Offsets 0x18/0x20/0x28/
// 0x30 are checked against recorded disassembly. Refuse every other MLX image.
// This is not a portable MLX API, and is used only for diagnostic snapshots.
struct Prefix {
  void* device; bool exiting; unsigned char padding[7];
  void* queue; void* buffer; int operations; unsigned padding2;
  uint64_t input_bytes, output_bytes;
};
static_assert(offsetof(Prefix, buffer) == 0x18);
static_assert(offsetof(Prefix, operations) == 0x20);
static_assert(offsetof(Prefix, input_bytes) == 0x28);
static_assert(offsetof(Prefix, output_bytes) == 0x30);
Prefix prefix(id buffer, bool& matched) {
  auto& encoder = mx::metal::get_command_encoder(mx::default_stream(mx::Device::gpu));
  Prefix result;
  std::memcpy(&result, &encoder, sizeof result);
  matched = result.buffer == (void*)buffer;
  return result;
}
long setting(const char* name, long fallback) {
  const char* s = std::getenv(name); return s ? std::strtol(s, nullptr, 10) : fallback;
}
void commit(id self, SEL selector) {
  using Function = void(*)(id, SEL);
  std::string active = phase();
  if (active.empty()) { original<Function>(self, selector)(self, selector); return; }
  bool matched = false;
  Prefix p = prefix(self, matched);
  uint64_t serial;
  Buffer recorded;
  {
    std::lock_guard<std::mutex> lock(state().mutex);
    serial = ++state().serial;
    ++state().pending;
    recorded = state().buffers[(uintptr_t)self];
    state().buffers.erase((uintptr_t)self);
  }
  const bool transient = setting("MLX_MAX_MB_TRANSIENT_ONLY", 0) != 0;
  const uint64_t charged = transient ? p.output_bytes : p.input_bytes;
  std::string reason = !matched ? "unmatched" :
      p.operations > setting("MLX_MAX_OPS_PER_BUFFER", 800) ? "operations_threshold" :
      (charged >> 20) > (uint64_t)setting("MLX_MAX_MB_PER_BUFFER", 512) ? "bytes_threshold" : "tail_or_explicit";
  uint64_t host = now();
  [self addCompletedHandler:^(id<MTLCommandBuffer> completed) {
    std::lock_guard<std::mutex> lock(state().mutex);
    if (state().output)
      std::fprintf(state().output, "{\"type\":\"buffer\",\"phase\":%s,\"id\":%llu,\"host_commit_ns\":%llu,\"gpu_start_ns\":%llu,\"gpu_end_ns\":%llu,\"dispatches\":%llu,\"matched\":%s,\"native_operations\":%d,\"input_bytes\":%llu,\"output_bytes\":%llu,\"threshold\":%s,\"host_encode_begin_ns\":%llu,\"kernel_counts\":{",
                   quote(active).c_str(), serial, host,
                   (uint64_t)(completed.GPUStartTime * 1e9),
                   (uint64_t)(completed.GPUEndTime * 1e9), recorded.dispatches,
                   matched ? "true" : "false", matched ? p.operations : -1,
                   matched ? p.input_bytes : 0, matched ? p.output_bytes : 0,
                   quote(reason).c_str(), recorded.encode_begin_ns);
    bool comma = false;
    for (const auto& [name, count] : recorded.kernels) {
      if (comma) std::fputc(',', state().output); comma = true;
      std::fprintf(state().output, "%s:%llu", quote(name).c_str(), count);
    }
    std::fputs("}}\n", state().output);
    --state().pending;
  }];
  original<Function>(self, selector)(self, selector);
}
id makeEncoder(id self, SEL selector, MTLDispatchType type) {
  using Function = id(*)(id, SEL, MTLDispatchType);
  const uint64_t start = now();
  id encoder = original<Function>(self, selector)(self, selector, type);
  if (!phase().empty()) {
    std::lock_guard<std::mutex> lock(state().mutex);
    state().encoders[(uintptr_t)encoder] = Encoder{(uintptr_t)self};
    auto& record = state().buffers[(uintptr_t)self];
    if (!record.encode_begin_ns) record.encode_begin_ns = start;
  }
  return encoder;
}
void pipeline(id self, SEL selector, id<MTLComputePipelineState> value) {
  using Function = void(*)(id, SEL, id<MTLComputePipelineState>);
  if (!phase().empty()) {
    std::lock_guard<std::mutex> lock(state().mutex);
    const char* label = value.label.UTF8String;
    auto& encoder = state().encoders[(uintptr_t)self];
    encoder.pipeline = label ? label : state().pipeline_names[(uintptr_t)value];
    if (encoder.pipeline.empty()) encoder.pipeline = "unlabelled";
    encoder.bindings = {};
    encoder.constants.clear();
    encoder.capture_copy = copies() && (encoder.pipeline.find("copy") != std::string::npos ||
        encoder.pipeline.find("steel_gemm") != std::string::npos ||
        encoder.pipeline.find("gemv") != std::string::npos);
    auto& record = state().kernels[{phase(), encoder.pipeline}];
    encoder.collect_sample = record.count < 4 || encoder.capture_copy;
    if (state().pipeline_arguments.count((uintptr_t)value)) record.arguments = state().pipeline_arguments[(uintptr_t)value];
  }
  original<Function>(self, selector)(self, selector, value);
}
void rememberPipeline(id result, NSString* name, MTLComputePipelineReflection* reflection) {
  if (!result) return;
  std::lock_guard<std::mutex> lock(state().mutex);
  state().pipeline_names[(uintptr_t)result] = name.UTF8String ?: "unlabelled";
  if (reflection) {
    std::ostringstream out; out << '['; bool comma = false;
    for (MTLArgument* a in reflection.arguments) {
      if (a.type != MTLArgumentTypeBuffer) continue;
      if (comma) out << ','; comma = true;
      out << "{\"name\":" << quote(a.name.UTF8String ?: "") << ",\"index\":" << a.index
          << ",\"access\":" << a.access << ",\"data_type\":" << a.bufferDataType
          << ",\"data_size\":" << a.bufferDataSize << ",\"members\":[";
      bool member_comma = false;
      if (a.bufferDataType == MTLDataTypeStruct) for (MTLStructMember* member in a.bufferStructType.members) {
        if (member_comma) out << ','; member_comma = true;
        out << "{\"name\":" << quote(member.name.UTF8String ?: "")
            << ",\"offset\":" << member.offset << ",\"data_type\":" << member.dataType << '}';
      }
      out << "]}";
    }
    out << ']'; state().pipeline_arguments[(uintptr_t)result] = out.str();
  }
}
id newPipeline(id self, SEL selector, id<MTLFunction> function, NSError** error) {
  using Function = id(*)(id, SEL, id<MTLFunction>, NSError**);
  MTLComputePipelineReflection* reflection = nil;
  id result = reflectionEnabled() ? [self newComputePipelineStateWithFunction:function
      options:(MTLPipelineOptionArgumentInfo | MTLPipelineOptionBufferTypeInfo) reflection:&reflection error:error]
      : original<Function>(self, selector)(self, selector, function, error);
  rememberPipeline(result, function.name, reflection);
  return result;
}
id newPipelineDescriptor(id self, SEL selector, MTLComputePipelineDescriptor* descriptor,
                         MTLPipelineOption options, MTLComputePipelineReflection** reflection, NSError** error) {
  using Function = id(*)(id, SEL, MTLComputePipelineDescriptor*, MTLPipelineOption, MTLComputePipelineReflection**, NSError**);
  MTLComputePipelineReflection* metadata = nil;
  id result = original<Function>(self, selector)(self, selector, descriptor,
      reflectionEnabled() ? options | MTLPipelineOptionArgumentInfo | MTLPipelineOptionBufferTypeInfo : options,
      reflectionEnabled() ? &metadata : reflection, error);
  if (reflectionEnabled() && reflection) *reflection = metadata;
  rememberPipeline(result, descriptor.computeFunction.name, metadata);
  return result;
}
void binding(id self, SEL selector, id<MTLBuffer> buffer, NSUInteger offset, NSUInteger index) {
  using Function = void(*)(id, SEL, id<MTLBuffer>, NSUInteger, NSUInteger);
  if (index < 32 && !phase().empty()) {
    std::lock_guard<std::mutex> lock(state().mutex);
    state().encoders[(uintptr_t)self].bindings[index] = {(uintptr_t)buffer, buffer.length, offset};
  }
  original<Function>(self, selector)(self, selector, buffer, offset, index);
}
void constants(id self, SEL selector, const void* bytes, NSUInteger length, NSUInteger index) {
  using Function = void(*)(id, SEL, const void*, NSUInteger, NSUInteger);
  bool sample = false;
  if (!phase().empty()) {
    std::lock_guard<std::mutex> lock(state().mutex);
    sample = state().encoders[(uintptr_t)self].collect_sample;
  }
  if (sample) {
    std::string hex; const auto* data = static_cast<const unsigned char*>(bytes);
    for (size_t i = 0; i < std::min<NSUInteger>(length, 128); ++i) {
      hex += "0123456789abcdef"[data[i] >> 4]; hex += "0123456789abcdef"[data[i] & 15];
    }
    std::lock_guard<std::mutex> lock(state().mutex);
    state().encoders[(uintptr_t)self].constants[index] = hex;
  }
  original<Function>(self, selector)(self, selector, bytes, length, index);
}
void dispatch(id self, SEL selector, MTLSize grid, MTLSize group) {
  using Function = void(*)(id, SEL, MTLSize, MTLSize);
  std::string active = phase();
  if (active.empty()) { original<Function>(self, selector)(self, selector, grid, group); return; }
  const uint64_t start = now();
  original<Function>(self, selector)(self, selector, grid, group);
  const uint64_t elapsed = now() - start;
  std::lock_guard<std::mutex> lock(state().mutex);
  auto& e = state().encoders[(uintptr_t)self];
  ++state().buffers[e.buffer].dispatches;
  ++state().buffers[e.buffer].kernels[e.pipeline];
  auto& k = state().kernels[{active, e.pipeline}];
  ++k.count; k.cpu_dispatch_ns += elapsed;
  uint64_t bound = 0;
  for (const auto& b : e.bindings) if (b.pointer && b.bytes >= b.offset) bound += b.bytes - b.offset;
  k.binding_span_bytes += bound;
  if (k.count <= 4 || e.capture_copy) {
    std::ostringstream sample;
    sample << "{\"grid\":[" << grid.width << ',' << grid.height << ',' << grid.depth
           << "],\"group\":[" << group.width << ',' << group.height << ',' << group.depth
           << "],\"grid_is_threadgroups\":" << (selector == @selector(dispatchThreadgroups:threadsPerThreadgroup:) ? "true" : "false") << ",\"bindings\":[";
    bool comma = false;
    for (size_t i = 0; i < e.bindings.size(); ++i) {
      const auto& b = e.bindings[i]; if (!b.pointer) continue;
      if (comma) sample << ','; comma = true;
      sample << "{\"index\":" << i << ",\"resource\":" << b.pointer << ",\"bytes\":" << b.bytes << ",\"offset\":" << b.offset << '}';
    }
    sample << "],\"constants_hex\":{"; comma = false;
    for (const auto& [index, hex] : e.constants) {
      if (comma) sample << ','; comma = true;
      sample << quote(std::to_string(index)) << ':' << quote(hex);
    }
    sample << "}}";
    if (k.count <= 4) k.samples.insert(sample.str());
    if (e.capture_copy && state().output)
      std::fprintf(state().output, "{\"type\":\"operation\",\"phase\":%s,\"pipeline\":%s,\"sample\":%s}\n",
                   quote(active).c_str(), quote(e.pipeline).c_str(), sample.str().c_str());
  }
}
void verifyMLX() {
  Dl_info info;
  if (!dladdr((void*)&mx::metal::get_command_encoder, &info)) std::abort();
  FILE* file = std::fopen(info.dli_fname, "rb"); if (!file) std::abort();
  CC_SHA256_CTX context; CC_SHA256_Init(&context);
  unsigned char buffer[65536], digest[32]; size_t n;
  while ((n = std::fread(buffer, 1, sizeof buffer, file))) CC_SHA256_Update(&context, buffer, (CC_LONG)n);
  std::fclose(file); CC_SHA256_Final(digest, &context);
  std::string hex;
  for (auto c : digest) { hex += "0123456789abcdef"[c >> 4]; hex += "0123456789abcdef"[c & 15]; }
  if (hex != "72dc8d9a4b1e07df2f6bc91383e0567acf3ce292398ae6ce4ed1d19fce0b1374") {
    std::fprintf(stderr, "Refusing unrecognized MLX ABI: %s\n", hex.c_str()); std::abort();
  }
}
}
// Read the SAME image that supplies the traced CommandEncoder. Loading an
// installed copy by pathname can create a second image with independent globals,
// even when its bytes and install name are identical. Prefix/count checks below
// additionally prove that this encoder is the one the plugin actually uses.
extern "C" void fbc_trace_stats(mlx::core::metal::DispatchStats* out) {
  *out = mlx::core::metal::dispatch_stats();
}
extern "C" const char* fbc_trace_mlx_image() {
  Dl_info info;
  if (!dladdr((void*)&mlx::core::metal::get_command_encoder, &info)) return nullptr;
  return info.dli_fname;
}
// Unlike MLX's completed counter, this also waits for OUR logging callback.
// This is diagnostic-only host waiting and must never be used for timings.
extern "C" uint64_t fbc_trace_wait() {
  for (int i = 0; i < 2000; ++i) {
    {
      std::lock_guard<std::mutex> lock(state().mutex);
      if (state().pending == 0) return 0;
    }
    std::this_thread::sleep_for(std::chrono::milliseconds(1));
  }
  std::lock_guard<std::mutex> lock(state().mutex);
  return state().pending;
}
extern "C" void fbc_trace_dump() {
  std::lock_guard<std::mutex> lock(state().mutex);
  if (!state().output) return;
  for (const auto& [key, k] : state().kernels) {
    std::fprintf(state().output, "{\"type\":\"kernel\",\"phase\":%s,\"pipeline\":%s,\"count\":%llu,\"cpu_dispatch_ns\":%llu,\"binding_span_bytes_upper_bound\":%llu,\"arguments\":%s,\"samples\":[",
                 quote(key.first).c_str(), quote(key.second).c_str(), k.count, k.cpu_dispatch_ns, k.binding_span_bytes, k.arguments.c_str());
    bool comma = false;
    for (const auto& sample : k.samples) { if (comma) std::fputc(',', state().output); comma = true; std::fputs(sample.c_str(), state().output); }
    std::fputs("]}\n", state().output);
  }
  std::fflush(state().output);
}
__attribute__((constructor)) static void install() {
  const char* path = std::getenv("FBC_MLX_TRACE_OUTPUT"); if (!path) return;
  verifyMLX(); state().output = std::fopen(path, "wx"); if (!state().output) std::abort();
  @autoreleasepool {
    id<MTLDevice> device = MTLCreateSystemDefaultDevice();
    hook(object_getClass(device), @selector(newComputePipelineStateWithFunction:error:), (IMP)newPipeline);
    hook(object_getClass(device), @selector(newComputePipelineStateWithDescriptor:options:reflection:error:), (IMP)newPipelineDescriptor);
    id<MTLCommandQueue> queue = [device newCommandQueue];
    id<MTLCommandBuffer> buffer = [queue commandBuffer];
    for (MTLDispatchType type : {MTLDispatchTypeSerial, MTLDispatchTypeConcurrent}) {
      id<MTLComputeCommandEncoder> encoder = [buffer computeCommandEncoderWithDispatchType:type];
      Class cls = object_getClass(encoder);
      hook(cls, @selector(setComputePipelineState:), (IMP)pipeline);
      hook(cls, @selector(setBuffer:offset:atIndex:), (IMP)binding);
      hook(cls, @selector(setBytes:length:atIndex:), (IMP)constants);
      hook(cls, @selector(dispatchThreads:threadsPerThreadgroup:), (IMP)dispatch);
      hook(cls, @selector(dispatchThreadgroups:threadsPerThreadgroup:), (IMP)dispatch);
      [encoder endEncoding];
    }
    hook(object_getClass(buffer), @selector(commit), (IMP)commit);
    hook(object_getClass(buffer), @selector(computeCommandEncoderWithDispatchType:), (IMP)makeEncoder);
    [queue release]; [device release];
  }
}
