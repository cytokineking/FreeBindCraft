# Fix for "Too Many Open Files" Error in BindCraft

## Problem Description

Users were encountering the error:
```
OSError: [Errno 24] Too many open files
```

This error occurs in the DSSP subprocess call within `calc_ss_percentage()` function, which is called frequently during the design loop (once per trajectory + once per MPNN design).

## Root Cause

The issue was caused by:
1. **High frequency of DSSP calls**: `calc_ss_percentage()` is called potentially hundreds of times per run
2. **Subprocess file descriptor leaks**: The BioPython DSSP class creates subprocesses that may not properly close file descriptors
3. **Insufficient garbage collection**: Python's garbage collector wasn't running frequently enough to clean up resources

## Solutions Implemented

### 1. Enhanced DSSP Error Handling and Retry Logic
- **File**: `functions/biopython_utils.py`
- **Function**: `safe_dssp_calculation()`
- Added retry logic with proper cleanup between attempts
- Added graceful fallback to default values if DSSP fails completely
- Improved exception handling with explicit cleanup

### 2. DSSP Result Caching
- **File**: `functions/biopython_utils.py`
- **Feature**: `_dssp_cache` global dictionary
- Caches DSSP results to avoid redundant calculations on the same PDB files
- Includes cache clearing function to prevent memory bloat
- Reduces the total number of DSSP subprocess calls significantly

### 3. Aggressive Garbage Collection
- **File**: `bindcraft.py`
- Added `gc.collect()` calls after each MPNN design iteration
- Added periodic DSSP cache clearing every 10 trajectories
- Enhanced memory cleanup throughout the design loop

### 4. User-Friendly Ulimit Helper Script
- **File**: `check_ulimit.sh`
- Checks current file descriptor limits
- Provides instructions for temporary and permanent fixes
- Offers interactive limit adjustment

## Files Modified

1. **`functions/biopython_utils.py`**:
   - Added `safe_dssp_calculation()` function
   - Added DSSP caching mechanism
   - Enhanced error handling in `calc_ss_percentage()`
   - Added explicit garbage collection

2. **`bindcraft.py`**:
   - Added `gc` import
   - Added frequent garbage collection calls
   - Added periodic DSSP cache clearing
   - Enhanced resource cleanup

3. **`check_ulimit.sh`** (new file):
   - Interactive script to check and fix file descriptor limits

## Usage Instructions

### Immediate Fix
1. Run the ulimit helper script:
   ```bash
   ./check_ulimit.sh
   ```

2. Or manually increase the limit:
   ```bash
   ulimit -n 65536
   ```

### Long-term Fix
The code changes should prevent the issue from occurring, but you can also make the ulimit change permanent by adding to your shell profile:
```bash
echo 'ulimit -n 65536' >> ~/.bashrc
source ~/.bashrc
```

## Technical Details

### DSSP Caching Strategy
- Cache key: PDB file path
- Cache stores both successful DSSP objects and failures (None)
- Cache is cleared every 10 trajectories to prevent memory bloat
- Failed calculations are cached to avoid repeated attempts

### Garbage Collection Strategy
- Called after each MPNN design iteration
- Called at the end of each trajectory
- DSSP cache clearing every 10 trajectories
- Explicit cleanup of large objects (models, structures)

### Error Recovery
- DSSP failures now return sensible default values instead of crashing
- Multiple retry attempts with cleanup between tries
- Graceful degradation when DSSP is unavailable

## Testing

To verify the fix is working:
1. Monitor file descriptor usage: `lsof -p <bindcraft_pid> | wc -l`
2. Check for DSSP error messages in the output
3. Observe cache clearing messages every 10 trajectories
4. Run longer design sessions that previously failed

## Performance Impact

- **Positive**: Reduced DSSP calls due to caching
- **Minimal**: Small overhead from garbage collection
- **Positive**: Better memory management prevents slowdowns
- **Positive**: Fewer subprocess creation/destruction cycles

## Future Improvements

Consider implementing:
1. More sophisticated cache eviction policies (LRU)
2. Persistent DSSP caching across runs
3. Alternative secondary structure prediction methods
4. Monitoring and alerting for resource usage
