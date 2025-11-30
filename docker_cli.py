#!/usr/bin/env python3

"""
docker_cli.py

Standalone interactive wrapper for running BindCraft in Docker.
Mirrors the interactive prompts and behavior in bindcraft.py without importing heavy deps.

Requirements: Docker installed on host and accessible from this process.
"""

import os
import sys
import json
import time
import subprocess


def input_with_default(prompt_text: str, default_value: str | None = None) -> str:
    if default_value is None:
        return input(prompt_text).strip()
    resp = input(f"{prompt_text} ").strip()
    return resp if resp else default_value


def yes_no(prompt_text: str, default_yes: bool = False) -> bool:
    default_hint = 'Y/n' if default_yes else 'y/N'
    resp = input(f"{prompt_text} ({default_hint}): ").strip().lower()
    if resp == '':
        return default_yes
    return resp in ('y', 'yes')


def require_executable(name: str) -> None:
    try:
        subprocess.run([name, "--version"], capture_output=True, text=True)
    except Exception:
        print(f"Error: '{name}' not found. Please install Docker and try again.")
        sys.exit(1)


def list_local_docker_images() -> list[str]:
    try:
        proc = subprocess.run(["docker", "images", "--format", "{{.Repository}}:{{.Tag}}"], capture_output=True, text=True)
        images = []
        for ln in proc.stdout.strip().splitlines():
            val = ln.strip()
            if not val or val == "<none>:<none>":
                continue
            images.append(val)
        # Deduplicate preserving order
        seen = set()
        return [x for x in images if not (x in seen or seen.add(x))]
    except Exception:
        return []


def list_json_basenames_in_image(image: str, container_path: str) -> list[str]:
    """Return json basenames (without .json) under container_path in the image."""
    # Use python in the container for robust listing
    code = (
        "import os,glob;"
        f"print('\n'.join(sorted(os.path.splitext(os.path.basename(p))[0] for p in glob.glob('{container_path}/*.json'))))"
    )
    try:
        proc = subprocess.run(
            ["docker", "run", "--rm", image, "python", "-c", code],
            capture_output=True, text=True
        )
        if proc.returncode != 0:
            return []
        names = [ln.strip() for ln in proc.stdout.splitlines() if ln.strip()]
        return names
    except Exception:
        return []


def main() -> None:
    require_executable("docker")

    print("\nBindCraft Docker Interactive Setup\n")

    while True:
        # Choose Docker image
        images = list_local_docker_images()
        selected_image = None
        if images:
            print("Local Docker images:")
            for idx, name in enumerate(images, 1):
                print(f"{idx}. {name}")
            print("0. Enter image name manually")
            choice = input_with_default("Choose image (press Enter for freebindcraft:gpu):", "")
            if not choice:
                selected_image = "freebindcraft:gpu"
            else:
                try:
                    cidx = int(choice)
                    if cidx == 0:
                        selected_image = None
                    elif 1 <= cidx <= len(images):
                        selected_image = images[cidx - 1]
                except Exception:
                    selected_image = None
        if not selected_image:
            selected_image = input_with_default("Docker image to use (default freebindcraft:gpu):", "freebindcraft:gpu")

        # GPU index
        gpu_idx = input_with_default("GPU index to expose (default 0):", "0")

        # Design type
        print("Design type:")
        print("1. Miniprotein (31+ aa)")
        print("2. Peptide (8-30 aa)")
        dtype_choice = input_with_default("Choose design type (press Enter for Miniprotein):", "")
        is_peptide = (dtype_choice.strip() == '2')

        # Required inputs
        project_name = input_with_default("Enter project/binder name:", None)
        while not project_name:
            print("Project name is required.")
            project_name = input_with_default("Enter project/binder name:", None)

        # PDB path (must exist)
        while True:
            pdb_raw = input_with_default("Enter path to PDB file:", None)
            if not pdb_raw:
                print("PDB path is required.")
                continue
            candidate = os.path.abspath(os.path.expanduser(pdb_raw))
            if os.path.isfile(candidate):
                pdb_path = candidate
                break
            print(f"Error: No PDB file found at '{candidate}'. Please re-enter.")

        # Output directory (auto-create)
        output_dir = input_with_default("Enter output directory:", os.path.join(os.getcwd(), f"{project_name}_bindcraft_out"))
        output_dir = os.path.abspath(os.path.expanduser(output_dir))
        os.makedirs(output_dir, exist_ok=True)

        chains = input_with_default("Enter target chains (e.g., A or A,B):", "A")
        hotspots = input_with_default("Enter hotspot residue(s) for BindCraft to target. Use format: chain letter + residue numbers (e.g., 'A1,B20-25'). Leave empty for no preference:", "")

        # Lengths
        if is_peptide:
            lengths_prompt_default = "8 25"
            lengths_input = input_with_default("Enter peptide min and max lengths (8-30) separated by space or comma (default 8 25):", lengths_prompt_default)
        else:
            lengths_prompt_default = "65 150"
            lengths_input = input_with_default("Enter miniprotein min and max lengths separated by space or comma (min>=31, default 65 150):", lengths_prompt_default)
        try:
            normalized = lengths_input.replace(',', ' ').split()
            min_len_val, max_len_val = int(normalized[0]), int(normalized[1])
        except Exception:
            min_len_val, max_len_val = (8, 25) if is_peptide else (65, 150)
        if is_peptide:
            min_len_val = max(8, min(min_len_val, 30))
            max_len_val = max(8, min(max_len_val, 30))
            if min_len_val > max_len_val:
                min_len_val, max_len_val = max_len_val, min_len_val
        else:
            if min_len_val < 31:
                min_len_val = 31
            if max_len_val < min_len_val:
                max_len_val = min_len_val
        lengths = [min_len_val, max_len_val]

        # Number of final designs
        num_designs_str = input_with_default("Enter number of final designs (default 100):", "100")
        try:
            num_designs = int(num_designs_str)
        except Exception:
            num_designs = 100

        # Fetch filters/advanced lists from image
        filter_names = list_json_basenames_in_image(selected_image, "/app/settings_filters")
        adv_names = list_json_basenames_in_image(selected_image, "/app/settings_advanced")

        # Filter menu per design type
        if is_peptide:
            filter_order = ['peptide_filters', 'peptide_relaxed_filters', 'no_filters']
            default_filter_name = 'peptide_filters'
        else:
            filter_order = ['default_filters', 'relaxed_filters', 'no_filters']
            default_filter_name = 'default_filters'
        ordered_filters = [n for n in filter_order if n in filter_names]
        if not ordered_filters:
            # Fallback in case image list failed
            ordered_filters = filter_order
        print("\nAvailable filter settings:")
        for i, name in enumerate(ordered_filters, 1):
            print(f"{i}. {name}")
        filter_idx = input_with_default(f"Choose filter (press Enter for {default_filter_name}):", "")
        if filter_idx:
            try:
                filter_idx_int = int(filter_idx)
                selected_filter_name = ordered_filters[filter_idx_int - 1]
            except Exception:
                selected_filter_name = default_filter_name
        else:
            selected_filter_name = default_filter_name

        # Advanced menu per design type
        if is_peptide:
            adv_order = [
                'peptide_3stage_multimer',
                'peptide_3stage_multimer_mpnn',
                'peptide_3stage_multimer_flexible',
                'peptide_3stage_multimer_mpnn_flexible'
            ]
            default_adv_name = 'peptide_3stage_multimer'
        else:
            adv_order = [
                'default_4stage_multimer',
                'default_4stage_multimer_mpnn',
                'default_4stage_multimer_flexible',
                'default_4stage_multimer_hardtarget',
                'default_4stage_multimer_flexible_hardtarget',
                'default_4stage_multimer_mpnn_flexible',
                'default_4stage_multimer_mpnn_hardtarget',
                'default_4stage_multimer_mpnn_flexible_hardtarget',
                'betasheet_4stage_multimer',
                'betasheet_4stage_multimer_mpnn',
                'betasheet_4stage_multimer_flexible',
                'betasheet_4stage_multimer_hardtarget',
                'betasheet_4stage_multimer_flexible_hardtarget',
                'betasheet_4stage_multimer_mpnn_flexible',
                'betasheet_4stage_multimer_mpnn_hardtarget',
                'betasheet_4stage_multimer_mpnn_flexible_hardtarget'
            ]
            default_adv_name = 'default_4stage_multimer'
        ordered_adv = [n for n in adv_order if n in adv_names]
        if not ordered_adv:
            ordered_adv = adv_order
        print("\nAvailable advanced settings:")
        for i, name in enumerate(ordered_adv, 1):
            print(f"{i}. {name}")
        advanced_idx = input_with_default(f"Choose advanced (press Enter for {default_adv_name}):", "")
        if advanced_idx:
            try:
                advanced_idx_int = int(advanced_idx)
                selected_adv_name = ordered_adv[advanced_idx_int - 1]
            except Exception:
                selected_adv_name = default_adv_name
        else:
            selected_adv_name = default_adv_name

        # Toggles
        verbose = yes_no("Enable verbose output?", default_yes=False)
        plots_on = yes_no("Enable saving plots?", default_yes=True)
        animations_on = yes_no("Enable saving animations?", default_yes=True)
        run_with_pyrosetta = yes_no("Run with PyRosetta?", default_yes=True)

        # Ranking method selection
        print("\nRanking method for final designs:")
        print("1. i_pTM (interface predicted TM-score)")
        print("2. ipSAE (interface predicted Structural Alignment Error)")
        rank_choice = input_with_default("Choose ranking method (press Enter for i_pTM):", "")
        if rank_choice.strip() == '2':
            rank_by_metric = 'ipSAE'
        else:
            rank_by_metric = 'i_pTM'

        # Summary
        print("\nConfiguration Summary:")
        print(f"Docker Image: {selected_image}")
        print(f"GPU Index: {gpu_idx}")
        print(f"Project Name: {project_name}")
        print(f"PDB File: {pdb_path}")
        print(f"Output Directory: {output_dir}")
        print(f"Design Type: {'Peptide' if is_peptide else 'Miniprotein'}")
        print(f"Chains: {chains}")
        print(f"Hotspots: {hotspots if hotspots else 'None'}")
        print(f"Length Range: {lengths}")
        print(f"Number of Final Designs: {num_designs}")
        print(f"Filter Setting: {selected_filter_name}")
        print(f"Advanced Setting: {selected_adv_name}")
        print(f"Verbose: {'Yes' if verbose else 'No'}")
        print(f"Plots: {'On' if plots_on else 'Off'}")
        print(f"Animations: {'On' if animations_on else 'Off'}")
        print(f"PyRosetta: {'On' if run_with_pyrosetta else 'Off'}")
        print(f"Ranking Method: {rank_by_metric}")

        if yes_no("Proceed with these settings?", default_yes=True):
            # Materialize target settings JSON
            target_settings = {
                "design_path": output_dir,
                "binder_name": project_name,
                "starting_pdb": pdb_path,
                "chains": chains,
                "target_hotspot_residues": hotspots,
                "lengths": lengths,
                "number_of_final_designs": num_designs
            }
            timestamp = time.strftime('%Y%m%d_%H%M%S')
            settings_filename = f"interactive_{project_name}_{timestamp}.json"
            settings_path_out = os.path.join(output_dir, settings_filename)
            try:
                with open(settings_path_out, 'w') as f:
                    json.dump(target_settings, f, indent=2)
            except Exception as e:
                print(f"Error writing settings JSON: {e}")
                sys.exit(1)

            # Compose container paths (mount host paths into identical locations)
            mounts = []
            cwd = os.getcwd()
            mounts.append((cwd, cwd))
            pdb_parent = os.path.dirname(pdb_path)
            if pdb_parent and os.path.abspath(pdb_parent) != cwd:
                mounts.append((pdb_parent, pdb_parent))
            mounts.append((output_dir, output_dir))

            # Assemble docker run command
            cmd = [
                "docker", "run", "--rm", "-it",
                "--gpus", f"device={gpu_idx}",
            ]
            for host, cont in mounts:
                cmd.extend(["-v", f"{host}:{cont}"])
            cmd.extend(["-w", cwd, selected_image, "python", "bindcraft.py",
                        "-s", settings_path_out,
                        "-f", f"/app/settings_filters/{selected_filter_name}.json",
                        "-a", f"/app/settings_advanced/{selected_adv_name}.json"])
            if not run_with_pyrosetta:
                cmd.append("--no-pyrosetta")
            if verbose:
                cmd.append("--verbose")
            if not plots_on:
                cmd.append("--no-plots")
            if not animations_on:
                cmd.append("--no-animations")
            if rank_by_metric != 'i_pTM':
                cmd.extend(["--rank-by", rank_by_metric])

            print("\nLaunching Docker:")
            print(" ".join(cmd))
            # Execute
            subprocess.run(cmd)
            return
        else:
            print("Let's re-enter the details.\n")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nAborted by user.")
        sys.exit(1)


