#!/usr/bin/env python3
"""
Physics validation test runner for metalquicha
Runs calculations and validates total energies against expected values
"""

import json
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional


class Colors:
    """ANSI color codes for terminal output"""
    GREEN = '\033[92m'
    RED = '\033[91m'
    YELLOW = '\033[93m'
    BLUE = '\033[94m'
    RESET = '\033[0m'
    BOLD = '\033[1m'


def run_calculation(input_file: str, exe_path: str = "./build/mqc", verbose: bool = False,
                   use_mpi: bool = False, nprocs: int = 4, test_type: str = "unfragmented") -> bool:
    """Run a metalquicha calculation and cache stdout/stderr"""
    import os

    # Set environment for reproducible calculations
    env = os.environ.copy()
    env['OMP_NUM_THREADS'] = '1'

    # Determine log filename from input file
    input_path = Path(input_file)
    log_file = f"validation_logs/{input_path.stem}.log"

    # Create logs directory if needed
    Path("validation_logs").mkdir(exist_ok=True)

    # Build command based on MPI settings
    if use_mpi:
        # Use 1 process for unfragmented, nprocs for fragmented
        np = 1 if test_type == "unfragmented" else nprocs
        cmd = ["mpirun", "-np", str(np), exe_path, input_file]
    else:
        cmd = [exe_path, input_file]

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300,  # 5 minute timeout
            env=env
        )

        # Save stdout and stderr to log file
        with open(log_file, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write(f"Calculation: {input_file}\n")
            f.write(f"Command: {' '.join(cmd)}\n")
            f.write("=" * 80 + "\n\n")
            f.write("STDOUT:\n")
            f.write("-" * 80 + "\n")
            f.write(result.stdout)
            f.write("\n" + "-" * 80 + "\n\n")
            f.write("STDERR:\n")
            f.write("-" * 80 + "\n")
            f.write(result.stderr)
            f.write("\n" + "-" * 80 + "\n\n")
            f.write(f"Return code: {result.returncode}\n")

        if verbose:
            print(f"  Log saved to: {log_file}")

        return result.returncode == 0
    except subprocess.TimeoutExpired:
        print(f"  {Colors.RED}✗ Timeout running {input_file}{Colors.RESET}")
        # Save timeout info to log
        with open(log_file, 'w') as f:
            f.write(f"TIMEOUT after 300 seconds\n")
            f.write(f"Calculation: {input_file}\n")
        return False
    except Exception as e:
        print(f"  {Colors.RED}✗ Error running {input_file}: {e}{Colors.RESET}")
        # Save error info to log
        with open(log_file, 'w') as f:
            f.write(f"ERROR: {e}\n")
            f.write(f"Calculation: {input_file}\n")
        return False


def get_output_filename(input_file: str) -> str:
    """Get the output JSON filename from the input filename"""
    basename = Path(input_file).stem  # Remove .mqc extension
    return f"output_{basename}.json"


def extract_top_level_energy(json_data: Dict, test_type: str) -> Optional[float]:
    """Extract the top-level total_energy from JSON output"""

    # The structure is: { "basename": { "total_energy": value, ... } }
    # or for multi-molecule: { "basename": { "mol1": {...}, "mol2": {...} } }

    # Get the first (and usually only) top-level key
    if not json_data:
        return None

    top_key = list(json_data.keys())[0]
    data = json_data[top_key]

    # For unfragmented or fragmented (not multi-molecule), energy is directly under top key
    if "total_energy" in data:
        return float(data["total_energy"])

    return None


def extract_gradient_norm(json_data: Dict) -> Optional[float]:
    """Extract gradient_norm from JSON output"""
    if not json_data:
        return None

    top_key = list(json_data.keys())[0]
    data = json_data[top_key]

    if "gradient_norm" in data:
        return float(data["gradient_norm"])

    return None


def extract_hessian_norm(json_data: Dict) -> Optional[float]:
    """Extract hessian_frobenius_norm from JSON output"""
    if not json_data:
        return None

    top_key = list(json_data.keys())[0]
    data = json_data[top_key]

    if "hessian_frobenius_norm" in data:
        return float(data["hessian_frobenius_norm"])

    return None


def extract_frequencies(json_data: Dict) -> Optional[List[float]]:
    """Extract vibrational frequencies from JSON output"""
    if not json_data:
        return None

    top_key = list(json_data.keys())[0]
    data = json_data[top_key]

    if "vibrational_analysis" in data:
        vib_data = data["vibrational_analysis"]
        if "frequencies_cm1" in vib_data:
            return [float(f) for f in vib_data["frequencies_cm1"]]

    return None


def extract_zpe(json_data: Dict) -> Optional[float]:
    """Extract zero-point energy in Hartree from JSON output"""
    if not json_data:
        return None

    top_key = list(json_data.keys())[0]
    data = json_data[top_key]

    if "thermochemistry" in data:
        thermo = data["thermochemistry"]
        if "zero_point_energy_hartree" in thermo:
            return float(thermo["zero_point_energy_hartree"])

    return None


def extract_gibbs_correction(json_data: Dict) -> Optional[float]:
    """Extract thermal correction to Gibbs free energy from JSON output"""
    if not json_data:
        return None

    top_key = list(json_data.keys())[0]
    data = json_data[top_key]

    if "thermochemistry" in data:
        thermo = data["thermochemistry"]
        if "thermal_corrections_hartree" in thermo:
            corrections = thermo["thermal_corrections_hartree"]
            if "to_gibbs" in corrections:
                return float(corrections["to_gibbs"])

    return None


def extract_multi_molecule_energies(json_data: Dict) -> Dict[str, float]:
    """Extract energies from multi-molecule calculation"""
    energies = {}

    if not json_data:
        return energies

    top_key = list(json_data.keys())[0]
    molecules = json_data[top_key]

    # Iterate through each molecule
    for mol_name, mol_data in molecules.items():
        if isinstance(mol_data, dict) and "total_energy" in mol_data:
            energies[mol_name] = float(mol_data["total_energy"])

    return energies


def validate_energy(calculated: float, expected: float, tolerance: float) -> bool:
    """Check if calculated energy matches expected within tolerance"""
    return abs(calculated - expected) < tolerance


def cleanup_output_files(verbose: bool = False) -> None:
    """Remove any existing output_*.json files to ensure clean test runs"""
    output_files = list(Path(".").glob("output_*.json"))

    if not output_files:
        if verbose:
            print("No output_*.json files to clean up.\n")
        return

    print(f"{Colors.BOLD}Cleaning up {len(output_files)} existing output file(s)...{Colors.RESET}")

    for output_file in output_files:
        try:
            output_file.unlink()
            if verbose:
                print(f"  Removed: {output_file}")
        except Exception as e:
            print(f"  {Colors.YELLOW}Warning: Could not remove {output_file}: {e}{Colors.RESET}")

    print()


def prepare_mqc_files(validation_dir: str = "validation", prep_script: str = "mqc_prep.py", verbose: bool = False) -> None:
    """Run mqc_prep.py on all JSON files in validation/inputs/"""
    inputs_dir = Path(validation_dir) / "inputs"

    if not inputs_dir.exists():
        print(f"{Colors.YELLOW}Warning: {inputs_dir} does not exist{Colors.RESET}")
        return

    json_files = list(inputs_dir.glob("*.json"))

    if not json_files:
        print(f"{Colors.YELLOW}No JSON files found in {inputs_dir}{Colors.RESET}")
        return

    print(f"\n{Colors.BOLD}Preparing .mqc files from {len(json_files)} JSON inputs...{Colors.RESET}\n")

    for json_file in json_files:
        output_mqc = inputs_dir / f"{json_file.stem}.mqc"

        if verbose:
            print(f"  {json_file.name} -> {output_mqc.name}")

        try:
            result = subprocess.run(
                ["python3", prep_script, str(json_file)],
                capture_output=True,
                text=True,
                timeout=30
            )
            if result.returncode != 0:
                print(f"  {Colors.RED}✗ Failed to prepare {json_file.name}{Colors.RESET}")
                if verbose:
                    print(f"    {result.stderr}")
        except Exception as e:
            print(f"  {Colors.RED}✗ Error preparing {json_file.name}: {e}{Colors.RESET}")

    print()


def run_validation_tests(manifest_file: str = "validation_tests.json",
                        exe_path: str = "./build/mqc",
                        verbose: bool = False,
                        use_mpi: bool = False,
                        nprocs: int = 4,
                        test_filter: Optional[str] = None) -> tuple:
    """Run all validation tests and return (passed, failed) counts"""

    # Load validation manifest
    with open(manifest_file, 'r') as f:
        manifest = json.load(f)

    tolerance = manifest.get("tolerance", 1.0e-9)
    tests = manifest.get("tests", [])

    # Filter tests if requested
    if test_filter:
        tests = [t for t in tests if test_filter.lower() in t["name"].lower() or
                 test_filter.lower() in t["input"].lower()]
        if not tests:
            print(f"{Colors.YELLOW}No tests matching '{test_filter}'{Colors.RESET}")
            return 0, 0, []

    passed = 0
    failed = 0
    errors = []

    print(f"\n{Colors.BOLD}Running {len(tests)} validation tests...{Colors.RESET}")
    print(f"Tolerance: {tolerance}")
    if use_mpi:
        print(f"MPI mode: enabled (nprocs={nprocs} for fragmented, np=1 for unfragmented)")
    print()

    for i, test in enumerate(tests, 1):
        test_name = test["name"]
        input_file = test["input"]
        test_type = test.get("type", "unfragmented")

        print(f"{i}/{len(tests)}: {Colors.BLUE}{test_name}{Colors.RESET} [{test_type}]")

        # Run calculation
        if verbose:
            if use_mpi:
                np = 1 if test_type == "unfragmented" else nprocs
                print(f"  Running: mpirun -np {np} {exe_path} {input_file}")
            else:
                print(f"  Running: {exe_path} {input_file}")

        if not run_calculation(input_file, exe_path, verbose=verbose,
                             use_mpi=use_mpi, nprocs=nprocs, test_type=test_type):
            print(f"  {Colors.RED}✗ FAILED{Colors.RESET} - Calculation error")
            print(f"  See validation_logs/{Path(input_file).stem}.log for details\n")
            failed += 1
            errors.append((test_name, "Calculation failed"))
            continue

        # Read output JSON
        output_file = get_output_filename(input_file)
        try:
            with open(output_file, 'r') as f:
                output_data = json.load(f)
        except FileNotFoundError:
            print(f"  {Colors.RED}✗ FAILED{Colors.RESET} - Output file not found: {output_file}\n")
            failed += 1
            errors.append((test_name, f"Output file not found: {output_file}"))
            continue
        except json.JSONDecodeError as e:
            print(f"  {Colors.RED}✗ FAILED{Colors.RESET} - Invalid JSON: {e}\n")
            failed += 1
            errors.append((test_name, f"Invalid JSON: {e}"))
            continue

        # Validate energy
        if test_type == "multi_molecule":
            # Handle multi-molecule separately
            expected_energies = test["expected_energies"]
            calculated_energies = extract_multi_molecule_energies(output_data)

            all_match = True
            for mol_name, expected in expected_energies.items():
                if mol_name not in calculated_energies:
                    print(f"  {Colors.RED}✗ FAILED{Colors.RESET} - Missing molecule: {mol_name}\n")
                    all_match = False
                    break

                calculated = calculated_energies[mol_name]
                diff = abs(calculated - expected)

                if not validate_energy(calculated, expected, tolerance):
                    print(f"  {Colors.RED}✗ FAILED{Colors.RESET} - {mol_name}")
                    print(f"    Expected:   {expected:.12f}")
                    print(f"    Calculated: {calculated:.12f}")
                    print(f"    Difference: {diff:.2e} (tolerance: {tolerance:.2e})\n")
                    all_match = False
                    break
                elif verbose:
                    print(f"    {mol_name}: {calculated:.12f} (diff: {diff:.2e})")

            if all_match:
                print(f"  {Colors.GREEN}✓ PASSED{Colors.RESET}\n")
                passed += 1
            else:
                failed += 1
                errors.append((test_name, "Energy mismatch"))

        else:
            # Single molecule (unfragmented or fragmented)
            expected = test["expected_energy"]
            calculated = extract_top_level_energy(output_data, test_type)

            if calculated is None:
                print(f"  {Colors.RED}✗ FAILED{Colors.RESET} - Could not extract energy from JSON\n")
                failed += 1
                errors.append((test_name, "Could not extract energy"))
                continue

            diff = abs(calculated - expected)
            test_passed = True
            failure_reasons = []

            # Validate energy
            if not validate_energy(calculated, expected, tolerance):
                print(f"  {Colors.RED}✗ FAILED{Colors.RESET} - Energy mismatch")
                print(f"    Expected:   {expected:.12f}")
                print(f"    Calculated: {calculated:.12f}")
                print(f"    Difference: {diff:.2e} (tolerance: {tolerance:.2e})")
                test_passed = False
                failure_reasons.append(f"Energy mismatch (diff: {diff:.2e})")

            # Validate gradient norm if present
            if "expected_gradient_norm" in test:
                expected_grad = test["expected_gradient_norm"]
                calculated_grad = extract_gradient_norm(output_data)

                if calculated_grad is None:
                    print(f"  {Colors.RED}✗ FAILED{Colors.RESET} - Could not extract gradient_norm from JSON")
                    test_passed = False
                    failure_reasons.append("Missing gradient_norm")
                else:
                    grad_diff = abs(calculated_grad - expected_grad)
                    if not validate_energy(calculated_grad, expected_grad, tolerance):
                        print(f"  {Colors.RED}✗ FAILED{Colors.RESET} - Gradient norm mismatch")
                        print(f"    Expected:   {expected_grad:.12f}")
                        print(f"    Calculated: {calculated_grad:.12f}")
                        print(f"    Difference: {grad_diff:.2e} (tolerance: {tolerance:.2e})")
                        test_passed = False
                        failure_reasons.append(f"Gradient mismatch (diff: {grad_diff:.2e})")
                    elif verbose:
                        print(f"    Gradient norm: {calculated_grad:.12f} (diff: {grad_diff:.2e})")

            # Validate Hessian Frobenius norm if present
            if "expected_hessian_frobenius_norm" in test:
                expected_hess = test["expected_hessian_frobenius_norm"]
                calculated_hess = extract_hessian_norm(output_data)

                if calculated_hess is None:
                    print(f"  {Colors.RED}✗ FAILED{Colors.RESET} - Could not extract hessian_frobenius_norm from JSON")
                    test_passed = False
                    failure_reasons.append("Missing hessian_frobenius_norm")
                else:
                    hess_diff = abs(calculated_hess - expected_hess)
                    if not validate_energy(calculated_hess, expected_hess, tolerance):
                        print(f"  {Colors.RED}✗ FAILED{Colors.RESET} - Hessian norm mismatch")
                        print(f"    Expected:   {expected_hess:.12f}")
                        print(f"    Calculated: {calculated_hess:.12f}")
                        print(f"    Difference: {hess_diff:.2e} (tolerance: {tolerance:.2e})")
                        test_passed = False
                        failure_reasons.append(f"Hessian mismatch (diff: {hess_diff:.2e})")
                    elif verbose:
                        print(f"    Hessian norm: {calculated_hess:.12f} (diff: {hess_diff:.2e})")

            # Validate frequencies if present
            if "expected_frequencies" in test:
                expected_freqs = test["expected_frequencies"]
                calculated_freqs = extract_frequencies(output_data)

                if calculated_freqs is None:
                    print(f"  {Colors.RED}✗ FAILED{Colors.RESET} - Could not extract frequencies from JSON")
                    test_passed = False
                    failure_reasons.append("Missing frequencies")
                elif len(calculated_freqs) != len(expected_freqs):
                    print(f"  {Colors.RED}✗ FAILED{Colors.RESET} - Frequency count mismatch")
                    print(f"    Expected:   {len(expected_freqs)} frequencies")
                    print(f"    Calculated: {len(calculated_freqs)} frequencies")
                    test_passed = False
                    failure_reasons.append(f"Frequency count mismatch ({len(calculated_freqs)} vs {len(expected_freqs)})")
                else:
                    freq_errors = []
                    # Use looser tolerance for near-zero frequencies (trans/rot modes)
                    # These are numerically noisy and vary between machines
                    near_zero_threshold = 10.0  # cm^-1
                    near_zero_tolerance = 1.0   # cm^-1 tolerance for near-zero modes
                    for idx, (calc_f, exp_f) in enumerate(zip(calculated_freqs, expected_freqs)):
                        freq_diff = abs(calc_f - exp_f)
                        # Use looser tolerance if both values are near zero
                        if abs(exp_f) < near_zero_threshold and abs(calc_f) < near_zero_threshold:
                            effective_tolerance = near_zero_tolerance
                        else:
                            effective_tolerance = tolerance
                        if freq_diff >= effective_tolerance:
                            freq_errors.append((idx, calc_f, exp_f, freq_diff))

                    if freq_errors:
                        print(f"  {Colors.RED}✗ FAILED{Colors.RESET} - Frequency mismatch ({len(freq_errors)} modes)")
                        for idx, calc_f, exp_f, freq_diff in freq_errors[:3]:  # Show first 3
                            print(f"    Mode {idx}: expected {exp_f:.4f}, got {calc_f:.4f} (diff: {freq_diff:.2e})")
                        if len(freq_errors) > 3:
                            print(f"    ... and {len(freq_errors) - 3} more")
                        test_passed = False
                        failure_reasons.append(f"Frequency mismatch ({len(freq_errors)} modes)")
                    elif verbose:
                        print(f"    Frequencies: {len(calculated_freqs)} modes validated")

            # Validate ZPE if present
            if "expected_zpe" in test:
                expected_zpe = test["expected_zpe"]
                calculated_zpe = extract_zpe(output_data)

                if calculated_zpe is None:
                    print(f"  {Colors.RED}✗ FAILED{Colors.RESET} - Could not extract ZPE from JSON")
                    test_passed = False
                    failure_reasons.append("Missing ZPE")
                else:
                    zpe_diff = abs(calculated_zpe - expected_zpe)
                    if not validate_energy(calculated_zpe, expected_zpe, tolerance):
                        print(f"  {Colors.RED}✗ FAILED{Colors.RESET} - ZPE mismatch")
                        print(f"    Expected:   {expected_zpe:.12f}")
                        print(f"    Calculated: {calculated_zpe:.12f}")
                        print(f"    Difference: {zpe_diff:.2e} (tolerance: {tolerance:.2e})")
                        test_passed = False
                        failure_reasons.append(f"ZPE mismatch (diff: {zpe_diff:.2e})")
                    elif verbose:
                        print(f"    ZPE: {calculated_zpe:.12f} (diff: {zpe_diff:.2e})")

            # Validate Gibbs correction if present
            if "expected_gibbs_correction" in test:
                expected_gibbs = test["expected_gibbs_correction"]
                calculated_gibbs = extract_gibbs_correction(output_data)

                if calculated_gibbs is None:
                    print(f"  {Colors.RED}✗ FAILED{Colors.RESET} - Could not extract Gibbs correction from JSON")
                    test_passed = False
                    failure_reasons.append("Missing Gibbs correction")
                else:
                    gibbs_diff = abs(calculated_gibbs - expected_gibbs)
                    if not validate_energy(calculated_gibbs, expected_gibbs, tolerance):
                        print(f"  {Colors.RED}✗ FAILED{Colors.RESET} - Gibbs correction mismatch")
                        print(f"    Expected:   {expected_gibbs:.12f}")
                        print(f"    Calculated: {calculated_gibbs:.12f}")
                        print(f"    Difference: {gibbs_diff:.2e} (tolerance: {tolerance:.2e})")
                        test_passed = False
                        failure_reasons.append(f"Gibbs mismatch (diff: {gibbs_diff:.2e})")
                    elif verbose:
                        print(f"    Gibbs correction: {calculated_gibbs:.12f} (diff: {gibbs_diff:.2e})")

            # Report final test status
            if test_passed:
                print(f"  {Colors.GREEN}✓ PASSED{Colors.RESET}")
                if verbose:
                    print(f"    Energy: {calculated:.12f} (diff: {diff:.2e})")
                print()
                passed += 1
            else:
                print()
                failed += 1
                errors.append((test_name, "; ".join(failure_reasons)))

    return passed, failed, errors


def main():
    """Main entry point"""
    import argparse

    parser = argparse.ArgumentParser(description="Run metalquicha validation tests")
    parser.add_argument("--manifest", default="validation_tests.json",
                       help="Path to validation manifest JSON")
    parser.add_argument("--exe", default="../build/mqc",
                       help="Path to metalquicha executable")
    parser.add_argument("--prep-script", default="../mqc_prep.py",
                       help="Path to mqc_prep.py script")
    parser.add_argument("--validation-dir", default=".",
                       help="Path to validation directory containing inputs/")
    parser.add_argument("--skip-prep", action="store_true",
                       help="Skip .mqc file preparation step")
    parser.add_argument("--mpi", action="store_true",
                       help="Run tests using mpirun (fragmented: -np N, unfragmented: -np 1)")
    parser.add_argument("--np", type=int, default=4,
                       help="Number of MPI processes for fragmented calculations (default: 4)")
    parser.add_argument("-v", "--verbose", action="store_true",
                       help="Verbose output")
    parser.add_argument("-t", "--test", type=str, default=None,
                       help="Run only tests matching this name (substring match)")

    args = parser.parse_args()

    # Clean up any existing output files to ensure fresh test runs
    cleanup_output_files(verbose=args.verbose)

    # Prepare .mqc files from JSON inputs (unless skipped)
    if not args.skip_prep:
        prepare_mqc_files(
            validation_dir=args.validation_dir,
            prep_script=args.prep_script,
            verbose=args.verbose
        )

    # Run tests
    passed, failed, errors = run_validation_tests(
        manifest_file=args.manifest,
        exe_path=args.exe,
        verbose=args.verbose,
        use_mpi=args.mpi,
        nprocs=args.np,
        test_filter=args.test
    )

    # Summary
    total = passed + failed
    print(f"{Colors.BOLD}{'='*60}{Colors.RESET}")
    print(f"{Colors.BOLD}Validation Summary{Colors.RESET}")
    print(f"{'='*60}")
    print(f"Total tests: {total}")
    print(f"{Colors.GREEN}Passed:      {passed}{Colors.RESET}")
    print(f"{Colors.RED}Failed:      {failed}{Colors.RESET}")

    if failed > 0:
        print(f"\n{Colors.BOLD}Failed tests:{Colors.RESET}")
        for test_name, reason in errors:
            print(f"  • {test_name}: {reason}")

    print(f"{'='*60}\n")

    # Exit with appropriate code for CI
    sys.exit(0 if failed == 0 else 1)


if __name__ == "__main__":
    main()
