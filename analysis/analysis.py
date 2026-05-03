import sys, subprocess
from pathlib import Path

energies = sys.argv[1] if len(sys.argv) > 1 else "energies.csv"
msd = sys.argv[2] if len(sys.argv) > 2 else "msd.csv"
plots = sys.argv[3] if len(sys.argv) > 3 else "plots"

base = Path(__file__).parent
for script, args in [
    ("energy_analysis.py", [energies, plots]),
    ("msd_analysis.py", [msd, plots]),
    ("temperature_check.py", [energies, plots]),
]:
    print(f"\n=== {script} ===")
    r = subprocess.run([sys.executable, base / script] + args)
    if r.returncode != 0:
        print(f"ERROR: {script} failed (code {r.returncode})")
