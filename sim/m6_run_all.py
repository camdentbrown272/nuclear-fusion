"""Regenerate every M6 number and figure: python3 sim/m6_run_all.py  (about 30-60 min on 4 CPUs;
m6_gratings.py dominates)."""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
for s in ("m6_spp.py", "m6_fields.py", "m6_modes.py", "m6_skins.py", "m6_sectors.py", "m6_psd.py", "m6_gratings.py"):
    print(f"=== {s}", flush=True)
    subprocess.run([sys.executable, os.path.join(HERE, s)], check=True, cwd=HERE)
