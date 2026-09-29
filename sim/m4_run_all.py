"""Regenerate every M4 number and figure:  python3 sim/m4_run_all.py"""
import os
import runpy

HERE = os.path.dirname(os.path.abspath(__file__))
for s in ["m4_thermal2d", "m4_power_meas", "m4_chem_he", "m4_noise", "m4_gasphase", "m4_dfm", "m4_budget"]:
    print(f"==== {s}")
    runpy.run_path(os.path.join(HERE, s + ".py"), run_name="__main__")
