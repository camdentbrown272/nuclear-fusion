"""Regenerate every M7 number and figure (about 1.5-2 h single-threaded).
Order matters: m7_wellcal reads the muon-shower rate written by m7_background."""
import os
import runpy

HERE = os.path.dirname(os.path.abspath(__file__))
for s in ['m7_common', 'm7_verify', 'm7_ipc', 'm7_transport', 'm7_detectors', 'm7_background',
          'm7_wellcal', 'm7_calibration', 'm7_modulation']:
    print(f'===== {s} =====', flush=True)
    runpy.run_path(os.path.join(HERE, s + '.py'), run_name='__main__')
