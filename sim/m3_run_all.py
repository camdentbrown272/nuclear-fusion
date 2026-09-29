"""Regenerate every M3 number and figure:  python3 sim/m3_run_all.py"""
import m3_loading
import m3_membrane
import m3_mechanics
import m3_cycling
import m3_sensitivity

if __name__ == "__main__":
    for mod in (m3_loading, m3_membrane, m3_mechanics, m3_cycling, m3_sensitivity):
        print(f"\n######## {mod.__name__} ########")
        mod.main()
