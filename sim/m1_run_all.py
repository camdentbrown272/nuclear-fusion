"""Regenerate every M1 number and figure: python3 sim/m1_run_all.py"""
import m1_sites
import m1_nonthermal
import m1_ladder
import m1_flux

if __name__ == "__main__":
    m1_sites.main()
    m1_nonthermal.lines.clear()
    m1_nonthermal.main()
    m1_ladder.main()
    m1_flux.main()
