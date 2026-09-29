"""
M4 — shared physical constants and material parameters (SI units).

Every value carries its source.  Where a web page could not be fetched from the
modelling container (several publisher domains were blocked by the egress
proxy), the value is marked "handbook value, not re-fetched" and given a
generous uncertainty; none of the M4 recommendations hinge on those values
to better than the stated uncertainty (see M4 doc, section 7).
"""

F = 96485.332          # C/mol, Faraday constant (CODATA 2018, exact via e*NA)
R_GAS = 8.314462       # J/(mol K) (CODATA 2018)
K_B = 1.380649e-23     # J/K
SIGMA_SB = 5.670374e-8  # W/(m^2 K^4)
EV = 1.602176634e-19   # J
N_A = 6.02214076e23

# ---------------------------------------------------------------- thermochemistry
# Standard enthalpy of formation of the liquids (298.15 K).
# H2O: -285.83 kJ/mol  (NIST WebBook, https://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185)
# D2O: -294.60 kJ/mol  (NIST WebBook, https://webbook.nist.gov/cgi/cbook.cgi?ID=C7789200;
#       also used by Miles, JCMNS 33 (2020) 74, https://jcmns.org/article/72550.pdf , E_H = 1.5267 V)
DHF_H2O = 285.83e3     # J/mol (magnitude)
DHF_D2O = 294.60e3     # J/mol (magnitude)
E_TN_H2O = DHF_H2O / (2 * F)   # 1.4812 V thermoneutral potential
E_TN_D2O = DHF_D2O / (2 * F)   # 1.5267 V
# Enthalpy of D absorption into beta-PdD (plateau), -34.6 kJ/mol D2 = -17.3 kJ/mol D
# Flanagan et al., calorimetric enthalpies, J. Less-Common Met. (1991)
# https://www.sciencedirect.com/science/article/abs/pii/0022508891904313  (+-1 kJ/mol D2)
DH_ABS_D = -17.3e3     # J/mol D
DH_ABS_H = -19.1e3     # J/mol H (-38.3 kJ/mol H2, same source family; +-1 kJ/mol H2)
# Vaporisation enthalpy of D2O at 25 C ~45.4 kJ/mol, H2O 44.0 kJ/mol
# (engineeringtoolbox heavy-water page https://www.engineeringtoolbox.com/heavy-water-thermodynamic-properties-d_2003.html)
DH_VAP_D2O = 45.4e3
DH_VAP_H2O = 44.0e3

# ---------------------------------------------------------------- nuclear (H2 hypothesis)
Q_DD_HE4 = 23.85e6 * EV   # J per 4He from D+D -> 4He (mass difference)
HE_AIR_FRAC = 5.24e-6     # volume fraction of He in dry air (Glueckauf 1946; standard atmosphere value)

# ---------------------------------------------------------------- materials (k W/mK, rho*c J/m^3K)
# Handbook values (CRC Handbook / engineeringtoolbox.com tables); uncertainty +-5 % unless noted.
MAT = {
    #              k        rho*c
    "ss316":   (16.3,   8000 * 500),     # https://www.engineeringtoolbox.com/thermal-conductivity-metals-d_858.html
    "al6061":  (167.0,  2700 * 896),
    "cu":      (398.0,  8960 * 385),
    "pd":      (71.8,   12020 * 244),
    "pt":      (71.6,   21450 * 133),
    "ptfe":    (0.25,   2200 * 1000),    # PTFE 0.25 W/mK (+-10 %)
    "peek":    (0.25,   1300 * 1340),
    "d2o":     (0.595,  1104 * 4210),    # D2O liquid 25 C, engineeringtoolbox heavy-water page (above)
    "h2o":     (0.607,  997 * 4181),
    "gas":     (0.14,   0.16 * 7000),    # D2 gas k=0.14 W/mK (CRC); rho*c tiny
    "air":     (0.026,  1.2 * 1005),
    "foam":    (0.030,  40 * 1400),      # PU/PIR foam
    "pad":     (3.0,    2.5e6),          # silicone gap pad, 3 W/mK class (vendor datasheets, +-30 %)
    "cat":     (0.5,    1.0e6),          # Pt/Al2O3 or Pd/C pellet bed in SS mesh (estimate, +-100 %)
}

# ---------------------------------------------------------------- thermoelectric module (heat-flux sensor)
# 40x40 mm, 127-couple Bi2Te3 module (e.g. TEC1-12706 class).  Typical datasheet values:
# module Seebeck alpha ~ 0.05 V/K, thermal conductance K_m ~ 0.5 W/K, R ~ 2 ohm.
# (Handbook/datasheet values, not re-fetched; calibrated in situ anyway.)
TEC_ALPHA = 0.050      # V/K per module (+-10 %)
TEC_K = 0.50           # W/K per module (+-20 %)
TEC_R = 2.0            # ohm
TEC_A = 0.040 ** 2     # m^2
TEC_T = 0.0040         # m, thickness incl. ceramic plates
TEC_TEMPCO = 2e-3      # 1/K, relative dS/dT of Seebeck sensitivity (Bi2Te3 near 300 K; +-1e-3)

# ---------------------------------------------------------------- He permeation
# Pyrex 7740, 25 C: K ~ 1e-11 cm3(STP) mm / (s cm2 cmHg)   (Norton 1953, J. Am. Ceram. Soc.;
#   quoted 0.91-1.5e-11 in https://tf.nist.gov/general/pdf/2830.pdf ), activation energy 6.4 kcal/mol
#   (Altemose, J. Appl. Phys. 25, 868; https://pubs.aip.org/aip/jap/article/25/7/868/160802)
K_HE_PYREX_25C = 1.2e-11          # cm3STP*mm/(s cm2 cmHg), +-30 %
EA_HE_PYREX = 6.4e3 * 4.184       # J/mol
# Viton O-ring, ~20 C: K_He = 15.1e-15 m^2 s^-1 hPa^-1 (Sturm et al., JGR 2004,
#   https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2003JD004073)
K_HE_VITON = 15.1e-15             # m3(STP)*m/(m^2 s hPa)
# Henry constant of He in water at 25 C: 3.7e-4 mol/(kg bar) (Sander 2015 compilation, ACP 15, 4399)
KH_HE_WATER = 3.7e-4              # mol/(kg bar)
