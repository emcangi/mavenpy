# Planck constant
hc_Jm = 1.98644568e-25  # J * m
hc_eVum = 1.239841  # eV um

eV_to_joules = 1.60217646e-19  # Joules per eV

# Constants for plasma moments:
c_kms = (2.99792458) * 1e5
# rest masses [eV/(km/s)^2]
mass = {"e-": (5.10998910) * 1e5 / c_kms ** 2,
        "H+": (938) * 1e6 / c_kms ** 2,
        "He++": (3.727) * 1e9 / c_kms ** 2}
