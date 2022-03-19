# Description
# the_vault stores all the constants taken from the plank paper to be used as coefficients for calculating
# probabilities. Found in the plank paper, Appendix B. Pages 180-181.


# Stored values
L = 0.05  # Length of the x domain (mm)
DP = 3.6 * 10**-6  # EC Diffusion Coefficient (mm^2/h)
DV = 3.6 * 10**-5  # VEGF Diffusion Coefficient (mm^2/h) plank=3.6 * 10**-5, skeletal muscle=0.374
DF = 3.6 * 10**-10  # Fibronectin Diffusion Coefficient (mm^2/h)
DA = 6.5 * 10**-5  # Angiostatin Diffusion Coefficient (mm^2/h)
DD = 3.6 * 10**-5  # Assumed PEDF diffusion coefficient (mm^2/h)
LAMDA1 = 74.769231  # Kinetic Parameters for VEGF (1/h)
LAMDA2 = 74.769231  # Assumed Kinetic Parameters for PEDF (1/h)
V1 = 0.007692309  # Kinetic Parameters for VEGF (1/uM)
D1 = 0.007692309  # Assumed Kinetic Parameters for PEDF (1/uM)
LAMDA3 = 19.277108  # Kinetic Parameters for Fibronectin (1/uM*h)
V3 = 1.2048193  # Kinetic Parameters for Fibronectin (1/uM)
VE = 1  # Protease Inhibitor Equilibrium Constant (1/uM)
MU = 4.56  # Protease Decay Rate (1/h)
TF = 18  # Fibronectin Production Time (h)
TA = 1  # Angiostatin Decay Time (h)
Q = 8 * 10**-4  # Proliferation Rates (1/h)
A = 44.13  # Proliferation Rates (1/uM)
LAMDA0 = 1.1 * 10**-9  # Proliferation Rates (uM^-2), don't know why the -2
M1 = 2  # Proliferation Rates (no units)
MU1 = 7.142857 * 10**-5  # Death Rate (1/h)
F0 = 6 * 10**-3  # Threshold Values (uM)
TIV = 4.5  # Threshold Values (h)
B1 = 1  # VEGF Transmission Rates (1/h)
B2 = 1  # Assumed PEDF Transmission Rates (1/h)
PSI = 2  # VEGF Transmission Rates (mm/h)
PSID = 2  # Assumed PEDF Transmission Rates (mm/h)
PSIS = 2  # Angiostatin Transmission Rate (mm/h)
V0 = 0.04  # VEGF Source Term Constants (uM*mm^2/h) plank=0.04 lab says 6.8442*10**-6
D0 = 0.04  # Assumed PEDF Source Term Constants (uM*mm^2/h)
SIGMA = 1.514705513 * 10**-3  # VEGF Source Term Constants plank=1.514705513 * 10**-3, lab=4.97512*10**-3
AR = 1700  # Angiostatin Source Term Constant (uM/h)
P0 = 10**-5  # Initial EC Density (uM)
F0 = 10**-2  # Initial Fibronectin Level (uM)
M0 = 12  # Exponent for vegf and pedf sources
RELAX1 = 1.45  # Over-relaxation Variable Used For VEGF and PEDF Iteration (Unit-less)
RELAX2 = 1  # Over-relaxation Variable Used For Fibronectin Iteration (Unit-less)
DELTAE = 10**5  # Number of Receptors per Cell
ALPHA1 = 0.1  # Protease Transition Probability Parameters
ALPHA2 = 1  # Protease Transition Probability Parameters
ALPHA3 = 0.1  # Protease Transition Probability Parameters
ALPHA4 = 1  # Protease Transition Probability Parameters
BETA1 = 1  # Fibronectin Transition Probability Parameters
BETA2 = 0.1  # Fibronectin Transition Probability Parameters
BETA3 = 1  # Fibronectin Transition Probability Parameters
BETA4 = 0.5  # Fibronectin Transition Probability Parameters
DELTA1 = 0.1  # VEGF Transition Probability Parameters
DELTA2 = 1  # VEGF Transition Probability Parameters
DELTA3 = 0.1  # VEGF Transition Probability Parameters
DELTA4 = 1  # VEGF Transition Probability Parameters
GAMMA1 = 100  # Transition Probability Exponents
GAMMA2 = 100  # Transition Probability Exponents
GAMMA3 = 40  # Transition Probability Exponents
GAMMA4 = 50  # Transition Probability Exponents
GAMMA5 = 37.5  # Transition Probability Exponents
GAMMA6 = 20  # Transition Probability Exponents


# K values, Non-dimensionalized
K1 = L**2 * LAMDA1 * DELTAE * P0 / DP
K2 = L**2 * B1 / DP
K3 = L**2 * MU / DP
K4 = 4 * L**2 / (DP * TF)
K5 = L**2 * LAMDA3 / (DP * V1)
K6 = F0 * V3
K7 = VE * AR * L**2 / DP
K8 = DP * TIV / L**2
K9 = L**2 / (DP * TA)
K10 = ALPHA1 * V1
K11 = ALPHA2 * V1
K12 = BETA1 / F0
K13 = BETA2 / F0
K14 = DELTA1 * V1
K15 = DELTA2 * V1
K17 = DP / DP
K18 = L**2 * Q / DP
K20 = L**2 * MU1 / DP
K21 = DV / DP
K22 = DF / DP
K23 = 4 * L**2 / (DP * TF)
K24 = DA / DP
K25 = A / V1
K26 = LAMDA0 / (V1 * M1)
K27 = ALPHA3 * V1
K28 = ALPHA4 * V1
K29 = BETA3 / F0
K30 = BETA4 / F0
K31 = DELTA3 * V1
K32 = DELTA4 * V1
K33 = DV / (L * PSI)
K34 = DA / (L * PSIS)
K35 = V0 * SIGMA * V1 / DV
K36 = L**2 * LAMDA2 * DELTAE * P0 / DP
K37 = L**2 * B2 / DP
K38 = D0 * SIGMA * D1 / DD
K39 = DD / DP
K40 = DD / (L * PSID)
