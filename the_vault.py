# Description
# the_vault stores all the constants taken from the plank paper to be used as coefficients for calculating
# probabilities. Found in the plank paper, Appendix B.


# Stored values
L = 0.02  # Length of Domain
DP = 3.6 * 10**-6  # EC Diffusion Coefficient
DV = 3.6 * 10**-5  # VEGF Diffusion Coefficient
DF = 3.6 * 10**-10  # Fibronectin Diffusion Coefficient
LAMDA1 = 74.769231  # Kinetic Parameters for VEGF
V1 = 0.007692308  # Kinetic Parameters for VEGF
LAMDA3 = 19.277108  # Kinetic Parameters for Fibronectin
V3 = 1.2048193  # Kinetic Parameters for Fibronectin
VE = 1  # Protease Inhibitor Equilibrium Constant
MU = 4.56  # Protease Decay Rate
TF = 18  # Fibronectin Production Times
Q = 8 * 10**-4  # Proliferation Rates
A = 44.13  # Proliferation Rates
LAMDA0 = 1.1 * 10**-9  # Proliferation Rates
M1 = 2  # Proliferation Rates
MU1 = 7.142857 * 10**-5  # Death Rate
F0 = 6 * 10**-3  # Thresehold Values
TIV = 4.5  # Thresehold Values
B1 = 1  # VEGF Transmission Rates
PSI = 2  # VEGF Transmission Rates
V0 = 0.04  # VEGF Source Term Constants
SIGMA = 1.514705513 * 10**-3  # VEGF Source Term Constants
P0 = 10**-5  # Initial EC Density
F0 = 10**-2  # Initial Fibronectin Level
M0 = 12  # Exponent for vegf source
RELAX1 = 1.45  # Overrelaxation variables
RELAX2 = 1  # Overrelaxation variables
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

# K values
K1 = L**2 * LAMDA1 * DELTAE * P0 / DP
K2 = L**2 * B1 / DP
K3 = L**2 * MU / DP
K4 = 4 * L**2 / (DP * TF)
K5 = L**2 * LAMDA3 / (DP * V1)
K6 = F0 * V3
K8 = DP * TIV / L**2
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
K25 = A / V1
K26 = LAMDA0 / (V1 * M1)
K27 = ALPHA3 * V1
K28 = ALPHA4 * V1
K29 = BETA3 / F0
K30 = BETA4 / F0
K31 = DELTA3 * V1
K32 = DELTA4 * V1
K33 = DV / (L * PSI)
K35 = V0 * SIGMA * V1 / DV
