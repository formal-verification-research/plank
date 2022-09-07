"""
Description:
The parameter_vault file stores all the parameters entered into the simulation model. Most of these parameters are
taken from the Plank Paper (Appendix B, pgs 180-181), but some have been decided by the Vargis lab.
"""


'''
Parameters (Ordered according to group and alphabet to make them simple to find and understand)
'''


# Anastomosis Toggle (True/False)
anastomotic = True


# AMD Source Parameters
# Exponent for VEGF and PEDF
M0 = 12
# PEDF
D0 = 3 * 10 ** -6  # uM/h
# VEGF
SIGMA = 1.514705513 * 10**-3  # plank = 1.514705513 * 10**-3; vargis = 4.97512 * 10**-3
V0 = 2.5 * 10 ** -6  # uM*mm^2/h, plank = 0.04; vargis = 6.8442 * 10**-6


# Basic Parameters
division = 1  # How many simulated hours must pass before an EC is permitted to divide again
graph_time = 450  # How often a graph is created in amount of time step
L = 0.05  # The physical length of the x domain in mm
max_cells_allowed = 100  # How many total EC are allowed in the simulation
number_of_cells = 5  # How many EC the model starts with in the parent blood vessel
simulation_time = 48  # The maximum time the simulation will last in hours
threshold = 0.6  # The amount that the fib must drop to in the blood vessel wall before the EC can escape
time_step_duration = 8  # How much time passes in each time step of the discrete model in seconds
tolerance = 0.001  # Used to determine when the iteration method is close enough to the true value
x_steps = 201  # How many mesh points lie along the x axis for the substrate and EC matrices


# Diffusion Parameters
# EC
DP = 3.6 * 10**-6  # mm^2/h
# Fibronectin
DF = 3.6 * 10**-10  # mm^2/h
# PEDF
DD = 5 * 10**-5  # mm^2/h
# VEGF
DV = 3.6 * 10**-5  # mm^2/h, plank = 3.6 * 10**-5; vargis = 0.374


# Initial Factors
DELTAE = 10**5  # Number of Receptors per Cell
f0 = 10**-2  # Initial Fibronectin Level uM
P0 = 10**-5  # Initial EC Density uM


# Kinetic Parameters
# Fibronectin
LAMDA3 = 19.277108  # 1/(uM*h)
V3 = 1.2048193  # 1/uM
# PEDF
D1 = 3.75 * 10 ** -4  # 1/uM
LAMDA2 = 74.769231  # 1/h
# Protease
MU = 4.56  # 1/h, Protease Decay Rate
VE = 1  # 1/uM, Inhibitor Equilibrium Constant
# VEGF
LAMDA1 = 74.769231  # 1/h
V1 = 3.75 * 10 ** -4  # 1/uM


# Over-relaxation Variables used for VEGF and PEDF Iteration
RELAX1 = 1.45
RELAX2 = 1.0


# Production Parameters
# Death Rate
MU1 = 7.142857 * 10**-5  # 1/h
# Fibronectin Production Time
TF = 18  # h
# Proliferation Rates
A = 44.13  # 1/uM
LAMDA0 = 1.1 * 10**-9  # uM^-2
M1 = 2
Q = 8 * 10**-4  # 1/h


# Transition Probabilities
# Exponents
GAMMA1 = 100
GAMMA2 = 100
GAMMA3 = 40
GAMMA4 = 50
GAMMA5 = 37.5
GAMMA6 = 20
# Fibronectin
BETA1 = 1
BETA2 = 0.1
BETA3 = 1
BETA4 = 0.5
# Protease
ALPHA1 = 0.1
ALPHA2 = 1
ALPHA3 = 0.1
ALPHA4 = 1
# VEGF
DELTA1 = 0.1
DELTA2 = 1
DELTA3 = 0.1
DELTA4 = 1


# Transmission Rate Parameters
# VEGF
B1 = 1  # 1/h
PSI = 2  # mm/h
# PEDF
B2 = 1  # 1/h
PSID = 2  # mm/h
# Threshold Values
F0 = 6 * 10**-3  # uM
TIV = 4.5  # h


'''
Non-Dimensionalization (See Plank Paper Page 150 and 153)
'''


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
K35 = V0 * SIGMA * L ** 2 / (V1 * DV)
K36 = L**2 * LAMDA2 * DELTAE * P0 / DP
K37 = L**2 * B2 / DP
K38 = D0 * SIGMA * L ** 2 / (DD * D1)
K39 = DD / DP
K40 = DD / (L * PSID)
