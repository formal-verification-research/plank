# Description
'''
The purpose of Values is to store all the determined values taken from the plank paper to be used as
coefficeints for calculating probabilities
'''


# Dimensionless Coefficients

# Length of Domain
L = 0.05
# EC Diffusion Coefficient
Dp = 3.6 * 10 ** -6
# VEGF Diffusion Coefficient
Dv = 3.6 * 10 ** -5
# Fibronectin Diffusion Coefficient
Df = 3.6 * 10 ** -10
# Kinetic Parameters for VEGF
lamda1 = 74.769231
v1 = 0.007692308
# Kinetic Parameters for Fibronectin
lamda3 = 19.277108
v3 = 1.2048193
# Protease Inhibitor Equilibrium Constant
ve = 1
# Protease Decay Rate
mu = 4.56
# Fibronectin Production Times
Tf = 18
# Proliferation Rates
Q = 8 * 10 ** -4
A = 44.13
lamda0 = 1.1 * 10 ** -9
m1 = 2
# Death Rate
mu1 = 7.142857 * 10 ** -5
# Thresehold Values
f0 = 6 * 10 ** -3
Tiv = 4.5
# VEGF Transmission Rates
B1 = 1
psi = 2
# VEGF Source Term Constants
v0 = 0.04
sigma = 1.514705513 * 10 ** -3
# Initial EC Density
p0 = 10 ** -5
# Initial Fibronectin Level
f0 = 10 ** -2
# Number of Receptors per Cell
deltae = 10 ** 5
# Protease Transition Probability Parameters
alpha1 = 0.1
alpha2 = 1
alpha3 = 0.1
alpha4 = 1
# Fibronectin Transition Probability Parameters
beta1 = 1
beta2 = 0.1
beta3 = 1
beta4 = 0.5
# VEGF Transition Probability Parameters
delta1 = 0.1
delta2 = 1
delta3 = 0.1
delta4 = 1
# Transition Probability Exponents
gamma1 = 100
gamma2 = 100
gamma3 = 40
gamma4 = 50
gamma5 = 37.5
gamma6 = 20
# K values
K1 = L ** 2 * lamda1 * deltae * p0 / Dp
K2 = L ** 2 * B1 / Dp
K3 = L ** 2 * mu / Dp
K4 = 4 * L ** 2 / (Dp * Tf)
K5 = L ** 2 * lamda3 / (Dp * v1)
K6 = f0 * v3
K8 = Dp * Tiv / L ** 2
K10 = alpha1 * v1
K11 = alpha2 * v1
K12 = beta1 / f0
K13 = beta2 / f0
K14 = delta1 * v1
K15 = delta2 * v1
K17 = Dp / Dp
K18 = L ** 2 * Q / Dp
K20 = L ** 2 * mu1 / Dp
K21 = Dv / Dp
K22 = Df / Dp
K23 = 4 * L ** 2 / (Dp * Tf)
K25 = A / v1
K26 = lamda0 / (v1 * m1)
K27 = alpha3 * v1
K28 = alpha4 * v1
K29 = beta3 / f0
K30 = beta4 / f0
K31 = delta3 * v1
K32 = delta4 * v1
K33 = Dv / (L * psi)
K35 = v0 * sigma * v1 / Dv
# Exponent for vegf source
m0 = 12
#Overrelaxation variables
relax1 = 1.45
relax2 = 1
