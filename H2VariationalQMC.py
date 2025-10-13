import numpy as np
import matplotlib.pyplot as plot

s = 0.74/0.529          # Angstroms
alpha = 2               # given
a = 0.840751            # Calculated in Mathematica
beta = np.linspace(0.2, 2.0, 19)
delta = 0.6             # Iterator
zeta = 2*13.6056981     # 2 times the Rydberg constant
Ntot = 10000               # Total number of runs
epsCoords = np.zeros(len(beta))

rL = np.array([-s/2, 0, 0])      # Position of left proton
rR = np.array([s/2, 0, 0])       # Position of right proton

r1 = np.array([np.random.uniform(-1/2, 1/2), np.random.uniform(-1/2, 1/2), np.random.uniform(-1/2, 1/2)])     # Position of left electron
r2 = np.array([np.random.uniform(-1/2, 1/2), np.random.uniform(-1/2, 1/2), np.random.uniform(-1/2, 1/2)])     # Position of right electron

def phi1(r1):
    r1L = rL + r1; magr1L = np.linalg.norm(r1L)  # Position of left electron relative to left proton
    r1R = rR + r1; magr1R = np.linalg.norm(r1R)  # Position of left electron relative to right proton

    Phi1 = np.exp(-magr1L/a) + np.exp(-magr1R/a)
    return Phi1, r1L, r1R, magr1L, magr1R

def phi2(r2):
    r2L = rL + r2; magr2L = np.linalg.norm(r2L)   # Position of right electron relative to left proton
    r2R = rR + r2; magr2R = np.linalg.norm(r2R)   # Position of right electron relative to right proton

    Phi2 = np.exp(-magr2L/a) + np.exp(-magr2R/a)
    return Phi2, r2L, r2R, magr2L, magr2R

def f(r1, r2):
    r12 = r1 - r2; magr12 = np.linalg.norm(r12)   # Position of electrons relative to each other

    F = np.exp(magr12/(alpha*(1 + beta*magr12)))
    return F, r12, magr12

def psi(r1, r2, s):    # Defining the wavefunction
    Psi = phi1(r1)[0]*phi2(r2)[0]*f(r1, r2)[0]
    return Psi

def epsilon(r1, r2, beta):      # Defining the energy function
    Phi1 = phi1(r1)[0]
    r1L = phi1(r1)[1]; magr1L = phi1(r1)[3]
    r1R = phi1(r1)[2]; magr1R = phi1(r1)[4]
    Phi2 = phi2(r2)[0]
    r2L = phi2(r2)[1]; magr2L = phi2(r2)[3]
    r2R = phi2(r2)[2]; magr2R = phi2(r2)[4]
    r12 = f(r1, r2)[1]; magr12 = f(r1, r2)[2]

    eta = 1/(1+beta*magr12)

    eps1 = -(1/a**2) - (eta**3)*((eta/4) + 1/magr12)
    eps2 = (1/a)*((np.exp(-magr1L/a)/magr1L) + np.exp(-magr1R/a)/magr1R)
    eps3 = ((eta**2)/2*a)*((np.exp(-magr1L/a)*(np.dot(r1L, r12)/(magr1L*magr12))) + (np.exp(-magr1R/a)*(np.dot(r1R, r12)/(magr1R*magr12))))
    eps4 = (1/a)*((np.exp(-magr2L/a)/magr2L) + np.exp(-magr2R/a)/magr2R)
    eps5 = ((eta**2)/2*a)*((np.exp(-magr2L/a)*(np.dot(r2L, r12)/(magr2L*magr12))) + (np.exp(-magr2R/a)*(np.dot(r2R, r12)/(magr2R*magr12))))
    eps6 = -(1/magr1L)-(1/magr1R)-(1/magr2L)-(1/magr2R)+(1/magr12)

    Epsilon = eps1 + (1/Phi1)*(eps2 + eps3) + (1/Phi2)*(eps4 - eps5) + eps6
    return Epsilon

for h in range(len(beta)):
    Na = 0
    Ne = 0
    sumE = 0
    sumE2 = 0
    AcceptanceCounter = 0
    for i in range(Ntot):           # Loop for x1
        rx = np.random.uniform(-1/2, 1/2)
        P0 = np.random.uniform(0, 1)
        x1T = r1[0] + rx*delta
        r1T = np.array([x1T, r1[1], r1[2]])
        ratio = (psi(r1T, r2, s)**2)/(psi(r1, r2, s)**2)
        if ratio[h] >= 1 or ratio[h] >= P0:
            AcceptanceCounter += 1
            if AcceptanceCounter % 11 == 0:
                r1 = r1T
                E = epsilon(r1, r2, beta[h])
                sumE += E
                sumE2 += E**2
                Ne += 1
                Na += 1
        else:
            r1 = r1
            E = epsilon(r1, r2, beta[h])
            sumE += E
            sumE2 += E**2
            Ne += 1
    for j in range(Ntot):                      # Loop for y1
        rx = np.random.uniform(-1/2, 1/2)
        P0 = np.random.uniform(0, 1)
        y1T = r1[0] + rx*delta
        r1T = np.array([r1[0], y1T, r1[2]])
        ratio = (psi(r1T, r2, s)**2)/(psi(r1, r2, s)**2)
        if ratio[h] >= 1 or ratio[h] >= P0:
            AcceptanceCounter += 1
            if AcceptanceCounter % 11 == 0:
                r1 = r1T
                E = epsilon(r1, r2, beta[h])
                sumE += E
                sumE2 += E ** 2
                Ne += 1
                Na += 1
        else:
            r1 = r1
            E = epsilon(r1, r2, beta[h])
            sumE += E
            sumE2 += E ** 2
            Ne += 1
    for k in range(Ntot):                  # Loop for z1
        rx = np.random.uniform(-1/2, 1/2)
        P0 = np.random.uniform(0, 1)
        z1T = r1[0] + rx*delta
        r1T = np.array([r1[0], r1[1], z1T])
        ratio = (psi(r1T, r2, s)**2)/(psi(r1, r2, s)**2)
        if ratio[h] >= 1 or ratio[h] >= P0:
            AcceptanceCounter += 1
            if AcceptanceCounter % 11 == 0:
                r1 = r1T
                E = epsilon(r1, r2, beta[h])
                sumE += E
                sumE2 += E ** 2
                Ne += 1
                Na += 1
        else:
            r1 = r1
            E = epsilon(r1, r2, beta[h])
            sumE += E
            sumE2 += E ** 2
            Ne += 1

    for l in range(Ntot):                    # Loop for x2
        rx = np.random.uniform(-1/2, 1/2)
        P0 = np.random.uniform(0, 1)
        x2T = r2[0] + rx*delta
        r2T = np.array([x2T, r2[1], r2[2]])
        ratio = (psi(r1, r2T, s)**2)/(psi(r1, r2, s)**2)
        if ratio[h] >= 1 or ratio[h] >= P0:
            AcceptanceCounter += 1
            if AcceptanceCounter % 11 == 0:
                r2 = r2T
                E = epsilon(r1, r2, beta[h])
                sumE += E
                sumE2 += E ** 2
                Ne += 1
                Na += 1
        else:
            r2 = r2
            E = epsilon(r1, r2, beta[h])
            sumE += E
            sumE2 += E ** 2
            Ne += 1
    for m in range(Ntot):                   # Loop for y2
        rx = np.random.uniform(-1/2, 1/2)
        P0 = np.random.uniform(0, 1)
        y2T = r2[0] + rx*delta
        r2T = np.array([r2[0], y2T, r2[2]])
        ratio = (psi(r1, r2T, s)**2)/(psi(r1, r2, s)**2)
        if ratio[h] >= 1 or ratio[h] >= P0:
            AcceptanceCounter += 1
            if AcceptanceCounter % 11 == 0:
                r2 = r2T
                E = epsilon(r1, r2, beta[h])
                sumE += E
                sumE2 += E ** 2
                Ne += 1
                Na += 1
        else:
            r2 = r2
            E = epsilon(r1, r2, beta[h])
            sumE += E
            sumE2 += E ** 2
            Ne += 1
    for n in range(Ntot):                       # Loop for z2
        rx = np.random.uniform(-1 / 2, 1 / 2)
        P0 = np.random.uniform(0, 1)
        z2T = r2[0] + rx*delta
        r2T = np.array([r2[0], r2[1], z2T])
        ratio = (psi(r1, r2T, s)**2)/(psi(r1, r2, s)**2)
        if ratio[h] >= 1 or ratio[h] >= P0:
            AcceptanceCounter += 1
            if AcceptanceCounter % 11 == 0:
                r2 = r2T
                E = epsilon(r1, r2, beta[h])
                sumE += E
                sumE2 += E ** 2
                Ne += 1
                Na += 1
        else:
            r2 = r2
            E = epsilon(r1, r2, beta[h])
            sumE += E
            sumE2 += E ** 2
            Ne += 1
    epsCoords[h] = epsilon(r1, r2, beta[h])
    sigma = np.sqrt(((sumE2 / Ne) - (sumE / Ne) ** 2) / Ne)       # Calculating uncertainty
    Pa = (Na / Ne) * 100                                          # Calculating acceptance probability
    E0 = sumE / Ne
    E0f = (E0 + 1 / s) * zeta                                     # Calculating final energy in eV
    Ebind = abs(E0f) - zeta                                       # Calculating binding energy in eV
    print(Pa)
    print(Ebind, '\u00B1', sigma)

plot.plot(beta, abs(epsCoords))
plot.xlabel('\u03B2')
plot.ylabel('Energy')
plot.show()