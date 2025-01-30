# canonical-long-range-field (clrf)

import numpy as np
from numba import jit
import scipy.fftpack

##### Functions #####

# Initial state
@jit(forceobj=True)
def inicial(N):
    state = 2*np.random.randint(2,size=N)-1
    return state

# Gaussian random field
@jit(forceobj=True)
def grf(N, delta, r = 0.001, flag_normalize = True):

    k_idx = np.mgrid[:N] - (N + 1)//2
    k_idx = scipy.fftpack.fftshift(k_idx)

    amplitude = np.power( k_idx**2 + 1e-10, - r/4.0 )
    amplitude[0] = 0

    noise = np.random.normal(size=(N)) \
        + 1j * np.random.normal(size=(N))

    gfield = np.fft.ifft(noise * amplitude).real

    if flag_normalize:
        gfield = gfield - np.mean(gfield)
        gfield = delta*gfield/np.std(gfield)

    return gfield     # campo aleatório com média zero
    # return abs(gfield)  # campo aleatório com média positiva

# Spin flip
@jit(forceobj=True)
def spin_flip(N,beta,config,alpha,delta):
    rindex = np.random.randint(N, size=N)
    ji = np.arange(1,N)**alpha
    H = grf(N,delta)
    for a in rindex:
        s = config[a]
        jr = ji[:a]
        nb = (config[:a]/jr[::-1]).sum() + (config[a+1:]/ji[:N-a-1]).sum()
        cost = 2*s*(nb + H[a]) # cost = 2*E (E1 - E0 = 2E)
        if cost < 0:
            s *= -1
        elif np.random.rand() < np.exp(-cost*beta):
            s *= -1
        config[a] = s
    return config

# Energy
@jit(forceobj=True)
def energy(N, config,alpha,delta):
    e = 0
    h = grf(N,delta)
    H = (h*config).sum()
    ji = np.arange(1,N)**alpha
    for i in range(N):
        # e += (-config[i]*config[i+1:]/ji[:N-i-1] - H[i]*config[i]).sum()
        e += (-config[i]*config[i+1:]/ji[:N-i-1]).sum()
    return e - H

# Magnetization
@jit(forceobj=True)
def magnetization(N,config):
    m = config.sum()
    return m


# Correlation function
@jit(forceobj=True)
def correlation_function(config,eqSteps):
    N = len(config)
    correlation = np.zeros(N//2)
    for r in range(N//2):
        for i in range(N):
            correlation[r] += config[i] * config[(i+r) % N]
        correlation[r] /= eqSteps
    return correlation
    
    
##### Parameters #####

nt     = 100        #  number of temperature points
N      = 2**10       #  size of the lattice, N
alpha = np.array([1.0,1.4,1.8,2.0])        #  exponent of the long range interaction.
delta = 0.5           #  external field
eqSteps = 1000       #  number of MC sweeps for equilibration
mcSteps = eqSteps       #  number of MC sweeps for calculation

# valores de temperatura
T       = np.linspace(0.5, 20.0, nt);

#  arrays para armazenar os valores médios
E,M,C,X,Cr,Rho = np.zeros((len(alpha),nt)), np.zeros((len(alpha),nt)), np.zeros((len(alpha),nt)), np.zeros((len(alpha),nt)), np.zeros(nt), np.zeros(nt)

n1, n2  = 1.0/(mcSteps*N), 1.0/(mcSteps*mcSteps*N)


##### Simulations #####

# Funcoes variando com a temperatura
for j in range(len(alpha)):

    for tt in range(nt):
        E1 = M1 = E2 = M2 = 0
        config = inicial(N)
        iT=1/T[tt]; iT2=iT*iT; # Termos referentes à temperatura => beta = 1/kT, k=1 é a constante de Boltzmann

        for i in range(eqSteps):         # equilibrate
            spin_flip(N,iT,config,alpha[j],delta)         # Monte Carlo moves

        for i in range(mcSteps):
            spin_flip(N,iT,config,alpha[j],delta)
            Ene = energy(N,config,alpha[j],delta)     # calculate the energy
            Mag = magnetization(N,config)        # calculate the magnetisation

            E1 = E1 + Ene
            M1 = M1 + Mag
            M2 = M2 + Mag*Mag
            E2 = E2 + Ene*Ene

        E[j][tt] = E1*n1
        M[j][tt] = M1*n1
        C[j][tt] = (E2*n1 - E1*E1*n2)*iT2
        X[j][tt] = (M2*n1 - M1*M1*n2)*iT


    ##### Output (.csv) #####

    data = {"energy": E, "magnetization": M, "specific_heat": C, "susceptibility": X, "temperature": T}
    with open("clrf{}-alpha{}.csv".format(N,alpha[j]), "w") as f:
        f.write("energy,magnetization,specific_heat,susceptibility,temperature\n")
        for i in range(nt):
            f.write("{},{},{},{},{}\n".format(data["energy"][j][i],data["magnetization"][j][i],data["specific_heat"][j][i],data["susceptibility"][j][i],data["temperature"][i]))

