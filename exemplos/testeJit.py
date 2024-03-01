import numpy.random as rdm
import numpy as np
import matplotlib.pyplot as plt
import scipy.fftpack


# estado inicial
def inicial(N):
    state = 2*rdm.randint(2,size=N)-1
    return state

# campo aleatório gaussiano
def grf(N):
    alpha = 1.0
    delta = 1.0
    flag_normalize = True

    k_idx = np.mgrid[:N] - int( (N + 1)/2 )
    k_idx = scipy.fftpack.fftshift(k_idx)

    amplitude = np.power( k_idx**2 + 1e-10, -alpha/4.0 )
    amplitude[0] = 0

    noise = rdm.normal(size=(N)) \
        + 1j * rdm.normal(size=(N))

    gfield = np.fft.ifft(noise * amplitude).real

    if flag_normalize:
        gfield = gfield - np.mean(gfield)
        gfield = delta*gfield/np.std(gfield)

    return gfield

# Spin flip
def spin_flip(beta,config):
    nb = 0
    for i in range(N):
        a = rdm.randint(N)
        s = config[a]
        for j in range(N):
            if j != a:
                nb += config[(j)%N]
        cost = 2*s*nb # cost = 2*E (E1 - E0 = 2E)
        if cost < 0:
            s *= -1
        elif rdm.rand() < np.exp(-cost*beta):
            s *= -1
        config[a] = s
    return config

# Energia
def energy(config,alpha):
    e = 0
    h = grf(N)
    for i in range(N):
        for j in range(i,N):
            if i != j:
                e += -config[i]*config[j]/(abs(i-j)**alpha) - h[i]*config[i]
    return e

# Magnetização
def magnetization(config):
    m = np.sum(config)
    return m

# Densidade de defeito
def defeito(config):
    d = 0
    h = grf(N)
    for i in range(N):
        if config[i] != h[i]:
            d += 1
    return d

rdm.seed(1234567890)
nt     = 100        #  number of temperature points
N      = 2**4       #  size of the lattice, N
alpha = 2.          #  exponent of the long range interaction
eqSteps = 500       #  number of MC sweeps for equilibration
mcSteps = 500       #  number of MC sweeps for calculation

# valores de temperatura
T       = np.linspace(1, 200, nt);


time    = np.linspace(0, 100, mcSteps)
E,M,C,X,Cr,Rho = np.zeros(nt), np.zeros(nt), np.zeros(nt), np.zeros(nt), np.zeros(nt), np.zeros(nt)
n1, n2  = 1.0/(mcSteps*N*N), 1.0/(mcSteps*mcSteps*N*N)

# Funções variando com a temperatura
config = inicial(N)
for tt in range(nt):
    E1 = M1 = E2 = M2 = 0
    # config = inicial(N)
    iT=1/T[tt]; iT2=iT*iT; # Termos referentes à temperatura => beta = 1/kT, k=1 é a constante de Boltzmann

    for i in range(eqSteps):         # equilibrate
        spin_flip(iT,config)         # Monte Carlo moves

    for i in range(mcSteps):
        spin_flip(iT,config)
        Ene = energy(config,alpha)     # calculate the energy
        Mag = magnetization(config)        # calculate the magnetization

        E1 = E1 + Ene
        M1 = M1 + Mag
        M2 = M2 + Mag*Mag
        E2 = E2 + Ene*Ene

        E[tt] = E1*n1
        M[tt] = M1*n1
        C[tt] = (E2*n1 - E1*E1*n2)*iT2
        X[tt] = (M2*n1 - M1*M1*n2)*iT

f = plt.figure(figsize=(18, 10)); # plot the calculated values

sp =  f.add_subplot(2, 2, 1 );
plt.plot(T, E, 'o', color="#A60628", label=' Energy');
plt.xlabel("Temperature (T)", fontsize=20);
plt.ylabel("Energy ", fontsize=20);
plt.legend(loc='best');

sp =  f.add_subplot(2, 2, 2 );
plt.plot(T, abs(M), '*', color="#348ABD", label='Magnetization');
plt.xlabel("Temperature (T)", fontsize=20);
plt.ylabel("Magnetization ", fontsize=20);
plt.legend(loc='best');

sp =  f.add_subplot(2, 2, 3 );
plt.plot(T, C, 'd', color="#A60628", label='Specific Heat');
plt.xlabel("Temperature (T)", fontsize=20);
plt.ylabel("Specific Heat ", fontsize=20);
plt.legend(loc='best');

sp =  f.add_subplot(2, 2, 4 );
plt.plot(T, X, 's', color="#348ABD", label='Susceptibility');
plt.xlabel("Temperature (T)", fontsize=20);
plt.ylabel("Susceptibility", fontsize=20);
plt.legend(loc='best');

plt.show()
