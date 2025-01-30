# canonical - long-range - nofield (clr)
import numpy as np
from numba import jit

##### Functions #####

# Initial state
@jit(forceobj=True)
def inicial(N):
    state = 2*np.random.randint(2,size=N)-1
    return state

# Spin flip
@jit(forceobj=True)
def spin_flip(N,beta,config,alpha):
    rindex = np.random.randint(N, size=N)
    ji = np.arange(1,N)**alpha
    for a in rindex:
        s = config[a]
        jr = ji[:a]
        nb = (config[:a]/jr[::-1]).sum() + (config[a+1:]/ji[:N-a-1]).sum()
        cost = 2*s*nb # cost = 2*E (E1 - E0 = 2E)
        if cost < 0:
            s *= -1
        elif np.random.rand() < np.exp(-cost*beta):
            s *= -1
        config[a] = s
    return config

# Energy
@jit(forceobj=True)
def energy(N, config,alpha):
    e = 0
    ji = np.arange(1,N)**alpha
    for i in range(N):
        e += (-config[i]*config[i+1:]/ji[:N-i-1]).sum()
    return e

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
            correlation[r] += config[i] * config[(r+i)%N]
    return correlation
    
##### Parameters #####

nt     = 120        #  number of temperature points
N      = 2**10     #  size of the lattice, N
alpha = np.array([1.,1.4,1.8,2.0])        #  exponent of the long range interaction
eqSteps = 1000       #  number of MC sweeps for equilibration
mcSteps = eqSteps     #  number of MC sweeps for calculation

# valores de temperatura
T       = np.linspace(0.5, 20.0, nt)

#  arrays para armazenar os valores médios
E_1,M_1,C_1,X_1= np.zeros((len(alpha),nt)), np.zeros((len(alpha),nt)), np.zeros((len(alpha),nt)), np.zeros((len(alpha),nt))
E_2,M_2,C_2,X_2= np.zeros((len(alpha),nt)), np.zeros((len(alpha),nt)), np.zeros((len(alpha),nt)), np.zeros((len(alpha),nt))
E_4,M_4,C_4,X_4= np.zeros((len(alpha),nt)), np.zeros((len(alpha),nt)), np.zeros((len(alpha),nt)), np.zeros((len(alpha),nt))


n1, n2  = 1.0/(mcSteps*N), 1.0/(mcSteps*mcSteps*N)


##### Simulations #####

# Funcoes variando com a temperatura
for j in range(len(alpha)):

    for tt in range(nt):
        E1 = M1 = E2 = M2 = 0
        config = inicial(N)
        iT=1/T[tt]; iT2=iT*iT; # Termos referentes à temperatura => beta = 1/kT, k=1 é a constante de Boltzmann
        for i in range(eqSteps):         # equilibrate
            spin_flip(N,iT,config,alpha[j])         # Monte Carlo moves

        for i in range(mcSteps):
            spin_flip(N,iT,config,alpha[j])
            Ene = energy(N,config,alpha[j])     # calculate the energy
            Mag = magnetization(N,config)        # calculate the magnetisation

            E1 = E1 + Ene
            M1 = M1 + Mag
            M2 = M2 + Mag*Mag
            E2 = E2 + Ene*Ene

        E_1[j][tt] = E1*n1
        M_1[j][tt] = M1*n1
        C_1[j][tt] = (E2*n1 - E1*E1*n2)*iT2
        X_1[j][tt] = (M2*n1 - M1*M1*n2)*iT
        #Considero que já fiz 2000 passos de Monte Carlo
        #E1 = M1 = E2 = M2 = 0
        mcSteps = 2000 # para recalcular n1, n2
        n1, n2  = 1.0/(mcSteps*N), 1.0/(mcSteps*mcSteps*N)
        mcSteps = 1000 #Mais 1000 passos de Monte Carlo
        for i in range(mcSteps):
            spin_flip(N,iT,config,alpha[j])
            Ene = energy(N,config,alpha[j])     # calculate the energy
            Mag = magnetization(N,config)        # calculate the magnetisation

            E1 = E1 + Ene
            M1 = M1 + Mag
            M2 = M2 + Mag*Mag
            E2 = E2 + Ene*Ene
        
        E_2[j][tt] = E1*n1
        M_2[j][tt] = M1*n1
        C_2[j][tt] = (E2*n1 - E1*E1*n2)*iT2
        X_2[j][tt] = (M2*n1 - M1*M1*n2)*iT
        #Considero que já fiz 4000 passos de Monte Carlo 
        #E1 = M1 = E2 = M2 = 0
        mcSteps = 4000 # para recalcular n1, n2
        n1, n2  = 1.0/(mcSteps*N), 1.0/(mcSteps*mcSteps*N)
        mcSteps = 2000 # mais 2000 passos de Monte Carlo
        for i in range(mcSteps):
            spin_flip(N,iT,config,alpha[j])
            Ene = energy(N,config,alpha[j])     # calculate the energy
            Mag = magnetization(N,config)        # calculate the magnetisation

            E1 = E1 + Ene
            M1 = M1 + Mag
            M2 = M2 + Mag*Mag
            E2 = E2 + Ene*Ene
        
        E_4[j][tt] = E1*n1
        M_4[j][tt] = M1*n1
        C_4[j][tt] = (E2*n1 - E1*E1*n2)*iT2
        X_4[j][tt] = (M2*n1 - M1*M1*n2)*iT
        

    ##### Output (.csv) #####

    data1 = {"energy": E_1, "magnetization": M_1, "specific_heat": C_1, "susceptibility": X_1,"temperature": T}
    data2 = {"energy": E_2, "magnetization": M_2, "specific_heat": C_2, "susceptibility": X_2,"temperature": T}
    data4 = {"energy": E_4, "magnetization": M_4, "specific_heat": C_4, "susceptibility": X_4,"temperature": T}
    with open("clr{}-alpha{}-1000.csv".format(N,alpha[j]), "w") as f:
        f.write("energy,magnetization,specific_heat,susceptibility,temperature\n")
        for i in range(nt):
            f.write("{},{},{},{},{}\n".format(data1["energy"][j][i], data1["magnetization"][j][i],
                                                data1["specific_heat"][j][i], data1["susceptibility"][j][i],
                                                data1["temperature"][i]))
    with open("clr{}-alpha{}-2000.csv".format(N,alpha[j]), "w") as f:
        f.write("energy,magnetization,specific_heat,susceptibility,temperature\n")
        for i in range(nt):
            f.write("{},{},{},{},{}\n".format(data2["energy"][j][i], data2["magnetization"][j][i],
                                                data2["specific_heat"][j][i], data2["susceptibility"][j][i],
                                                data2["temperature"][i]))
    with open("clr{}-alpha{}-4000.csv".format(N,alpha[j]), "w") as f:
        f.write("energy,magnetization,specific_heat,susceptibility,temperature\n")
        for i in range(nt):
            f.write("{},{},{},{},{}\n".format(data4["energy"][j][i], data4["magnetization"][j][i],
                                                data4["specific_heat"][j][i], data4["susceptibility"][j][i],
                                                data4["temperature"][i]))


print("Done!")
   
