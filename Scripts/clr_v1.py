# canonical - long-range - nofield (clr)
import numpy as np
from numba import jit
import argparse

##### Functions #####

# Initial state
@jit(forceobj=True)
def inicial(N):
    state = 2*np.random.randint(2, size=N, dtype=np.int8) - 1
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
    
def args_input():
    p = argparse.ArgumentParser(description='cmd parameters')
    p.add_argument('--spins', type=int, default=1024)
    p.add_argument('--eqs', type=int, default=1000)
    p.add_argument('--mcs', type=int, default=1000)
    p.add_argument('--rseed', type=int, default=0) 
    p.add_argument('--nTemp', type=int, default=120)
    p.add_argument('--tMin', type=float, default=0.5)
    p.add_argument('--tMax', type=float, default=20.0)
    return p.parse_args()

def main():
    args = args_input()
    ##### Parameters #####

    #definindo a semente do gerador de números aleatórios
    if args.rseed != 0:
        np.random.seed(args.rseed)
    
    nt     = args.nTemp                 #  number of temperature points
    N      = args.spins                 #  size of the lattice, N
    alpha = np.array([1.,1.4,1.8,2.0])  #  exponent of the long range interaction
    eqSteps = args.eqs                  #  number of MC sweeps for equilibration
    mcSteps = args.mcs                  #  number of MC sweeps for calculation

    # valores de temperatura
    T       = np.linspace(args.tMin, args.tMax, nt)
    lAlpha = len(alpha)

    #  arrays para armazenar os valores médios
    E,M,C,X= np.zeros((lAlpha,nt)), np.zeros((lAlpha,nt)), np.zeros((lAlpha,nt)), np.zeros((lAlpha,nt))
    
    n1, n2  = 1.0/(mcSteps*N), 1.0/(mcSteps*mcSteps*N)

    ##### Simulations #####

    # Funcoes variando com a temperatura
    for j in range(lAlpha):
        print(f"alpha = {alpha[j]}") # So para acompanhar o andamento do código
        for tt in range(nt):
            E1 = M1 = E2 = M2 = 0
            config = inicial(N)
            iT=1.0/T[tt]; iT2=iT**2; # Termos referentes à temperatura => beta = 1/kT, k=1 é a constante de Boltzmann
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

            E[j][tt] = E1*n1
            M[j][tt] = M1*n1
            C[j][tt] = (E2*n1 - E1*E1*n2)*iT2
            X[j][tt] = (M2*n1 - M1*M1*n2)*iT
        
        ##### Output (.csv) #####

        data = {"energy": E, "magnetization": M, "specific_heat": C, "susceptibility": X,"temperature": T}
        
        with open("clr{}-alpha{}-{}-{}.csv".format(N,alpha[j],eqSteps,mcSteps), "w") as f:
            f.write("energy,magnetization,specific_heat,susceptibility,temperature\n")
            for i in range(nt):
                f.write("{},{},{},{},{}\n".format(data["energy"][j][i], data["magnetization"][j][i],
                                                    data["specific_heat"][j][i], data["susceptibility"][j][i],
                                                    data["temperature"][i]))
        

    print("Done!")


if __name__ == "__main__":
    main()


