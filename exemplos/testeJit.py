import numpy as np
import numpy.random as rdm
import matplotlib.pyplot as plt
import time
import metro as mt
import metro_plus as mp

def testeMetro():
    start = time.perf_counter()
    
    nt     = 100        #  number of temperature points
    N      = 2**4       #  size of the lattice, N
    alpha = 2.          #  exponent of the long range interaction
    eqSteps = 500       #  number of MC sweeps for equilibration
    mcSteps = 500       #  number of MC sweeps for calculation

    # valores de temperatura
    T       = np.linspace(1, 200, nt);


    #time    = np.linspace(0, 100, mcSteps)
    E,M,C,X,Cr,Rho = np.zeros(nt), np.zeros(nt), np.zeros(nt), np.zeros(nt), np.zeros(nt), np.zeros(nt)
    n1, n2  = 1.0/(mcSteps*N*N), 1.0/(mcSteps*mcSteps*N*N)

    rdm.seed(1234567890)

    # Funções variando com a temperatura
    #config = mt.inicial(N)
    
    for tt in range(nt):
        E1 = M1 = E2 = M2 = 0
        config = mt.inicial(N)  # Este ponto precisa ser discutido
        iT=1/T[tt]; iT2=iT*iT; # Termos referentes à temperatura => beta = 1/kT, k=1 é a constante de Boltzmann

        for i in range(eqSteps):         # equilibrate
            mt.spin_flip(N, iT,config)         # Monte Carlo moves

        for i in range(mcSteps):
            mt.spin_flip(N, iT,config)
            Ene = mt.energy(N, config, alpha)     # calculate the energy
            Mag = mt.magnetization(config)        # calculate the magnetisation

            E1 = E1 + Ene
            M1 = M1 + Mag
            M2 = M2 + Mag*Mag
            E2 = E2 + Ene*Ene

        E[tt] = E1*n1
        M[tt] = M1*n1
        C[tt] = (E2*n1 - E1*E1*n2)*iT2
        X[tt] = (M2*n1 - M1*M1*n2)*iT
    
    
    stop = time.perf_counter()
    return start, stop, T, E, M, C, X

    
def testeMetroPlus():
    start = time.perf_counter()
    
    nt     = 100        #  number of temperature points
    N      = 2**4       #  size of the lattice, N
    alpha = 2.          #  exponent of the long range interaction
    eqSteps = 500       #  number of MC sweeps for equilibration
    mcSteps = 500       #  number of MC sweeps for calculation

    # valores de temperatura
    T       = np.linspace(1, 200, nt);


    #time    = np.linspace(0, 100, mcSteps)
    E,M,C,X,Cr,Rho = np.zeros(nt), np.zeros(nt), np.zeros(nt), np.zeros(nt), np.zeros(nt), np.zeros(nt)
    n1, n2  = 1.0/(mcSteps*N*N), 1.0/(mcSteps*mcSteps*N*N)

    rdm.seed(1234567890)

    # Funções variando com a temperatura
    #config = mp.inicial(N)
    
    for tt in range(nt):
        E1 = M1 = E2 = M2 = 0
        config = mp.inicial(N)  # Este ponto precisa ser discutido
        iT=1/T[tt]; iT2=iT*iT; # Termos referentes à temperatura => beta = 1/kT, k=1 é a constante de Boltzmann

        for i in range(eqSteps):         # equilibrate
            mp.spin_flip(N, iT,config)         # Monte Carlo moves

        for i in range(mcSteps):
            mp.spin_flip(N, iT,config)
            Ene = mp.energy(N, config,alpha)     # calculate the energy
            Mag = mp.magnetization(config)        # calculate the magnetisation

            E1 = E1 + Ene
            M1 = M1 + Mag
            M2 = M2 + Mag*Mag
            E2 = E2 + Ene*Ene

        E[tt] = E1*n1
        M[tt] = M1*n1
        C[tt] = (E2*n1 - E1*E1*n2)*iT2
        X[tt] = (M2*n1 - M1*M1*n2)*iT

    stop = time.perf_counter()
    return start, stop, T, E, M, C, X

def main():
    
    perf = np.zeros(10)
    print("1, ", end="")
    
    start, stop, T, E, M, C, X = testeMetro()
    perf[0] = stop - start
    for i in range(9):
        print("{}, ".format(i+2), end="")
        start, stop, T, E, M, C, X = testeMetro()
        perf[i+1] = stop - start
        

    print("Elapsed time = {}s std = {}".format(perf.mean(), perf.std()))

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

    print("1, ", end="")
    start, stop, T, E, M, C, X = testeMetroPlus()
    perf[0] = stop - start
    for i in range(9):
        print("{}, ".format(i+2), end="")
        start, stop, T, E, M, C, X = testeMetroPlus()
        perf[i+1] = stop - start

    print("Elapsed time = {}s std = {}".format(perf.mean(), perf.std()))

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

if __name__ == "__main__":
    main()

