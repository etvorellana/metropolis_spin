
import numpy.random as rdm
import numpy as np
import scipy.fftpack

# estado inicial
def inicial(N):
    ''' 
        Gera o estado inicial do sistema. Um array, de tamanho N, preenchidos com valores
        1 ou -1 
    '''
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
def spin_flip(N, beta,config):
    #nb = 0
    for i in range(N):
        a = rdm.randint(N)
        s = config[a]
        nb = 0
        for j in range(N):
            if j != a:
                nb += config[(j)%N]
        cost = 2*s*nb # cost = 2*E (E1 - E0 = 2E)
        if cost < 0:
            s *= -1
        #elif rdm.rand() < math.exp(-cost*beta):    #antes usa a função exp do pacote math
        elif rdm.rand() < np.exp(-cost*beta):       #depois usa a função exp do pacote numpy  
            s *= -1
        config[a] = s
    return config

# Energia
def energy(N, config,alpha):
    e = 0
    h = grf(N)
    for i in range(N):
        for j in range(i,N):
            if i != j:
                e += -config[i]*config[j]/(abs(i-j)**alpha) - h[i]*config[i]
    return e

# Magnetizacao
def magnetization(config):
    m = np.sum(config)
    return m

# Densidade de defeito
def defeito(config):
    d = 0
    h = grf(N)   # Aqui se pode usar grfM(N) para testar
    for i in range(N):
        if config[i] != h[i]:
            d += 1
    return d

