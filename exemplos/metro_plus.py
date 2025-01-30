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

def grf(N, alpha=1.0, delta=1.0, flag_normalize=True):
    
    '''
    shift the zero frequency component to the center of the spectrum
    numpy.mgrid = <numpy.lib.index_tricks.MGridClass object>
    An instance which returns a dense multi-dimensional “meshgrid”.
    '''
    
    #k_idx = np.mgrid[:N] - int( (N + 1)/2 )    # antes
    k_idx = np.mgrid[:N] - (N + 1)//2           # Gera o mesmo resultado
    k_idx = scipy.fftpack.fftshift(k_idx)

    amplitude = np.power( k_idx**2 + 1e-10, -alpha/4.0 )
    amplitude[0] = 0

    noise = rdm.normal(size=(N)) + 1j * rdm.normal(size=(N)) # complex noise

    gfield = np.fft.ifft(noise * amplitude).real # inverse Fourier transform

    # normalize the field to have zero mean and unit standard deviation
    if flag_normalize:
        gfield = gfield - np.mean(gfield)
        gfield = delta*gfield/np.std(gfield)

    return gfield

# Spin flip
def spin_flip(N, beta, config): # N é o tamanho do array, beta é a temperatura e config é o array de spins
    nb = 0
    rindex = rdm.randint(N, size=N) #Gera todos os índices aleatórios de uma vez
    #for i in range(N):
    for a in rindex:   # Percorre o array de índices aleatórios
        #a = rdm.randint(N)
        s = config[a]
        #for j in range(N):
        #    if j != a:
        #        nb += config[(j)%N]
        #nb += config[:a].sum() + config[(a+1):].sum()
        nb = config[:a].sum() + config[(a+1):].sum()
        cost = 2*s*nb # cost = 2*E (E1 - E0 = 2E)
        if cost < 0:
            s *= -1
        elif rdm.rand() < np.exp(-cost*beta):
            s *= -1
        config[a] = s
    return config

# Energia
def energy(N, config,alpha):
    e = 0
    h = grf(N)
    h = h*config
    ji = np.arange(1,N)**alpha
    for i in range(N):
        #for j in range(i+1,N):
        #    e += -config[i]*config[j]/((j-i)**alpha) - h[i]
        e += (-config[i]*config[i+1:]/ji[:N-i-1] - h[i]).sum()
    return e

# Magnetizacao
def magnetization(config):
    m = np.sum(config)
    return m