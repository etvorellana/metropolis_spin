import numpy.random as rdm
import numpy as np
import matplotlib.pyplot as plt
import scipy.fftpack

# estado inicial
def inicial(N):
    state = 2*rdm.randint(2,size=N)-1
    return state
 
# campo aleatório gaussiano
#def grf(N):
def grf(N, alpha=1.0, delta=1.0, flag_normalize=True):
    #alpha = 1.0
    #delta = 1.0
    #flag_normalize = True

    #k_idx = np.mgrid[:N] - int( (N + 1)/2 )
    k_idx = np.mgrid[:N] - (N + 1)//2
    k_idx = scipy.fftpack.fftshift(k_idx)   # shift the zero frequency component to the center of the spectrum

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
#def spin_flip(beta,config):
def spin_flip(N, beta,config):
    nb = 0
    #for i in range(N):
    rindex = rdm.randint(N, size=N)
    for a in rindex:
        #a = rdm.randint(N)
        s = config[a]
        #for j in range(N):
        #    if j != a:
        #        nb += config[(j)%N]
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
    h = h * config
    for i in range(N):
        for j in range(i,N):
            if i != j:
                e += -config[i]*config[j]/(abs(i-j)**alpha) - h[i]*config[i]
    return e

def main():
    pass

if __name__ == "__main__":
    main()