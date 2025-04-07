#!/bin/bash

# Run the application clr



# Interação de longo alcance

#    alpha = 0.25 a 1.0, em intervalos de 0.25
    
#    T = np.linspace(0.5, 500, nt)

# OBSERVAÇÃO: Testar primeiro para alpha = 0.25 e 1.25, para as quatro condições de campo, 
# para ter certeza de que o intervalo do eixo de temperatura está de acordo.

# Experimento 1. Tetando a escala de temperatura do longo alcance

#    p.add_argument('--spins', type=int, default=1024)
#    p.add_argument('--eqs', type=int, default=1000)
#    p.add_argument('--mcs', type=int, default=1000)
#    p.add_argument('--rseed', type=int, default=0) 
#    p.add_argument('--nTemp', type=int, default=120)
#    p.add_argument('--tMin', type=float, default=0.5)
#    p.add_argument('--tMax', type=float, default=20.0)
#    p.add_argument('--alpha', type=str, default='short') # short, long, shortT, longT 

time python clr_v2.py --spins 4096 --eqs 4000 --mcs 4000 --rseed 1234 --nTemp 120 --tMin 0.5 --tMax 500.0 --alpha longT

#   p.add_argument('--spins', type=int, default=1024)
#   p.add_argument('--eqs', type=int, default=1000)
#   p.add_argument('--mcs', type=int, default=1000)
#   p.add_argument('--rseed', type=int, default=0) 
#   p.add_argument('--nTemp', type=int, default=120)
#   p.add_argument('--tMin', type=float, default=0.5)
#   p.add_argument('--tMax', type=float, default=20.0)
#   p.add_argument('--alpha', type=str, default='short') # short, long, shortT, longT
#   p.add_argument('--field', type=str, default='zero') # zero, positive
#   p.add_argument('--delta', type=float, default=0.5)

#Campo externo aleatório positivo (com abs) e média zero (sem abs)
    
    #delta = 0.1, 1.0, 10.0

time python clrf_v2.py --spins 4096 --eqs 4000 --mcs 4000 --rseed 1234 --nTemp 120 --tMin 0.5 --tMax 500.0 --alpha longT --field 0 --delta 0.1
time python clrf_v2.py --spins 4096 --eqs 4000 --mcs 4000 --rseed 1234 --nTemp 120 --tMin 0.5 --tMax 500.0 --alpha longT --field 0 --delta 1.0
time python clrf_v2.py --spins 4096 --eqs 4000 --mcs 4000 --rseed 1234 --nTemp 120 --tMin 0.5 --tMax 500.0 --alpha longT --field 0 --delta 100.0

time python clrf_v2.py --spins 4096 --eqs 4000 --mcs 4000 --rseed 1234 --nTemp 120 --tMin 0.5 --tMax 500.0 --alpha longT --field 1 --delta 0.1
time python clrf_v2.py --spins 4096 --eqs 4000 --mcs 4000 --rseed 1234 --nTemp 120 --tMin 0.5 --tMax 500.0 --alpha longT --field 1 --delta 1.0
time python clrf_v2.py --spins 4096 --eqs 4000 --mcs 4000 --rseed 1234 --nTemp 120 --tMin 0.5 --tMax 500.0 --alpha longT --field 1 --delta 100.0
