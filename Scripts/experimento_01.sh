#!/bin/bash

#    p.add_argument('--spins', type=int, default=1024)
#    p.add_argument('--eqs', type=int, default=1000)
#    p.add_argument('--mcs', type=int, default=1000)
#    p.add_argument('--rseed', type=int, default=0) 
#    p.add_argument('--nTemp', type=int, default=120)
#    p.add_argument('--tMin', type=float, default=0.5)
#    p.add_argument('--tMax', type=float, default=20.0)
#    p.add_argument('--alpha', type=str, default='short') # short, long, shortT, longT 


#time python clr_v2.py --spins 4096 --eqs 4096 --mcs 4096 --rseed 1234 --nTemp 120 --tMin 0.5 --tMax 20.0 --alpha short
#time python clr_v2.py --spins 4096 --eqs 4096 --mcs 4096 --rseed 1234 --nTemp 120 --tMin 0.5 --tMax 1000.0 --alpha long
#time python clr_v2.py --spins 4096 --eqs 4096 --mcs 4096 --rseed 1234 --nTemp 120 --tMin 0.5 --tMax 6000.0 --alpha zero

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

#time python clrf_v2.py --spins 4096 --eqs 4096 --mcs 4096 --rseed 1234 --nTemp 120 --tMin 0.5 --tMax 20.0 --alpha short --field 0 --delta 0.1
#time python clrf_v2.py --spins 4096 --eqs 4096 --mcs 4096 --rseed 1234 --nTemp 120 --tMin 0.5 --tMax 20.0 --alpha short --field 0 --delta 1.0
#time python clrf_v2.py --spins 4096 --eqs 4096 --mcs 4096 --rseed 1234 --nTemp 120 --tMin 0.5 --tMax 20.0 --alpha short --field 0 --delta 10.0

#time python clrf_v2.py --spins 4096 --eqs 4096 --mcs 4096 --rseed 1234 --nTemp 120 --tMin 0.5 --tMax 20.0 --alpha short --field 1 --delta 0.1
#time python clrf_v2.py --spins 4096 --eqs 4096 --mcs 4096 --rseed 1234 --nTemp 120 --tMin 0.5 --tMax 20.0 --alpha short --field 1 --delta 1.0
#time python clrf_v2.py --spins 4096 --eqs 4096 --mcs 4096 --rseed 1234 --nTemp 120 --tMin 0.5 --tMax 20.0 --alpha short --field 1 --delta 10.0

#time python clrf_v2.py --spins 4096 --eqs 4096 --mcs 4096 --rseed 1234 --nTemp 120 --tMin 0.5 --tMax 1000.0 --alpha long --field 0 --delta 0.1
#time python clrf_v2.py --spins 4096 --eqs 4096 --mcs 4096 --rseed 1234 --nTemp 120 --tMin 0.5 --tMax 1000.0 --alpha long --field 0 --delta 1.0
time python clrf_v2.py --spins 4096 --eqs 4096 --mcs 4096 --rseed 1234 --nTemp 120 --tMin 0.5 --tMax 1000.0 --alpha long --field 0 --delta 10.0

time python clrf_v2.py --spins 4096 --eqs 4096 --mcs 4096 --rseed 1234 --nTemp 120 --tMin 0.5 --tMax 1000.0 --alpha long --field 1 --delta 0.1
time python clrf_v2.py --spins 4096 --eqs 4096 --mcs 4096 --rseed 1234 --nTemp 120 --tMin 0.5 --tMax 1000.0 --alpha long --field 1 --delta 1.0
time python clrf_v2.py --spins 4096 --eqs 4096 --mcs 4096 --rseed 1234 --nTemp 120 --tMin 0.5 --tMax 1000.0 --alpha long --field 1 --delta 10.0

time python clrf_v2.py --spins 4096 --eqs 4096 --mcs 4096 --rseed 1234 --nTemp 120 --tMin 0.5 --tMax 6000.0 --alpha zero --field 0 --delta 0.1
time python clrf_v2.py --spins 4096 --eqs 4096 --mcs 4096 --rseed 1234 --nTemp 120 --tMin 0.5 --tMax 6000.0 --alpha zero --field 0 --delta 1.0
time python clrf_v2.py --spins 4096 --eqs 4096 --mcs 4096 --rseed 1234 --nTemp 120 --tMin 0.5 --tMax 6000.0 --alpha zero --field 0 --delta 10.0

time python clrf_v2.py --spins 4096 --eqs 4096 --mcs 4096 --rseed 1234 --nTemp 120 --tMin 0.5 --tMax 6000.0 --alpha zero --field 1 --delta 0.1
time python clrf_v2.py --spins 4096 --eqs 4096 --mcs 4096 --rseed 1234 --nTemp 120 --tMin 0.5 --tMax 6000.0 --alpha zero --field 1 --delta 1.0
time python clrf_v2.py --spins 4096 --eqs 4096 --mcs 4096 --rseed 1234 --nTemp 120 --tMin 0.5 --tMax 6000.0 --alpha zero --field 1 --delta 10.0



#time python clrf_v2.py --spins 4096 --eqs 4000 --mcs 4000 --rseed 1234 --nTemp 20 --tMin 0.5 --tMax 1500.0 --alpha longT --field 1 --delta 0.1
#time python clrf_v2.py --spins 4096 --eqs 4000 --mcs 4000 --rseed 1234 --nTemp 20 --tMin 0.5 --tMax 1500.0 --alpha longT --field 1 --delta 1.0
#time python clrf_v2.py --spins 4096 --eqs 4000 --mcs 4000 --rseed 1234 --nTemp 20 --tMin 0.5 --tMax 1500.0 --alpha longT --field 1 --delta 10.0
