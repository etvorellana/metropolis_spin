#!/bin/bash

# Run the application clr

# Experimento 1. Variando a quantidade de eqs e mcs da mesma forma em cada teste, para a mesma quantidade de spins
time python clr_v1.py --spins 2048 --eqs 2000 --mcs 2000 --rseed 1234 --nTemp 120 --tMin 0.5 --tMax 20.0
time python clr_v1.py --spins 4096 --eqs 4000 --mcs 4000 --rseed 1234 --nTemp 120 --tMin 0.5 --tMax 20.0
#time python clr_v1.py --spins 1024 --eqs 4000 --mcs 4000 --rseed 1234 --nTemp 120 --tMin 0.5 --tMax 20.0

# Experimento 2. Todos os testes com eqs=1000, variando mcs, para a mesma quantidade de spins
#time python clr_v1.py --spins 1024 --eqs 1000 --mcs 1000 --rseed 1234 --nTemp 120 --tMin 0.5 --tMax 20.0
#time python clr_v1.py --spins 1024 --eqs 1000 --mcs 2000 --rseed 1234 --nTemp 120 --tMin 0.5 --tMax 20.0
#time python clr_v1.py --spins 1024 --eqs 1000 --mcs 4000 --rseed 1234 --nTemp 120 --tMin 0.5 --tMax 20.0

# Run the application clrf

# Experimento 1. Variando a quantidade de eqs e mcs da mesma forma em cada teste, para a mesma quantidade de spins
time python clrf_v1.py --spins 2048 --eqs 2000 --mcs 2000 --rseed 1234 --nTemp 120 --tMin 0.5 --tMax 20.0 --delta 0.5
time python clrf_v1.py --spins 4096 --eqs 4000 --mcs 4000 --rseed 1234 --nTemp 120 --tMin 0.5 --tMax 20.0 --delta 0.5
#time python clrf_v1.py --spins 1024 --eqs 4000 --mcs 4000 --rseed 1234 --nTemp 120 --tMin 0.5 --tMax 20.0 --delta 0.5

# Experimento 2. Todos os testes com eqs=1000, variando mcs, para a mesma quantidade de spins
#time python clrf_v1.py --spins 1024 --eqs 1000 --mcs 1000 --rseed 1234 --nTemp 120 --tMin 0.5 --tMax 20.0 --delta 0.5
#time python clrf_v1.py --spins 1024 --eqs 1000 --mcs 2000 --rseed 1234 --nTemp 120 --tMin 0.5 --tMax 20.0 --delta 0.5
#time python clrf_v1.py --spins 1024 --eqs 1000 --mcs 4000 --rseed 1234 --nTemp 120 --tMin 0.5 --tMax 20.0 --delta 0.5
