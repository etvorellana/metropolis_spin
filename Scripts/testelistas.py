lista = []

for i in range(5):
    subLista = []
    for j in range(i,2*i+1):
        subLista.append(j)
    lista.append(subLista)
for subLista in lista:
    print(subLista)

print(5*'-')

novaLista = []
for i in range(len(lista)):
    for j in range(len(lista[i])):
        if len(novaLista) < j+1:
            novaLista.append([])
        novaLista[j].append(lista[i][j])

for subLista in novaLista:
    print(subLista)

import numpy as np

matrix = np.array(lista)
print (matrix)



 