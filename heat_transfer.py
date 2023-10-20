#heat_tranfer.py

import numpy as np

rc  = (29+.75)/2 /1000 # m - raio externo do carbono
rpi = .250 # m - raio interno do pvc
rpe = rpi+.003 # m - raio externo do pvc
L   = 1.5 # m - comprimento do tubo

P = 3.18 # W - potência esquentando o carbono
Tinf = 25 # °C - temperatura ambiente

ka = .03 # W/mK - condutividade térmica do ar
kp = .2 # W/mK - condutividade térmica do pvc
h = 1 # W/m²K - coeficiente de convecção do pvc com o ar externo

sigma = 5.67e-8 # W/m²K⁴ - constante de Stefan-Boltzmann

dra = .0001 # m - espaçamento dos pontos no ar interno
drp = .00001 # m - espaçamento dos pontos no pvc

# Variáveis auxiliares
Na = int(np.round((rpi-rc )/dra))
Np = int(np.round((rpe-rpi)/drp))
N = Na+Np+1
print(Na)
print(Np)
print(N)

Ac  = 2*np.pi*rc*L
Api = 2*np.pi*rpi*L

Tc = Tinf
Tp = Tinf
PR = 0

# A*T = B
T = np.zeros((N))
A = np.zeros((N,N))
B = np.zeros((N))

Tc_ant = Tc+1

cont = 0
while abs(Tc-Tc_ant)>1e-3 and cont<100:

  # Calcula a matriz A e o vetor B
  for i in range(N):
    if i==0:
      A[i][i  ] = 1
      A[i][i+1] = -1
      B[i] = dra/Ac/ka*(P-PR)
    elif 0<i<Na:
      r = rc + dra*i
      A[i][i-1] = 2*r - dra
      A[i][i  ] = -4*r
      A[i][i+1] = 2*r + dra
      B[i] = 0
    elif i==Na:
      A[i][i-1] = ka/dra
      A[i][i  ] = -ka/dra-kp/drp
      A[i][i+1] = kp/drp
      B[i] = -PR/Api
    elif Na<i<N-1:
      r = rpi + drp*(i-Na+1)
      A[i][i-1] = 2*r - drp
      A[i][i  ] = -4*r
      A[i][i+1] = 2*r + drp
      B[i] = 0
    elif i==N-1:
      A[i][i-1] = -1
      A[i][i  ] = h*drp/kp + 1
      B[i] = h*drp/kp * Tinf
  # Resolve o sistema de equações
  T = np.linalg.solve(A, B)

  Tc = T[0]
  Tp = T[Na-1]
  PR = sigma*(Ac*Tc**4 - Api*Tp**4)

  cont+=1

print("P  = "+str(P)+" W")
print("PR = "+str(PR)+" W")
