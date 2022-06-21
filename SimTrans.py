# Importando bibliotecas

import numpy as np
from numpy import sqrt, arcsin, sin, cos, pi, matlib
from math import ceil

from astropy.constants import M_sun, R_sun, G
import astropy.units as u

import matplotlib.pyplot as plt
import seaborn as sns

from stochastic.processes.noise import FractionalGaussianNoise   # pip install stochastic

# Funções criadas

def linspace_float(low, up, leng):
    step = ((up-low) * 1.0 / leng)
    return [low+i*step for i in range(leng)]

# Estilo dos gráficos

sns.set_style('darkgrid')
sns.set_context(context = 'paper', font_scale= 1.5)

# Constantes

M_sun = M_sun.value
R_sun = R_sun.value

# Conversão de unidades

Gday = G.to(u.m**3 / (u.kg * u.d**2))
Gday = Gday.value

# Parâmetros de entrada

Rstar = 1 * R_sun               # Raio da estrela em unidades solares
Rplanet = 0.1 * R_sun           # Fração de raios solares. Raio da Terra = 0.00915 Rsun
depth = (Rplanet/Rstar)**2      # Profundidade de trânsito
Porb = 15                       # Período orbital em dias
ai = 90                         # Ângulo de inclinação da estrela (entre 0 e 90 graus)
cadence_kepler = 0.020833       # Cadência em dias (Missão Kepler)
cadence_plato = 0.0002895       # Cadência em dias (Missão Plato)

# Escolha da cadência

cadence_question = str(input(' Cadência utilizada, PLATO ou Kepler: '))

if cadence_question.lower() in ['plato']:
    cadence = cadence_plato
else:
    cadence = cadence_kepler
    

# Semi-eixo maior (aa)

aa = (Gday * M_sun * Porb ** 2/ (4 * pi **2))**(1/3)

# Parâmetro de impacto, em a/Rstar

b = (aa / Rstar) * cos(ai * pi / 180)

# Tempo total

part1 = Rstar / aa
part2 = (1 + Rplanet/Rstar)**2 - b**2
part3 = (1 - Rplanet/Rstar)**2 - b**2
part4 = 1 - (cos(ai * pi / 180))**2

tTotal = (Porb / pi) * arcsin(part1 * sqrt(part2/part4)) # em dias

tT = tTotal * 24 # em horas

tF = (Porb / pi) * arcsin(sin(tTotal * (pi / Porb))) * sqrt(part3)/sqrt(part2) # em dias

tFh = tF * 24 # em horas

print(f'\n duração do trânsito total: {tT:.4f} hr')
print(f'\n duração entre o segundo e terceiro contato: {tFh:.4f} hr')
print(f'\n diferença entre as durações: {(tT - tFh):.4f} hr')

# Criando um trânsito planetário (onda trapezoidal)

tramp = (tTotal - tF)/2                    # Tempo de rampa
ntTF = tramp/cadence                
nt = ceil(ntTF)                            # Número de pontos entre o 1º e o 2º contato 
ntF = ceil(tF/cadence)                     # Número de pontos entre o 2º e o 3º contato

interTransit = ceil(Porb/cadence - tramp)  # Número de pontos entre os dois trânsitos

rampdow = linspace_float(1, 1-depth, nt)
rampup = linspace_float(1-depth,1,nt)

x1 = (1 - depth) * np.ones(ntF)
x2 = np.ones(interTransit)

rep = 10
pulse = np.concatenate([rampdow,x1,rampup,x2])
pulse_train = np.matlib.repmat(pulse, 1, rep)

PT = pulse_train[0]
time = [float(x)*cadence for x in range(0,len(PT))]


fig, ax = plt.subplots(figsize = (20,10))
plt.plot(time, PT, c = 'black')

plt.xlabel('tempo (em dias baricêntricos Julianos)')
plt.ylabel('Fluxo normalizado')
plt.title('Trânsito planetário sem ruído')

# Criação do ruído 

fs = len(time) - 1
fGn = FractionalGaussianNoise(hurst = 0.5)
noise = fGn.sample(n = fs + 1)

nPT = PT + noise

# Plot do ruído

fig, ax = plt.subplots(figsize = (20,10))

plt.scatter(time,nPT, c = 'black', s = 1)
plt.xlabel('tempo (em dias baricêntricos Julianos)')
plt.ylabel('Fluxo normalizado')
plt.title('Trânsito planetário com ruído')

### Modelo da modulação rotacional ###

# Parâmetros de entrada

# Paper base: BONOMO, Aldo S.; LANZA, Antonino F. Starspot activity and rotation 
# of the planet-hosting star Kepler-17. Astronomy & Astrophysics, v. 547, p. A37, 2012.

gamma1 = 0.3985             # gamma 1 e 2 são parâmetros do escurecimento do limbo
gamma2 = 0.2586  
           
ap = 1 - gamma1 - gamma2    # ap, bp e cp são coeficientes que nós dá a magnitude
bp = gamma1 + 2 * gamma2      # bolométrica da fotosfera como uma função do ângulo
cp = - gamma2                # do limbo adotado como uma lei quadrátca do escurecimento
                            # do limbo

cs = 0.850                  # coeficiente que especifica o contraste do spot bolométrico
cf0 = 1.115                 
Q = 0.5                     # é uma função da fase do ciclo solar
As = 0.01                   # é a área do spot da região ativa em unidades de ...

ii = ai * pi / 180            # ângulo de inclinação em reseito da linha de visada
theta = 0.17 * pi           # é a posição da região ativa na estrela 
             
Peq = 15                    # em dias
AmpliP = 0.2                # alpa = DeltaPer/Per

P2 = Peq * (1 - AmpliP * (sin(theta) ** 2)) ** (-1)
lambda_ = 2 * pi / 5


# Spot 1

omega_1 = 2 * pi / Peq        # velocidade angular da estrela

c = (ap + 2 * bp / 3 + cp / 2) ** (-1)

# Ângulo do limbo:
mii_1 = cos(ii) * cos(theta) + sin(ii) * sin(theta) * cos(lambda_ + omega_1 * np.linspace(0, time[-1], len(time))) 

cf = cf0 * (1 - mii_1)        # Intensidade do contraste da facúla

s1 = 1 + As * (ap + bp * mii_1 + cp * mii_1 ** 2) * c * (Q * cf - cs) * mii_1

# Spot 2

omega_2 = 2 * pi / P2

mii_2 = cos(ii) * cos(theta) + sin(ii) * sin(theta) * cos(lambda_ + omega_2 * np.linspace(0, time[-1], len(time))) 

cf = cf0 * (1 - mii_2)

s2 = 1 + As * (ap + bp * mii_2 + cp * mii_2 ** 2) * c * (Q * cf - cs) * mii_2

spots = (s1 + s2) - np.mean(s1 + s2)

# Plot

wts = spots + nPT

fig, ax = plt.subplots(figsize = (20,10))
plt.scatter(time, wts, c = 'black', s = 1)
plt.xlabel('tempo (dias baricêntricos Julianos)')
plt.ylabel('Fluxo normalizado')
plt.title('Trânsito planetário com ruído e rotação')