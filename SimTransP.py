import numpy as np
from numpy import sqrt, arcsin, sin, cos, pi, matlib
from math import ceil

from scipy.fft import fft, ifft

from PyAstronomy.pyTiming.pyPeriod import Gls

import matplotlib.pyplot as plt
import seaborn as sns

class SimTransP:
    """
    Esta classe gera curva de luz e suas caracterizações
    """


    def __init__(self, tmin = 0, tmax = 90,cadencia = 0.020833, ai = 90):
        
        '''
        Os valores padrões correspondem a uma curva de luz típica observada 
        pela missão Kepler
        
        ---------------
        
        tmax = quarter da missão Kepler
        
        cadencia = cadência da missão Kepler
        
        ai = Ângulo de inclinação da estrela (entre 0 e 90 graus).
        '''

        #Adicionando atributos
        self.cadencia = cadencia
        self.ai = ai
        self.tmin = tmin
        self.tmax = tmax
        
        #Criando um vetor tempo
        self.tempo = np.arange(tmin,tmax,cadencia)

        #Vetor do fluxo de uma estrela sem spot e sem rotação
        self.fluxo = np.array([1.]*len(self.tempo))
        
    
    def plot(self,eixo_x=20,eixo_y = 8,jump=1, titulo = 'Curva de luz', **kwargs):
        
        '''
        Plota o gráfico
        
        ---------------------
        
        eixo_x = Largura em polegadas da figura
        
        eixo_y = altura em polegadas da figura
        
        jump = Remove pontos periódicos
        
        titulo = título do gráfico
        '''
        
        #Adicionando um estilo ao gráfico
        sns.set_style('darkgrid')
        sns.set_context(context = 'paper', font_scale= 2.5)
        
        plt.figure(figsize=(eixo_x,eixo_y))
        
        plt.plot(self.tempo[::jump],self.fluxo[::jump],'.', c = 'black', markersize=1, **kwargs)
        
        plt.xlabel('Tempo [dias]')
        plt.ylabel('Fluxo normalizado')

        plt.title(titulo)
        
    def add_spots(self, gamma1 = 0.4, gamma2 = 0.262, cs = 0.677, cf0 = 0.115, Q = 1.6, As = 0.01,
                 Peq = 15, AmpliP = 0.2):  
        
        '''
        Adicionando 2 spots de acordo com os parâmetros de entrada descritos no paper: BONOMO, Aldo S.; 
        LANZA, Antonino F. Starspot activity and rotation of the planet-hosting star Kepler-17. Astronomy 
        & Astrophysics, v. 547, p. A37, 2012.
        
        ------------------
        
        gamma1, gamma2 = parâmetros do escurecimento do limbo
        
        cs = coeficiente que especifica o contraste do spot bolométrico
        
        cf0 = contraste da fácula assumido como solar
        
        Q = a razão entre a área das fáculas e os pontos frios em regiões ativas
        
        As = Área da estrela (entre 1% e 2%)
        
        Peq = período de rotação da estrela (escolha para subsolar, solar e supersolar)
        
        AmpliP = medida da rotação diferencial (entre 0-0.5)
        '''
        
        self.Peq = Peq
        
        ap = 1 - gamma1 - gamma2
        bp = gamma1 + 2 * gamma2
        cp = - gamma2
        
        ii = self.ai * pi / 180     # ângulo de inclinação em respeito da linha de visada
        theta = 0.17 * pi           # é a posição da região ativa na estrela
        
        P2 = Peq * (1 - AmpliP * (sin(theta) ** 2)) ** (-1)
        lambda_ = 2 * pi / 5
        
        # Spot 1
        
        omega_1 = 2 * pi / Peq      # velocidade angular da estrela
        c = (ap + 2 * bp / 3 + cp / 2) ** (-1)
        
        # Ângulo do limbo:
        
        mii_1 = cos(ii) * cos(theta) + sin(ii) * sin(theta) * cos(lambda_ + omega_1 * np.linspace(0, self.tempo[-1], len(self.tempo))) 
        
        cf = cf0 * (1 - mii_1)      # Intensidade do contraste da facúla
        
        s1 = 1 + As * (ap + bp * mii_1 + cp * mii_1 ** 2) * c * (Q * cf - cs) * mii_1
        
        # Spot 2
        
        omega_2 = 2 * pi / P2
        
        mii_2 = cos(ii) * cos(theta) + sin(ii) * sin(theta) * cos(lambda_ + omega_2 * np.linspace(0, self.tempo[-1], len(self.tempo))) 
        
        cf = cf0 * (1 - mii_2)
        
        s2 = 1 + As * (ap + bp * mii_2 + cp * mii_2 ** 2) * c * (Q * cf - cs) * mii_2
        
        spots = (s1 + s2) - np.mean(s1 + s2)
        
        self.fluxo = self.fluxo * spots
    
    def add_noise(self, std = 0.0006):
        
        '''
        Adiciona ruído na curva de luz
        
        ------------------
        
        std = desvio padrão
        '''
        
        self.fluxo_noise = self.fluxo[:]
        self.fluxo_noise = np.random.normal(self.fluxo_noise,std)
        self.fluxo = self.fluxo_noise    

    
    def add_transito(self, Re = 1, Rp = 0.1, Porb = 15, b = 0):
        
        '''
        Adiciona trânsito planetário na curva de luz da estrela.
        
        ------------------
        
        Re = Raio da estrela em unidade solares
        
        Rp = Fração de raios solares. Ex: Raio da Terra = 0.00915
        
        Porb = Período orbital em dias
        
        b = parâmetro de impacto
        '''
        
        # Constantes
        
        from astropy.constants import M_sun, R_sun, G
        import astropy.units as u
        
        M_sun = M_sun.value
        R_sun = R_sun.value
        Gday = G.to(u.m**3 / (u.kg * u.d**2))
        Gday = Gday.value
        
        # Parâmetros de entrada atribuidos
        
        Rstar = Re * R_sun
        Rplanet = Rp * R_sun
        depth = (Rplanet/Rstar)**2
        
        # Semi-eixo maior (aa)
        
        aa = (Gday * M_sun * Porb ** 2/ (4 * pi **2))**(1/3)
          
        # Tempo total
        
        part1 = Rstar / aa
        part2 = (1 + Rplanet/Rstar)**2 - b**2
        part3 = (1 - Rplanet/Rstar)**2 - b**2
        part4 = 1 - (cos(self.ai * pi / 180))**2
        
        tTotal = (Porb / pi) * arcsin(part1 * sqrt(part2/part4)) # em dias
        
        tF = (Porb / pi) * arcsin(sin(tTotal * (pi / Porb))) * sqrt(part3)/sqrt(part2) # em dias
        
        tramp = (tTotal - tF)/2                         # Tempo de rampa
        
        ntTF = tramp/self.cadencia               
        nt = ceil(ntTF)                                 # Número de pontos entre o 1º e o 2º contato 
        ntF = ceil(tF/self.cadencia)                    # Número de pontos entre o 2º e o 3º contato
        
        interTransit = ceil(Porb/self.cadencia - tramp) # Número de pontos entre os dois trânsitos
        
        def linspace_float(low, up, leng):
            step = ((up-low) * 1.0 / leng)
            return [low+i*step for i in range(leng)]
    
        rampdow = linspace_float(1, 1-depth, nt)
        rampup = linspace_float(1-depth,1,nt)
        
        x1 = (1 - depth) * np.ones(ntF)
        x2 = np.ones(interTransit)
        
        pulse = np.concatenate([rampdow,x1,rampup,x2])
        rep = int(len(self.fluxo) / len(pulse))
        pulse_train = np.matlib.repmat(pulse, 1, rep)
        
        PT = pulse_train[0]
        
        fill = len(self.fluxo) - len(PT)
        PT = np.append(PT, [1]*fill)
        
        self.fluxo += PT  
    
    def periodograma(self):
        
        '''
        Recupera o período de rotação pelo peridograma de Lomb-Scargle generalizado. 
        A função GLS fornece uma implementação do periodograma Generalized Lomb-Scargle 
        conforme descrito por echmeister & Kuerster 2009 (A&A 496, 577) e anteriormente 
        por Ferraz-Mello 1981 (AJ 86, 691). 
        '''
        
        gls = Gls((self.tempo, self.fluxo))
        
        # Obtenha o índice associado à maior potência
        ifmax = np.argmax(gls.power)
        
        # e maior potência e frequência associada
        fmax = gls.freq[ifmax]
        
        # Converter frequência em período
        self.hpp = 1./fmax
        print(f"Período de maior potência Lomb-Scargle: {self.hpp} dias")
        
        sns.set_context(context = 'paper', font_scale= 1.5)
        
        plt.plot(1./gls.freq, gls.power, '.k')
        
        plt.xlim(self.tmin, self.tmax)
        
        plt.xlabel('Tempo (dias)')
        plt.ylabel('Potência Lomb-Scargle')
    
    def erro_relativo(self):
        
        '''
        Retorna o erro relativo os períodos teórico e observado pelo peridograma de 
        Lomb-Scargle generalizado.
        '''
        erro_relativo = abs(self.Peq - self.hpp) / self.Peq
        
        erro_relativo_per = erro_relativo * 100
        
        print(f'Erro relativo: {erro_relativo} ou {erro_relativo_per}%')