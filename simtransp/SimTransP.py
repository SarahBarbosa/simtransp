from PyAstronomy.pyTiming.pyPeriod import Gls
from astropy.constants import M_sun, R_sun, G
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import numpy.matlib as matlib
import math

class STP:
    """
    Esta classe gera curva de luz e suas caracterizações.
    """

    def __init__(self, tmin: float = 0, tmax: float = 90, cadencia: float = 0.020833, ai: float = 90) -> None:
        """
        Inicializa a classe STP com os parâmetros especificados.

        Parameters
        ----------
        tmin : float, optional
            Tempo mínimo em dias, por padrão 0.
        tmax : float, optional
            Tempo máximo em dias, por padrão 90.
        cadencia : float, optional
            Cadência da missão Kepler em dias, por padrão 0.020833.
        ai : float, optional
            Ângulo de inclinação da estrela em graus (entre 0 e 90), por padrão 90.
        """

        #Adicionando atributos
        self.cadencia = cadencia
        self.ai = ai
        self.tmin = tmin
        self.tmax = tmax
        
        #Criando um vetor tempo
        self.tempo = np.arange(tmin,tmax,cadencia)

        #Vetor do fluxo de uma estrela sem spot e sem rotação
        self.fluxo = np.array([1.]*len(self.tempo))
    
        # Configurações padrão do matplotlib
        self._configurar_matplotlib()
        
    def _configurar_matplotlib(self) -> None:
        """
        Configura os parâmetros do matplotlib.
        """
        plt.rcParams.update({
            "axes.spines.right": False,
            "axes.spines.top": False,
            "font.size": 12,
            "axes.labelsize": 12,
            "axes.titlesize": 12,
            "legend.fontsize": 10,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "figure.figsize": (10, 4)
        })
        
    
    def plot(self, titulo: str = 'CURVA DE LUZ', **kwargs) -> plt.axes:
        """
        Plota o gráfico da curva de luz.

        Parameters
        ----------
        titulo : str, optional
            Título do gráfico, por padrão 'CURVA DE LUZ'.
        
        Returns
        -------
        plt.axes
            O objeto de eixo do gráfico.
        """    
        ax = plt.gca()
        ax.scatter(self.tempo, self.fluxo, c='k', s = 1, **kwargs)
        
        ax.set_xlabel('Tempo [dias]')
        ax.set_ylabel('Fluxo normalizado')

        ax.set_title(titulo)
        
        return ax
        
    def add_spots(self, gamma1: float = 0.4, gamma2: float = 0.262, cs: float = 0.677, cf0: float = 0.115,
                  Q: float = 1.6, As: float = 0.01, Peq: float = 15, AmpliP: float = 0.2) -> None:
        
        """
        Adicionando 2 spots de acordo com os parâmetros de entrada descritos no paper: BONOMO, Aldo S.; 
        LANZA, Antonino F. Starspot activity and rotation of the planet-hosting star Kepler-17. Astronomy 
        & Astrophysics, v. 547, p. A37, 2012.
        
        Parameters
        ----------
        gamma1 : float, optional
            Parâmetro do escurecimento do limbo, por padrão 0.4.
        gamma2 : float, optional
            Parâmetro do escurecimento do limbo, por padrão 0.262.
        cs : float, optional
            Coeficiente que especifica o contraste do spot bolométrico, por padrão 0.677.
        cf0 : float, optional
            Contraste da fácula assumido como solar, por padrão 0.115.
        Q : float, optional
            A razão entre a área das fáculas e os pontos frios em regiões ativas, por padrão 1.6.
        As : float, optional
            Área da estrela em percentual (entre 1% e 2%), por padrão 0.01.
        Peq : float, optional
            Período de rotação da estrela em dias, por padrão 15.
        AmpliP : float, optional
            Medida da rotação diferencial (entre 0 e 0.5), por padrão 0.2.
        """
        self.Peq = Peq
        ap = 1 - gamma1 - gamma2
        bp = gamma1 + 2 * gamma2
        cp = - gamma2
        
        ii = self.ai * np.pi / 180     # ângulo de inclinação em respeito da linha de visada
        theta = 0.17 * np.pi           # é a posição da região ativa na estrela
        
        P2 = Peq * (1 - AmpliP * (np.sin(theta) ** 2)) ** (-1)
        lambda_ = 2 * np.pi / 5
        
        ## Spot 1
        
        omega_1 = 2 * np.pi / Peq      # velocidade angular da estrela
        c = (ap + 2 * bp / 3 + cp / 2) ** (-1)
        
        ## Ângulo do limbo
        
        mii_1 = np.cos(ii) * np.cos(theta) + \
            np.sin(ii) * np.sin(theta) * np.cos(lambda_ + omega_1 * np.linspace(0, self.tempo[-1], len(self.tempo))) 
        
        cf = cf0 * (1 - mii_1)      # Intensidade do contraste da facúla
        
        s1 = 1 + As * (ap + bp * mii_1 + cp * mii_1 ** 2) * c * (Q * cf - cs) * mii_1
        
        ## Spot 2
        
        omega_2 = 2 * np.pi / P2
        
        mii_2 = np.cos(ii) * np.cos(theta) + \
            np.sin(ii) * np.sin(theta) * np.cos(lambda_ + omega_2 * np.linspace(0, self.tempo[-1], len(self.tempo))) 
        
        cf = cf0 * (1 - mii_2)
        
        s2 = 1 + As * (ap + bp * mii_2 + cp * mii_2 ** 2) * c * (Q * cf - cs) * mii_2
        
        spots = (s1 + s2) - np.mean(s1 + s2)
        
        self.fluxo = self.fluxo * spots
    
    def add_noise(self, std: float = 0.0006) -> None:
        """
        Adiciona ruído na curva de luz.

        Parameters
        ----------
        std : float, optional
            Desvio padrão do ruído, por padrão 0.0006.
        """        
        self.fluxo_noise = self.fluxo[:]
        self.fluxo_noise = np.random.normal(self.fluxo_noise,std)
        self.fluxo = self.fluxo_noise    
    
    def _linspace_float(self, low, up, leng):
        step = ((up-low) * 1.0 / leng)
        return [low+i*step for i in range(leng)]

    def add_transito(self, Re: float = 1, Rp: float = 0.1, Porb: float = 15, b: float = 0) -> None:
        """
        Adiciona trânsito planetário na curva de luz da estrela.

        Parameters
        ----------
        Re : float, optional
            Raio da estrela em unidades solares, por padrão 1.
        Rp : float, optional
            Fração de raios solares, por padrão 0.1.
        Porb : float, optional
            Período orbital em dias, por padrão 15.
        b : float, optional
            Parâmetro de impacto, por padrão 0.
        """
        # Constantes
        m_sun = M_sun.value
        r_sun = R_sun.value
        Gday = G.to(u.m**3 / (u.kg * u.d**2))
        Gday = Gday.value
        
        # Parâmetros de entrada atribuidos
        Rstar = Re * r_sun
        Rplanet = Rp * r_sun
        depth = (Rplanet/Rstar)**2
        
        # Semi-eixo maior (aa)
        aa = (Gday * m_sun * Porb ** 2/ (4 * np.pi **2))**(1/3)
          
        # Tempo total
        part1 = Rstar / aa
        part2 = (1 + Rplanet/Rstar)**2 - b**2
        part3 = (1 - Rplanet/Rstar)**2 - b**2
        part4 = 1 - (np.cos(self.ai * np.pi / 180))**2
        
        tTotal = (Porb / np.pi) * np.arcsin(part1 * np.sqrt(part2/part4))                                # em dias
        tF = (Porb / np.pi) * np.arcsin(np.sin(tTotal * (np.pi / Porb))) * np.sqrt(part3)/np.sqrt(part2) # em dias
        
        tramp = (tTotal - tF)/2                              # Tempo de rampa
        
        ntTF = tramp/self.cadencia               
        nt = math.ceil(ntTF)                                 # Número de pontos entre o 1º e o 2º contato 
        ntF = math.ceil(tF/self.cadencia)                    # Número de pontos entre o 2º e o 3º contato
        
        interTransit = math.ceil(Porb/self.cadencia - tramp) # Número de pontos entre os dois trânsitos
    
        rampdow = self._linspace_float(1, 1-depth, nt)
        rampup = self._linspace_float(1-depth,1,nt)
        
        x1 = (1 - depth) * np.ones(ntF)
        x2 = np.ones(interTransit)
        
        pulse = np.concatenate([rampdow,x1,rampup,x2])
        rep = int(len(self.fluxo) / len(pulse))
        pulse_train = matlib.repmat(pulse, 1, rep)
        
        PT = pulse_train[0]
        
        fill = len(self.fluxo) - len(PT)
        PT = np.append(PT, [1]*fill)
        
        self.fluxo += PT  
    
    @property
    def periodograma(self) -> float:
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

        return self.hpp
    
    @property
    def erro_relativo(self) -> float:
        """
        Retorna o erro relativo entre os períodos teórico e observado pelo periodograma de Lomb-Scargle generalizado.
        
        Returns
        -------
        float
            O erro relativo entre os períodos.
        """
        erro_relativo = abs(self.Peq - self.hpp) / self.Peq
        return erro_relativo
