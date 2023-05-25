# Simulador de Trânsito Planetário (SimTransP)

O Simulador de Trânsito Planetário (SimTransP) é um projeto que implementa a classe `SimTransP`, a qual gera curvas de luz simuladas e permite a análise de suas características. As curvas de luz são geradas utilizando valores padrões típicos observados pela missão Kepler. Além disso, a classe também adiciona manchas estelares de acordo com os parâmetros de entrada descritos em [Bonomo e Lanza (2012)](https://www.aanda.org/articles/aa/abs/2012/11/aa19999-12/aa19999-12.html).

## Dependências

Para utilizar o `SimTransP`, é necessário ter as seguintes dependências instaladas:

- numpy: Biblioteca para cálculos numéricos.
- PyAstronomy: Biblioteca para cálculos astronômicos.
- matplotlib: Biblioteca para criação de gráficos e visualizações.
- seaborn: Biblioteca para melhoria estética.
- astropy: Biblioteca para cálculos astrofísicos e manipulação de dados.

Certifique-se de que essas dependências estão instaladas antes de utilizar o SimTransP (ou use `pip install -r requisitos.txt`).

## Uso

Um exemplo de como utilizar o `SimTransP` pode ser encontrado no notebook [Exemplo_trânsito_planetário](https://github.com/SarahBarbosa/SimTransP/blob/main/Exemplo_tr%C3%A2nsito_planet%C3%A1rio.ipynb). Nesse exemplo, são demonstrados os passos para gerar uma curva de luz simulada, adicionar ruído e manchas estelares, e realizar análises sobre a curva gerada.

## Conclusão

A classe `SimTransP` oferece uma maneira conveniente e simples de gerar curvas de luz com diversos recursos e realizar caracterizações. Ela permite adicionar pontos, ruídos e trânsitos, além de analisar o periodograma e avaliar o erro relativo do período observado. Sinta-se à vontade para explorar o SimTransP e adicionar novos recursos conforme necessário.

> Status do Projeto: Concluido :heavy_check_mark:
