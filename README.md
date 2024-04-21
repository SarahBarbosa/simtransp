# Simulador de Trânsito Planetário (SimTransP)

O Simulador de Trânsito Planetário (SimTransP) é um projeto desenvolvido na disciplina de Astrofísica Observacional (do curso de pós-graduação em Física da Universidade Federal do Ceará, Brasil), que gera curvas de luz simuladas e permite a análise de suas características de forma simples. As curvas de luz são geradas utilizando valores padrões típicos observados pela missão Kepler e adicionando duas manchas estelares de acordo com os parâmetros de entrada descritos em [Bonomo e Lanza (2012)](https://www.aanda.org/articles/aa/abs/2012/11/aa19999-12/aa19999-12.html).

## Instalação com Git

```bash
git clone https://github.com/SarahBarbosa/simtransp.git
cd simtransp
pip install -e .
```

## Uso

Um exemplo de como utilizar o `SimTransP` pode ser encontrado no notebook [Exemplo.ipynb](https://github.com/SarahBarbosa/SimTransP/blob/main/Exemplo_tr%C3%A2nsito_planet%C3%A1rio.ipynb). Neste exemplo, são demonstrados os passos para gerar uma curva de luz simulada, adicionar ruído e manchas estelares, e realizar análises sobre a curva gerada.

**Observação:** Esta classe não deve ser usada para fins científicos, sendo considerada apenas para propósitos educacionais na construção de curvas de luz simuladas.

> Status do Projeto: Concluído :heavy_check_mark: