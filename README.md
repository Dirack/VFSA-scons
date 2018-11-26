# VFSA-scons
Projeto VFSA adaptado para o pacote de processamento sísmico MADAGASCAR

Site do MADAGASCAR: http://www.ahay.org/wiki/Main_Page

Estes dois programas realizam a inversão dos parâmetros do CRS zero offset (RN, RNIP e BETA) para uma aproximação de tempo de trânsito CRS escolhida. Para realizar tal procedimento, necessita de uma superfície de tempo de trânsito CRS extraída dos dados (modelada). A inversão é realizada a partir da coerência (semblance) entre uma superfície de tempo de trânsito aproximada e a superfície de tempo de trânsito modelada. Em resumo, obtem-se os parâmetros do CRS zero offset otimizados, quando a superfície de tempo de trânsito aproximada se ajusta da melhor maneira possível a superfície de tempo de trânsito modelada.

A convergência entre as superfícies de tempo de trânsito é atingida por um algoritmo de otimização global, chamado Very Fast Simulated Aneeling (VFSA), descrito em detalhes em: INGBER, L. Very fast simulated re-annealing. Math1. Comput. Modelling, v. 12, p.967–973, 1989.
