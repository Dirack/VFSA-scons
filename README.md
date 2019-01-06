# VFSA-scons
## Projeto da tese de doutorado em geofísica

>O VFSA-scons é uma adaptação do algoritmo de otimização global Very Fast Simulated Aneeling (VFSA) para a inversão 
>dos parâmetros do CRS zero-offset. Escrito em linguagem C, utiliza API's e SConstruct para funcionar no 
>pacote de processamento sísmico MADAGASCAR.

Site oficial do MADAGASCAR: http://www.ahay.org/wiki/Main_Page

Estes dois programas (sfvfsa e sfcrs) realizam a inversão dos parâmetros do CRS zero offset (RN, RNIP e BETA) 
para uma aproximação de tempo de trânsito CRS escolhida. 
Para realizar tal procedimento, necessita de uma superfície de tempo de trânsito CRS extraída dos dados (modelada). 
A inversão é realizada a partir da coerência (semblance) entre uma superfície de tempo de trânsito aproximada e 
a superfície de tempo de trânsito modelada. Em resumo, obtem-se os parâmetros do CRS zero offset otimizados, 
quando a superfície de tempo de trânsito aproximada se ajusta da melhor maneira possível 
à superfície de tempo de trânsito modelada.

O procedimento descrito acima é realizado automaticamente com o SConstruct (script do MADAGASCAR), bastando a adição dos 
programas ao pacote MADAGASCAR já instalado na máquina do usuário.

### Resumo:

### Referências:

A convergência entre as superfícies de tempo de trânsito 
é atingida pelo algoritmo de otimização global Very Fast Simulated Aneeling (VFSA), descrito em detalhes em: 

INGBER, L. Very fast simulated re-annealing. Math1. Comput. Modelling, v. 12, p.967–973, 1989.
