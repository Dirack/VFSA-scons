# VFSA-scons
## Projeto da tese de doutorado em geofísica

>O VFSA-scons é uma adaptação do algoritmo de otimização global Very Fast Simulated Aneeling (VFSA) para a inversão 
>dos parâmetros do CRS zero-offset. Escrito em linguagem C, utiliza API's e SConstruct para funcionar no 
>pacote de processamento sísmico MADAGASCAR.

>Site oficial do MADAGASCAR: http://www.ahay.org/wiki/Main_Page

Os programas _sfvfsa_ e _sfcrs_ realizam a inversão dos parâmetros do CRS zero offset (RN, RNIP e BETA) 
para uma aproximação de tempo de trânsito CRS escolhida. 
Tal procedimento necessita de uma superfície de tempo de trânsito CRS extraída dos dados (modelada). 
A inversão é realizada a partir da coerência (semblance) entre esta superfície de tempo de trânsito modelada 
e uma superfície de tempo de trânsito aproximada escolhida entre as seguintes: 

* CRS Hiperbólico 
* CRS NÃO hiperbólico
* Aproximação Parabólica do CRS quarta ordem
* Aproximação Hiperbólica do CRS quarta ordem
* Aproximação Deslocada do CRS quarta ordem
* Aproximação Parabólica CRS-Padé com expansão em h
* Aproximação Parabólica CRS-Padé com expansão em m
* Aproximação Hiperbólica CRS-Padé com expansão em h
* Aproximação Hiperbólica CRS-Padé com expansão em m
* Aproximação Deslocada CRS-Padé com expansão em h
* Aproximação Deslocada CRS-Padé com expansão em m

Em resumo, os parâmetros do CRS zero offset serão otimizados, 
quando a superfície de tempo de trânsito aproximada se ajusta da melhor maneira possível 
à superfície de tempo de trânsito modelada.

O procedimento descrito acima é realizado automaticamente com o SConstruct (script do MADAGASCAR), bastando a adição dos 
programas ao pacote MADAGASCAR já instalado na máquina do usuário.

### Resumo:

### Referências:

* A convergência entre as superfícies de tempo de trânsito 
é atingida pelo algoritmo de otimização global Very Fast Simulated Aneeling (VFSA), descrito em detalhes em: 

		INGBER, L. Very fast simulated re-annealing. Math1. Comput. Modelling, v. 12, p.967–973, 1989.

* Referências sobre as aproximações de tempo de trânsito CRS utilizadas:

		FOMEL, S.; KAZINNIK, R. Nonhyperbolic common reflection surface. Geophysical Prospecting, v. 61, p. 21–27, 2013.

		HöCHT, G. Traveltime approximations for 2D and 3D media and kinematic wavefield attributes. Tese (Doutorado) — Faculdade de Física Karlsruhe (TH) genehmigte, 2002.

		NEVES, R. Aproximações não hiperbólicas do tempo de trânsito utilizando aproximantes de Padé. Dissertação (Mestrado) — Universidade Federal do Pará - UFPa, Belém - PA, 2017.

		JAGER, R. et al. Common-reflection-surface stack: image and attributes. Geophysics,v. 66, p. 97–109, 2001.
