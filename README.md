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
programas ao pacote MADAGASCAR já instalado na máquina do usuário. A vantagem do SConstruct é que ele pode ser readaptado
para outros problemas, retirando partes do script que funcionam de maneira independente, como a etapa de modelagem e
simulação da aquisição sísmica.

### Resumo:

Na etapa de modelagem, simulamos a aquisição de dados sísmicos a partir do modelo de um refletor gaussiano 
em um meio de variação linear de velocidade (a velocidade crece linearmente com a profundidade). A aquisição
é simulada diretamente no domínio meio-afastamento x CMP x tempo (um domínio 3D chamado genericamente de cubo de dados).

![](https://raw.githubusercontent.com/Dirack/Images/master/dome.jpeg)

![](https://raw.githubusercontent.com/Dirack/Images/master/data.jpeg)

Depois, extraímos a superfície de tempo de trânsito modelada do cubo de dados, obtida na etapa anterior, e a representamos 
em um mapa 2D.

![](https://raw.githubusercontent.com/Dirack/Images/master/pick.jpeg)

Finalmente, o algoritmo VFSA produz o melhor ajuste entre a superfície de tempo de trânsito modelada e uma superfície de
tempo de trânsito aproximada definida pelo usuário: Cada aproximação de tempo de trânsito CRS é uma equação que define
a superfície de tempo de trânsito CRS aproximada dados as coordenadas CMP e meio-afastamento (que sao as mesmas para a 
superfície modelada e aproximada) e três parâmetros do CRS (RN, RNIP e BETA), de modo que ao variar estes parâmetros,
a superfície CRS aproximada se torna semelhante ou difere da superfície modelada.

![](https://raw.githubusercontent.com/Dirack/Images/master/err-0.jpeg)

Assim, o algoritmo VFSA busca obter o melhor ajuste possível entre a superfície de tempo de trânsito modelada e a superfície
de tempo de trânsito aproximada, e isso implica em obter os parâmetros RN, RNIP e BETA que melhor produzem tal ajuste. O
critério de convergência e a formulação matemática do algoritmo e das aproximações são detalhados nas referências a seguir.

### Referências:

* A convergência entre as superfícies de tempo de trânsito 
é atingida pelo algoritmo de otimização global Very Fast Simulated Aneeling (VFSA), descrito em detalhes em: 

		INGBER, L. Very fast simulated re-annealing. Math1. Comput. Modelling, v. 12, p.967–973, 1989.

* Referências sobre as aproximações de tempo de trânsito CRS utilizadas:

		FOMEL, S.; KAZINNIK, R. Nonhyperbolic common reflection surface. Geophysical Prospecting, v. 61, p. 21–27, 2013.

		HöCHT, G. Traveltime approximations for 2D and 3D media and kinematic wavefield attributes. Tese (Doutorado) — Faculdade de Física Karlsruhe (TH) genehmigte, 2002.

		NEVES, R. Aproximações não hiperbólicas do tempo de trânsito utilizando aproximantes de Padé. Dissertação (Mestrado) — Universidade Federal do Pará - UFPa, Belém - PA, 2017.

		JAGER, R. et al. Common-reflection-surface stack: image and attributes. Geophysics,v. 66, p. 97–109, 2001.
