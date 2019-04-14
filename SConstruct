#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Madagascar-VFSA-SEMB  (Script do Madagascar)
#
# Objetivo: Realizar a inversão dos parâmetros do CRS zero offset (RN, RNIP e BETA) a 
# partir do algoritmo Very Fast Simulated Aneeling (VFSA). A inversão utiliza como 
# critério de convergência o semblance entre a superfície modelada, extraída dos dados,
# e a superfície aproximada (escolhida uma aproximação de tempo de trânsito CRS). 
#
# Site: http://www.dirackslounge.online
# 
# Versão 1.1 - temp0 e c0 são fornecidos pelo usuário
#
# Programador: Rodolfo A. C. Neves 17/11/2018
#
# email: rodolfo_profissional@hotmail.com
#
# Funcionamento: digite 'scons' no terminal e aperte enter.
#
# Resumo e referências:
#
# 1-Modelagem:
#		Utilizamos como base o scons desenvolvido por Sergey Fomel e Roman Kazinik em (2013)
#		que realiza a modelagem Kirchoff sobre um modelo de um refletor gaussiano para simular
#		a aquisição de dados sísmicos em uma família CMP. Este script também faz o picking e
#		obtém a superfície de tempo de trânsito extraída dos dados que será útil nas próximas
#		etapas.
# Referências:	
#		http://www.reproducibility.org/RSF/book/tccs/crs/paper_html/paper.html (Fomel e Kazinik, 2013)
#		http://www.reproducibility.org/RSF/book/tccs/crs/dome2.html (scons)
#
# 2-Picking:
#		O script da etapa anterior também realiza o 'picking' da superfície de tempo de trânsito
#		modelada (extrai a superfície de tempo de trânsito t(m,h) a partir dos dados modelados).
#
# 3-VFSA:
#		Lê a superfície de tempo de trânsito CRS obtida na etapa de modelagem, e realiza a
#		Inversão dos parâmetros do CRS zero offset (RN, RNIP, BETA) através do algoritmo
#		Very fast simulated Aneeling (VFSA). Utilizando como critério de convergência, a 
#		coerência (Semblance) entre a superfície de tempo de trânsito modelada e uma superfície
#		de tempo de trânsito aproximada, escolhida a aproximação CRS.
#
#		A escolha da aproximação se dará através da determinação da variável de configuração 'app'
#		que representa uma aproximação de tempo de trânsito CRS, dada por:
#				-fomel [ app=1 ]
#		 		-Jager [ app=2 ]
#				-Germam-t [ app=3 ]
#				-Germam-t2 [ app=4 ]
#				-Germam-tshift [ app=5 ]
#				-Padé-t-h [ app=6 ]
#				-Padé-t-m [ app=7 ]
#				-Padé-t2-h [ app=8 ]
#				-Padé-t2-m [ app=9 ]
#				-Padé-tshift-h [ app=10 ]
#				-Padé-tshift-m [ app=11 ]
# 
# Referências sobre as aproximações de tempo de trânsito CRS utilizadas:
#
#	-FOMEL, S.; KAZINNIK, R. Nonhyperbolic common reflection surface. Geophysical Prospecting, v. 61, 
#	p. 21–27, 2013.
#
#	-HöCHT, G. Traveltime approximations for 2D and 3D media and kinematic wavefield attributes. 
#	Tese (Doutorado) — Faculdade de Física Karlsruhe (TH) genehmigte, 2002.
#
#	-NEVES, R. Aproximações não hiperbólicas do tempo de trânsito utilizando aproximantes de Padé. 
#	Dissertação (Mestrado) — Universidade Federal do Pará - UFPa, Belém - PA, 2017.
#
#	-JAGER, R. et al. Common-reflection-surface stack: image and attributes. Geophysics,v. 66, p. 97–109, 2001.
#
# Referências sobre o algoritmo VFSA:
#
#	-INGBER, L. Very fast simulated re-annealing. Math1. Comput. Modelling, v. 12, p.967–973, 1989.
#
# Licensa:
#  		Copyright (C) 2018 grupo de programação aplicada à geofísica (GPGEOF)
#  		da Universidade Federal do Pará (UFPA); Belém, Pará, Brasil. 
#
#  		Esse programa é um software livre; você pode redistribuir e/ou modificar
#  		sobre os termos da licensa pública geral (LPG) publicada pela Free 
#  		Software Foundation; na versão 2 da licensa, ou (a seu critério) qualquer
#  		versão posterior.
#
#  		Este programa é distribuído na esperança que será útil, mas SEM NENHUMA
#  		GARANTIA; nem mesmo a garantia implícita de MERCANTILIDADE ou SERVENTIA
#  		A UM PROPÒSITO DETERMINADO. veja a LPG licensa pública geral para mais
#  		detalhes.
#
#  		Você deve ter recebido uma cópia da LPG licensa pública geral junto com
#  		este programa; se não, escreva para a Free Software Foundation, Inc., 
#  		51 Franklin Street, Quinquagésimo andar, Boston, MA  02110-1301, USA.


# Importe a Biblioteca do madagascar 
from rsf.proj import *

# Biblioteca matemática do python
import math 

#------------------------{ Modelagem }-----------------------
# Cria um refletor gaussiano e o modelo de velocidades v(z)
# em que a velocidade cresce linearmente com a profundidade.
# Também cria o arquivo dip.rsf que é a inclinação local do
# refletor (necessária para a modelagem kirchoff).
#------------------------------------------------------------

# Criar o refletor a partir de uma função gaussiana
Flow('dome',None,
     '''
     math d1=0.01 n1=2001 o1=-5 unit1=km label1=Afastamento
     output="4-3*exp(-(x1-5)^2/9)"
     ''')

# TODO Comentários
for g in range(2):
    dome = 'dome%d' % g
    Plot(dome,'dome',
         '''
         graph min2=0 max2=4 min1=0 max1=10
         yreverse=y plotcol=%d plotfat=%d
         wantaxis=n wanttitle=n scalebar=y pad=n
         ''' % ((7,0)[g],(7,3)[g]))

#--------{ Modelo de velocidades } ---------
# Montar um modelo de velocidades v(z), 
# velocidade cresce com a profundidade
# v_0 é 1.5 (Km/s)
# 0.5 Gradiente de velocidade
#-------------------------------------------
Flow('vel','dome',
     '''
     window min1=0 max1=10 |
     spray axis=1 n=451 d=0.01 o=0 label=Profundidade unit=km |
     math output="1.5+0.5*x1+0.0*x2"
     ''')

Plot('vel',
     '''
     grey color=j allpos=y bias=1.5 scalebar=y wanttitle=n
     barreverse=y barlabel=Velocidade barunit=km/s
     ''')

Result('dome','vel dome0 dome1','Overlay')

# dip.rsf é a inclinação do refletor gaussiano (derivada)
Flow('dip','dome','math output="2/3*(x1-5)*input" ')

#---------------{ Modelagem Kirchoff }--------------------
# Unidades de medida (Km, s, Km/s)
#	h é o meio afastamento, m é o CMP
# 	nh número de receptores 
# 	dh incremento entre os receptores
# 	h0 orige dos receptores
# 	ns número de fontes
# 	ds incremento entre as fontes
# 	s0 origem das fontes
#---------------------------------------------------------
Flow('data','dome dip',
     '''
     kirmod cmp=y dip=${SOURCES[1]} 
     nh=161 dh=0.025 h0=0
     ns=401 ds=0.025 s0=0
     freq=10 dt=0.004 nt=1001
     vel=1.5 gradz=0.5 gradx=0.0 verb=y |
     put d2=0.0125
     ''')

Result('data',
       '''
       byte |
       transp plane=23 |
       grey3 flat=n frame1=500 frame3=80 frame2=200
       label1=Tempo unit1=s 
       label3=Meio-afastamento unit3=km 
       label2=PMC unit2=km
       title='Dados Modelo refletor gaussiano' point1=0.8 point2=0.8 
       ''')

Plot('data',
       '''
       byte |
       transp plane=23 memsize=1000 |
       grey3 flat=y frame1=500 frame3=80 frame2=200
       label1=Tempo unit1=s 
       label3=Meio-afastamento unit3=km 
       label2=PMC unit2=km
       title='Dados modelo refletor gaussiano' point1=0.8 point2=0.8 
       ''')

#---------------------{ Picking }-----------------------
# Extrai a superfície de tempo de trânsito de reflexão
# dos dados modelados
#-------------------------------------------------------

gplot = '''
transp |
graph3 frame1=2 frame3=80 frame2=200
wanttitle=n wantaxis=n min=4 max=0 plotfat=3
point1=0.8 point2=0.8 plotcol=3
'''

def gplot2(title,col,x0):
    return '''
    window min2=%g max2=%g |
    transp |
    graph3 frame1=2 frame3=80 frame2=%d
    title="%s" wantaxis=n min=3 max=1 plotfat=3 plotcol=%d
    ''' % (x0-dx,x0+dx,fx,title,4-col)

# TODO Explicar aqui a metodologia de picking!!!!
Flow('pick','data','envelope | max1 | window n1=1 | real | put label="Tempo" unit="s" ')

#---------------------------------{ VFSA }----------------------------------------
# 	Inversão dos parâmetros do CRS zero ofsset utilizando o algoritmo
#	Very Fast Simulated Aneeling (VFSA), descrito no livro de 
#	Stoffa & SEN, 1995 (Global optimization Methods in Geophysical Inversion) 
#	na página 106
#----------------------------------------------------------------------------------

v0=1.5 # Velocidade (Km/s)
app=1 # Índice da aproximação CRS a ser utilizada (veja o cabeçalho deste arquivo)
m0=5 # CMP central m0
verb=1 # Modo ativo, Manter assim! (Informa ao usuário sobre a aproximação utilizada)
temp0=10 # Temperatura inicial VFSA
c0=0.8 # Fator de amortecimento VFSA


for iter in range(2):

	par = 'param-%i' % (iter)

	otm = 'otm-%i' % (iter)

	erro = 'erro-%i' % (iter)

	err = 'err-%i' % (iter)

	out = 'out-%i'% (iter)

	# VFSA
	Flow([out,par],'data',
	'''
	vfsaSemb param=${TARGETS[1]} verb=%d app=%d m0=%g v0=%g temp0=%g c0=%g
	''' % (verb,app,m0,v0,temp0,c0))
	
	Plot(['pick','bar'],'pick','grey color=j allpos=y title="SRC-Modelada" label3="Km" unit3="s" label2="PMC" unit2="Km" label1="Afastamento" unit1="Km" scalebar=y verb=y clip=3.65 bias=0 minval=0 maxval=3.65')
	
	# Gerar superfície aproximada e superfície de erro relativo absoluto
	Flow(otm,['pick',par],
	'''
	crs param=${SOURCES[1]} verb=%d app=%d m0=%g v0=%g
	''' % (verb,app,m0,v0))
	Plot(otm,'grey color=j allpos=y clip=3.65 title="SRC-Otimizada" label3="Km" unit3="s" label2="PMC" unit2="Km" label1="Afastamento" unit1="Km" scalebar=y verb=y bias=0 minval=0 maxval=3.65')

	## Superficie de erro relativo absoluto (Aproximada - Modelada)
	Flow(erro,[otm,'pick'],'add scale=1,-1 ${SOURCES[1]} | math output="abs(input)" ')
	Plot(erro,'grey color=j allpos=y title="Erro Relativo Absoluto" label3="Km" unit3="s" label2="PMC" unit2="Km" label1="Afastamento" unit1="Km" scalebar=y verb=y clip=1 bias=0 minval=0 maxval=1')
	Result(err,[erro, otm, 'pick'],'SideBySideAniso',vppen='txscale=1.5')

End()
