/* Versão 1.0 - Gerar as coordenadas m de uma curva de tempo de trânsito CRE a partir de uma superfície de tempo de trânsito CRE e dos parâmetro do CRS de afastamento nulo (RN, RNIP e BETA)

Este programa depende dos parâmetros do CRS zero offset fornecidos no arquivo 'param.rsf' na variável 'par='. 
Este arquivo contém os parâmetros do CRS zero offset organizados na ordem RN, RNIP e BETA. Este arquivo pode ser gerado com o programa sfvfsa.
O arquivo de entrada padrão 'in.rsf' deste programa é uma superfície de tempo de trânsito CRS t(m,h), também gerada com sfvfsa.
Onde:
	x1	é o meio afastamento h
	x2	é a coordenada do cmp

A curva CRE necessita também da coordenada do CMP central m0 fornecida como parâmetro.

Exemplo de uso:

	<in.rsf sfplotcre m0=5 verb=1 par=param.rsf > out.rsf
	< out.rsf sfgraph > out.vpl
	sfpen out.vpl
	 
Site: http://www.dirackslounge.online
	 
Programador: Rodolfo A. C. Neves (Dirack) 17/03/2019
	 
Email: rodolfo_profissional@hotmail.com
	 
Referências sobre o método CRE:

	-CRUZ, J. et al. The common reflecting element (cre) method revisited. Geophysics, v. 65, p. 979–993, 2000.

*/
/*
  Copyright (C) 2018 grupo de programação aplicada à geofísica (GPGEOF)
  da Universidade Federal do Pará (UFPA); Belém, Pará, Brasil. 

  Esse programa é um software livre; você pode redistribuir e/ou modificar
  sobre os termos da licensa pública geral (LPG) publicada pela Free 
  Software Foundation; na versão 2 da licensa, ou (a seu critério) qualquer
  versão posterior.

  Este programa é distribuído na esperança que será útil, mas SEM NENHUMA
  GARANTIA; nem mesmo a garantia implícita de MERCANTILIDADE ou SERVENTIA
  A UM PROPÒSITO DETERMINADO. veja a LPG licensa pública geral para mais
  detalhes.

  Você deve ter recebido uma cópia da LPG licensa pública geral junto com
  este programa; se não, escreva para a Free Software Foundation, Inc., 
  51 Franklin Street, Quinquagésimo andar, Boston, MA  02110-1301, USA.
*/


#include <rsf.h> //incluindo a biblioteca do Madagascar
#include <stdio.h>
#include <math.h>

int main(int argc, char* argv[])
{

	int i1; // contador de laço
	int nh; // número de amostras no domínio meio afastamento
	int nm; // número de amostras no domínio do CMP
	float dh // Intervalo de amostragem no domínio meio afastamento
	float h0; // Origem no domínio meio afastamento
	float dm; // intervalo de amostragem no domínio CMP
	float x0; // origem no domínio CMP
	float m0; // CMP central (origem da curva CRE)
	float alpha, fatorAlpha; // Fator de assimetria do método CRE
	float RNIP, BETA; // Parâmetros do CRS zero offset
	float *cre_m; // Coordenadas m(h) da curva CRE
	float h; // Meio afastamento
	float *otm; // Parâmetros CRS zero offset (RN, RNIP, BETA)
	bool verb; // Selecionar modo ativo ou silencioso
	sf_axis ax, ay, az; // Eixos do arquivo de saída

	/* Inicializa arquivos rsf*/
	sf_file in, param, out; 

	/* Permite receber variáveis pela linhas de comando */
	sf_init(argc,argv); 

	in = sf_input("in"); // Superfície de tempo de trânsito CRS
	param = sf_input("param"); // Parâmetros CRS zero offset (RN, RNIP, BETA)
	out = sf_output("out"); // Coordenadas m,h da curva CRE m(h)

	if(! sf_getbool("verb",&verb)) verb=0;
	/* Modo= 1: modo ativo;	0: modo silencioso */

	if (verb) {

		sf_warning("Modo ativo ligado!!!");

	}

	/* Obtenha iformações sobre os eixos do arquivo de entrada 'in.rsf' */
	if (!sf_histint(in,"n1",&nh)) sf_error("Sem n1= no arquivo de entrada");
	if (!sf_histfloat(in,"d1",&dh)) sf_error("Sem d1= no arquivo de entrada");
	if (!sf_histfloat(in,"o1",&h0)) sf_error("Sem o1= no arquivo de entrada");

	if (!sf_histint(in,"n2",&nm)) sf_error("Sem n2= no arquivo de entrada");
	if (!sf_histfloat(in,"d2",&dm)) sf_error("Sem d2= no arquivo de entrada");
	if (!sf_histfloat(in,"o2",&x0)) sf_error("Sem o2= no arquivo de entrada");

	if (!sf_getfloat("m0",&m0)) m0=0.;
	/* Coordenada CMP central (Km)*/

	cre_m=sf_floatalloc(nh); // alocar memória vetor cre_m
	
	/* Ler os parâmetros do CRS do arquivo param */
	otm=sf_floatalloc(3);

	sf_floatread(otm,3,param);	
	     
	/* Parâmetros do CRS zero offset */
	RNIP=otm[1];
	BETA=otm[2];

	/* Calcular alpha - Fator de assimetria método CRE */
	alpha = sin(BETA)/RNIP;
	fatorAlpha = 1/(2*alpha);

	/**********[ Preparar o cabeçalho do arquivo de saída ]*************
	* Isso será feito com as seguintes funções:
	*
	* void sf_putint (sf_file file, const char* key, int par)
	* < Coloca um parâmetro int no cabeçalho do arquivo >
	*
	* void sf_putfloat (sf_file file, const char* key,float par)
	* < Coloca um parâmetro float no cabeçalho do arquivo >
	********************************************************************/
	sf_putint(out,"n1",nh); 
	sf_putfloat(out,"d1",dh); 
	sf_putfloat(out,"o1",h0);
	
	/* eixo = sf_maxa(n,o,d)*/
	ax = sf_maxa(nh, h0, dh);
	ay = sf_maxa(1,0,1);
	az = sf_maxa(1,0,1);
	sf_oaxa(out,ax,1);
	sf_oaxa(out,ay,2);
	sf_oaxa(out,az,3);

		/*****[ laço sobre as amostras ]****************
		* As coordenadas m da curva CRE são calculadas
		***********************************************/
		for (i1=0; i1 < nh; i1++) {

			h=h0+dh*i1;

			cre_m[i1]=m0+(fatorAlpha)*(1-sqrt(1+4*alpha*alpha*h*h));
		}
	
	sf_floatwrite(cre_m,nh,out); // Escreve m(h) no arquivo de saída

	exit(0);
}
