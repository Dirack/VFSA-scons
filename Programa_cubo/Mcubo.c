/* Versão 1.0 - Gerar cubo de coerência dos parâmetros do CRS zero offset (RN, RNIP, BETA)

Montar o cubo de coerência (semblance) para os parâmetros do CRS zero offset. Utiliza os arquivos 
gerados pelo programa 'sfvfsaSemb' com os valores de semblance, RN, RNIP e BETA.

Exemplo de uso:

		TODO

Programador: Rodolfo A. C. Neves 15/04/2019

Email:  rodolfo_profissional@hotmail.com  

Acesse conteúdo exclusivo e tire dúvidas no nosso site:
	http://www.dirackslounge.online

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

#include <rsf.h> //Biblioteca do MADAGASCAR

int main(int argc, char* argv[]){

	float *rn_;
	float *rnip_;
	float *beta_;
	float *semb_;
	float ***cubo_;
	int n,i,j,k;
	int RN, RNIP, BETA;

	/* ar é o eixo de RN e RNIP no cubo de coerência */
	/* ab é o eixo de BETA no cubo de coerência */
	sf_axis arn, arnip, ab;

	/* Inicializa arquivos rsf*/
	sf_file in, rn_in, rnip_in, beta_in, out;

	/* Permite receber variáveis pela linhas de comando */
	sf_init(argc,argv);

	/* Verificar se os arquivos de entrada foram fornecidos */
	if(!sf_getstring("rn_in")) sf_error("Arquivo 'rn_in=' não foi fornecido pelo usuário.");
	if(!sf_getstring("rnip_in")) sf_error("Arquivo 'rnip_in=' não foi fornecido pelo usuário.");
	if(!sf_getstring("beta_in")) sf_error("Arquivo 'beta_in=' não foi fornecido pelo usuário.");

	/* Arquivos .rsf de I/O */
	in = sf_input("in");
	rn_in = sf_input("rn_in");
	rnip_in = sf_input("rnip_in");
	beta_in = sf_input("beta_in");
	out = sf_output("out");

	/* Receber tamanho dos vetores */
	if (!sf_histint(in,"n1",&n)) sf_error("Sem n1= no arquivo de entrada");

	/* Alocar memória e ler os vetores */
	semb_ = sf_floatalloc(n);
	rn_ = sf_floatalloc(n);
	rnip_ = sf_floatalloc(n); 
	beta_ = sf_floatalloc(n);	

	sf_floatread(semb_ ,n,in);
	sf_floatread(rn_,n,rn_in);
	sf_floatread(rnip_,n,rnip_in); 
	sf_floatread(beta_,n,beta_in);

	/* Montar o grid do cubo de coerência */
	// RN e RNIP irão de 0 a 5Km
	// BETA varia de -3.14 a 3.14
	
	/* eixo = sf_maxa(n,o,d)*/
	arn = sf_maxa(250,0, 0.02);
	arnip = sf_maxa(250,0, 0.02);
	ab = sf_maxa(314, -3.14, 0.02);

	sf_setlabel(arn,"RN");
	sf_setlabel(arnip,"RNIP");
	sf_setlabel(ab,"BETA");

	/* sf_oaxa(arquivo, eixo, índice do eixo) */
	sf_oaxa(out,arn,1);
	sf_oaxa(out,arnip,2);
	sf_oaxa(out,ab,3);

	cubo_ = sf_floatalloc3(250,250,314);	

	for(i=0;i<314;i++){

		BETA = i*0.02 -3.14;

		for(j=0;j<250;j++){

 			RNIP = j*0.02;

			for(k=0;k<250;k++){

				RN = k*0.02;
				cubo_[i][j][k]=1.15E-12;
				//sf_warning("%f;%f;%f",RN,RNIP,BETA);

			}

		}

	}

	/* Inserir dados no cubo de coerência */
	for(i=0;i<n;i++){

		RN = rn_[i]/0.02;
		RNIP = rnip_[i]/0.02;
		BETA = (beta_[i]+3.14)/0.02;

		//sf_warning("%f\n",semb_[i]);

		cubo_[BETA][RNIP][RN] = semb_[i];
	}

	/* Escrever a superfície otimizada no arquivo 'out'*/
	sf_floatwrite(cubo_[0][0],250*250*314,out);
	
	exit(0);

}
