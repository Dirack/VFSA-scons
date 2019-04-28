/* Versão 2.0 - Gerar cubo de coerência dos parâmetros do CRS zero offset (RN, RNIP, BETA)

Montar o cubo de coerência (semblance) para os parâmetros do CRS zero offset. 
Utiliza o cubo de dados A(m,h,t).

Modificação: Agora o programa utiliza a função Semblance da biblioteca 'f_vfsa.h' para calcular o 
semblance em vários pontos do cubo de coerência.

Exemplo de uso:

		TODO

Programador: Rodolfo A. C. Neves 15/04/2019 (1.0 - original)
	     Rodolfo A. C. Neves 27/04/2019 (2.0)

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
#include "f_vfsaSemb.h" // Calcular Semblance

int main(int argc, char* argv[]){

	float ***t;
	float ***cubo_;
	int i,j,k;
	float RN, RNIP, BETA;
	float m0;
	float v0;
	float t0;
	int im0;
	int nt, nh, nm;
	float dm, dt, dh, x0, h0, ot;
	bool verb;
	int app;
	char* app_s;
	int nBETA=50, nRN=50, nRNIP=50;
	float dBETA=0.04, dRN=0.02, dRNIP=0.02;
	float oBETA=-1, oRN=1.2, oRNIP=1.2;

	/* ar é o eixo de RN e RNIP no cubo de coerência */
	/* ab é o eixo de BETA no cubo de coerência */
	sf_axis arn, arnip, ab;

	/* Inicializa arquivos rsf*/
	sf_file in, out;

	/* Permite receber variáveis pela linhas de comando */
	sf_init(argc,argv);


	/* Arquivos .rsf de I/O */
	in = sf_input("in"); /* cubo de dados A(m,h,t) */
	out = sf_output("out"); /* cubo de coerência */


	if (!sf_getint("app",&app)) app=1;
	/* Aproximação de tempo de trânsito CRS:  
 		-fomel [ app=1 ]
 		-Jager [ app=2 ]
		-Germam-t [ app=3 ]
		-Germam-t2 [ app=4 ]
		-Germam-tshift [ app=5 ]
		-Padé-t-h [ app=6 ]
		-Padé-t-m [ app=7 ]
		-Padé-t2-h [ app=8 ]
		-Padé-t2-m [ app=9 ]
		-Padé-tshift-h [ app=10 ]
		-Padé-tshift-m [ app=11 ]
	*/

	if(! sf_getbool("verb",&verb)) verb=0;
	/* Modo= 1: modo ativo;	0: modo silencioso */

	if (verb) {

		sf_warning("Modo ativo ligado!!!");
		
		/* Avisar o usuário qual aproximação de tempo de 
		trânsito CRS está sendo utilizada */
		f_vfsa_aviso(app);
	}


	/* Obtenha a informação sobre os eixos do arquivo de entrada 'in.rsf' */
	if (!sf_histint(in,"n1",&nt)) sf_error("Sem n1= no arquivo de entrada");
	if (!sf_histfloat(in,"d1",&dt)) sf_error("Sem d1= no arquivo de entrada");
	if (!sf_histfloat(in,"o1",&ot)) sf_error("Sem o1= no arquivo de entrada");

	if (!sf_histint(in,"n2",&nh)) sf_error("Sem n2= no arquivo de entrada");
	if (!sf_histfloat(in,"d2",&dh)) sf_error("Sem d2= no arquivo de entrada");
	if (!sf_histfloat(in,"o2",&h0)) sf_error("Sem o2= no arquivo de entrada");

	if (!sf_histint(in,"n3",&nm)) sf_error("Sem n3= no arquivo de entrada");
	if (!sf_histfloat(in,"d3",&dm)) sf_error("Sem d3= no arquivo de entrada");
	if (!sf_histfloat(in,"o3",&x0)) sf_error("Sem o3= no arquivo de entrada");

	if (!sf_getfloat("m0",&m0)) m0=0.;
	/* CMP central (Km) */

	if (!sf_getfloat("v0",&v0)) v0=1.5;
	/* Velocidade próxima a superfície (Km/s) */
              
	/* Ler o cubo de dados */
	t=sf_floatalloc3(nt,nh,nm);
	sf_floatread(t[0][0],nh*nm*nt,in);

	/* Indice da amostra do CMP central m0*/
	im0=(m0/dm);	

	//sf_warning("im0=%i;dt=%f;nt=%i;m0=%f;dm=%f",im0,dt,nt,m0,dm);

	/* t0 tempo de trânsito do raio normal */
	t0=encontrarT0(t, im0, dt, nt);//im0*dt;
	//sf_error("t0=%f",t0);

	/* Montar o grid do cubo de coerência */
	// RN e RNIP irão de 0 a 5Km
	// BETA varia de -3.14 a 3.14
	
	/* eixo = sf_maxa(n,o,d)*/
	arn = sf_maxa(nRN,oRN,dRN);
	arnip = sf_maxa(nRNIP,oRNIP, dRNIP);
	ab = sf_maxa(nBETA, oBETA, dBETA);

	sf_setlabel(arn,"RN");
	sf_setlabel(arnip,"RNIP");
	sf_setlabel(ab,"BETA");

	/* sf_oaxa(arquivo, eixo, índice do eixo) */
	sf_oaxa(out,arn,1);
	sf_oaxa(out,arnip,2);
	sf_oaxa(out,ab,3);

	cubo_ = sf_floatalloc3(nRN,nRNIP,nBETA);	


	/* Inserir dados no cubo de coerência */
	for(i=0;i<nBETA;i++){

		BETA = i*dBETA + oBETA;

		for(j=0;j<nRNIP;j++){

 			RNIP = j*dRNIP + oRNIP;

			for(k=0;k<nRN;k++){

				RN = k*dRN + oRN;
	
				/* Escolhendo a aproximação de tempo de trânsito CRS */
				switch (app)
				{
				   case 1: //aproximação fomel (CRS NÃO hiperbólico)
					app_s="Fomel";
					cubo_[i][j][k]=1.2E-11;
					cubo_[i][j][k]=fomel(dt,t0, m0,  h0, x0, v0,  RN,  RNIP,  BETA,  nh,  dh,  nm,  dm, t);
					sf_warning("(%i,%i,%i)RN=%f;RNIP=%f;BETA=%f;SEMB=%f",i,j,k,RN,RNIP,BETA,cubo_[i][j][k]);
				   break;

				   case 2: //aproximação jager (CRS hiperbólico)
					app_s="Jager";
					cubo_[i][j][k]=jager(t0, m0,  h0, x0, v0,  RN,  RNIP,  BETA,  nh,  dh,  nm,  dm, t);
				   break;
			   
				   case 3: //aproximação germam-t (CRS quarta ordem - parabólico)
					app_s="Germam t";
					cubo_[i][j][k]=germam_t(t0, m0,  h0, x0, v0,  RN,  RNIP,  BETA,  nh,  dh,  nm,  dm, t);
				   break;
				   
				   case 4: //aproximação germam-t2 (CRS quarta ordem - quadrático)
					app_s="Germam t2";
					cubo_[i][j][k]=germam_t2(t0, m0,  h0, x0, v0,  RN,  RNIP,  BETA,  nh,  dh,  nm,  dm, t);
				   break;
				   
				   case 5: //aproximação germam-tshift (CRS quarta ordem - hipérbole deslocada)
					app_s="Germam tshift";
					cubo_[i][j][k]=germam_tshift(t0, m0,  h0, x0, v0,  RN,  RNIP,  BETA,  nh,  dh,  nm,  dm, t);
				   break;
				   
				   case 6: //aproximação Padé parabólico (CRS Padé parabólico expansão em h)
					app_s="Padé th";
					cubo_[i][j][k]=pade_th(t0, m0,  h0, x0, v0,  RN,  RNIP,  BETA,  nh,  dh,  nm,  dm, t);
				   break;
				   
				   case 7: //aproximação Padé parabólico (CRS Padé parabólico expansão em m)
					app_s="Padé tm";
					cubo_[i][j][k]=pade_tm(t0, m0,  h0, x0, v0,  RN,  RNIP,  BETA,  nh,  dh,  nm,  dm, t);
				   break;
				   
				   case 8: //aproximação Padé hiperbólico (CRS Padé hiperbólico expansão em h)
					app_s="Padé t2h";
					cubo_[i][j][k]=pade_t2h(t0, m0,  h0, x0, v0,  RN,  RNIP,  BETA,  nh,  dh,  nm,  dm, t);
				   break;
				   
				   case 9: //aproximação Padé hiperbólico (CRS Padé hiperbólico expansão em m)
					app_s="Padé t2m";
					cubo_[i][j][k]=pade_t2m(t0, m0,  h0, x0, v0,  RN,  RNIP,  BETA,  nh,  dh,  nm,  dm, t);
				   break;
				   
				   case 10: //aproximação Padé Deslocado (CRS Padé hipérbole deslocada expansão em h)
					app_s="Padé tsh";			 
					cubo_[i][j][k]=pade_tsh(t0, m0,  h0, x0, v0,  RN,  RNIP,  BETA,  nh,  dh,  nm,  dm, t);
				   break;
				   
				   case 11: //aproximação Padé Deslocado (CRS Padé hipérbole deslocada expansão em m)
					app_s="Padé tsm";
					cubo_[i][j][k]=pade_tsm(t0, m0,  h0, x0, v0,  RN,  RNIP,  BETA,  nh,  dh,  nm,  dm, t);
				   break;

				   default:
					 sf_error("Opção app=%i Não disponível", app);
				}

			}

		}

	}

	
	/* Escrever a superfície otimizada no arquivo 'out'*/
	sf_floatwrite(cubo_[0][0],nBETA*nRN*nRNIP,out);
	
	exit(0);

}
