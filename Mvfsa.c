/* Versão 1.1 - Inversão dos parâmetros do CRS zero offset (RN, RNIP, BETA) utilizando o VFSA 

- Modificação: C0 e temp0 são Fornecidos pelo usuário.

A partir de uma superfície de tempo de trânsito de reflexão, extraída dos dados, se
obtém os parâmetros que melhor ajustam uma aproximação de tempo de trânsito CRS 
(Dentre as 11 opções a serem escolhidas) à superfície fornecida.
A otimização dos parâmetros é feita através do algoritmo Very Fast Simulated Aneeling (VFSA),
usando como critério de convergência o semblance entre a superfície de tempo de trânsito
fornecida com a superfície de tempo de trânsito aproximada dada uma aproximação de tempo
de trânsito CRS.
Para facilitar a convergência, forneça os valores iniciais de busca de RN, RNIP e BETA. 
A variável app identifica a aproximação de tempo de trânsito CRS a ser utilizada.
O arquivo de entrada padrão 'in.rsf' é uma superfície de tempo de trânsito crs t(m,h). Onde:
	x1	é o meio afastamento h
	x2	é a coordenada do cmp

A aproximação é dada em função de um CMP central m0 e da velocidade v0 próxima a superfície.

Exemplo de uso:

	<in.rsf sfvfsa param=param.rsf verb=1 app=1 m0=5 v0=1.5 rn=2000 rnip=1000 beta=0 > out.rsf
	< out.rsf sfgrey > out.vpl
	sfpen out.vpl

Referências sobre as aproximações de tempo de trânsito CRS utilizadas:

	-FOMEL, S.; KAZINNIK, R. Nonhyperbolic common reflection surface. Geophysical Prospecting, v. 61, p. 21–27, 2013.

	-HöCHT, G. Traveltime approximations for 2D and 3D media and kinematic wavefield attributes. Tese (Doutorado) — Faculdade de Física Karlsruhe (TH) genehmigte, 2002.

	-NEVES, R. Aproximações não hiperbólicas do tempo de trânsito utilizando aproximantes de Padé. Dissertação (Mestrado) — Universidade Federal do Pará - UFPa, Belém - PA, 2017.

	-JAGER, R. et al. Common-reflection-surface stack: image and attributes. Geophysics,v. 66, p. 97–109, 2001.

Referências sobre o algoritmo VFSA:

	-INGBER, L. Very fast simulated re-annealing. Math1. Comput. Modelling, v. 12, p.967–973, 1989.

Programador: Rodolfo A. C. Neves 12/11/2018

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

#include <rsf.h> //incluindo a biblioteca do Madagascar

#include "f_vfsa.h" //biblioteca de funções do SFVFSA

int main(int argc, char* argv[])
{
   
	float x0; // origem do eixo dos CMP's
	float m0; // CMP central 
	float h0; // origem do eixo do meio-afastamento
	float dh; // intervalo de amostragem no domínio do meio afastamento
	float dm; // intervalo de amostragem entre os CMP's
	float *otm; //Parâmetros otimizados
	int nh; // número de amostras no eixo do meio afastamento
	int nm; // número de amostras no eixo do CMP
	bool verb; // Chave modo silencioso 0 e modo ativo 1
	int app; // Id da aproximação de tempo de trânsito escolhida
	float **t; // Superfície de tempo de trânsito modelada (in.rsf)
	int im0; // índice da amostra do CMP central m0
	float t0; // tempo de trânsito do raio normal
	float temp; // temperatura na iteração 'q' do VFSA
	float temp0; // temperatura inicial do VFSA
	float c0; // Parâmetro do VFSA
	float y; // perturbação dos parâmetros (incremento)
	int q; // contador de iteração
	int p; // contador de parâmetro
	float c[3]; // vetor de parâmetros do crs invertidos
	float mmin[3], mmax[3]; // vetor q define janela de busca
	float cnew[3]; //vetor temporario de parâmetros
	int itmax=25000; // número máximo de interações no VFSA
	float u; // numero aleatório entre 0 e 1
	float deltaE; // Semblamce iteração atual - iteração anterior (critério de convergência)
	float Em0=0; // semblance da iteração anterior
	float PM; // Parâmetro do critério de metrópolis
	float R_N, R_NIP, BETA; // parâmetros do crs zero offset
	float semb;// valor de semblance em uma iteração qualquer
	float otsemb, otrn, otrnip, otbeta; // parâmetros ótimos (resultado da inversão)
	float t_R_NIP; // Variável temporária
	float v0; // velocidade próxima a superfície
	char *app_s; // Armazena o nome da aproximação

	/* Parâmetros da linha de comando, valores iniciais de busca */
	float RN_in;
	float RNIP_in;
	float BETA_in;
	float semb_in; 

/**************************[ Configuração ]*******************************************************
*
*	Abrir arquivos para armazenar os resultados das etapas das iterações. Inicializar variáveis
*	e informar o usuário caso haja algum erro.
*
***************************************************************************************************/
	
	/* Inicializa arquivos rsf*/
	sf_file in, out,param;

	/* Inicializa eixos */
	sf_axis ax,ay,az;

	/* Permite receber variáveis pela linhas de comando */
	sf_init(argc,argv); 

	/* Arquivos .rsf de I/O */
	in = sf_input("in"); // Superficie de tempo de trânsito CRS modelada 
	out = sf_output("out"); // Superfície de tempo de trânsito CRS aproximada
	param = sf_output("param"); // Parâmetros do CRS otimizados (RN, RNIP, BETA)

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
	if (!sf_histint(in,"n1",&nh)) sf_error("Sem n1= no arquivo de entrada");
	if (!sf_histfloat(in,"d1",&dh)) sf_error("Sem d1= no arquivo de entrada");
	if (!sf_histfloat(in,"o1",&h0)) sf_error("Sem o1= no arquivo de entrada");

	if (!sf_histint(in,"n2",&nm)) sf_error("Sem n2= no arquivo de entrada");
	if (!sf_histfloat(in,"d2",&dm)) sf_error("Sem d2= no arquivo de entrada");
	if (!sf_histfloat(in,"o2",&x0)) sf_error("Sem o2= no arquivo de entrada");

	if (!sf_getfloat("m0",&m0)) m0=0.;
	/* CMP central (Km) */

	if (!sf_getfloat("v0",&v0)) v0=1.5;
	/* Velocidade próxima a superfície (Km/s) */
              
	if (!sf_getfloat("rn",&RN_in)) RN_in=1.;
	/* raio de curvatura N inicial */

	if (!sf_getfloat("rnip",&RNIP_in)) RNIP_in=1.;
	/* raio de curvatura NIP inicial */

	if (!sf_getfloat("beta",&BETA_in)) BETA_in=0.;
	/* ângulo de emergência do raio normal inicial */

	if (!sf_getfloat("semb",&semb_in)) semb_in=0.;
	/* Semblamce inicial */

	if (!sf_getfloat("c0",&c0)) c0=0.5;
	/* Fator de amortecimento */

	if (!sf_getfloat("temp0",&temp0)) temp0=10.;
	/* Temperatura inicial */
    
	/* Ler a superfície de tempo de trânsito modelada */
	t=sf_floatalloc2(nh,nm);

	sf_floatread(t[0],nh*nm,in);

	/* Indice da amostra do CMP central m0*/
	im0=(m0/dm);	

	/* t0 tempo de trânsito do raio normal */
	t0=t[im0][0];
	
	/* Inicializando parâmetros a inverter */
	c[0]=RN_in;  // RN inicial
	c[1]=RNIP_in;  // RNIP inicial
	c[2]=BETA_in;  // BETA inicial
	
	/* Definindo o intervalo de busca dos Parâmetros */
	mmin[0]=0.5;
	mmax[0]=5;
	mmin[1]=0.5;
	mmax[1]=5;
	mmin[2]=-3.14159;
	mmax[2]=3.14159;
	
///*****************************[ Inversão dos parâmetros do CRS utilizando VFSA ]**********************
//*
//*	Aqui inicia o algoritmo de inversão Very Fast Simulated Aneeling (VFSA), descrito no livro
//*	Stoffa & SEN, 1995 (Global optimization Methods in Geophysical Inversion) pg. 106
//*
//******************************************************************************************************/
		
	/*loop sobre as iterações*/
	for (q=0; q <itmax; q++){
		
		/* Temperatura da iteração */
		temp=temp0*expf(-c0*pow(q,0.25));
		
		
		/*loop sobre os parâmetros a inverter (3 parâmetros)*/
		for (p=0; p < 3; p++){
			
			/* gerar semente de número aleatório */
			srand(time(NULL)*q*p);
			
			/* número aleatório entre 0 e 1 */
			u=(float)(rand()%1000)/1000;
			
			/* perturbação no parâmetro */
			y = sinal(u - 0.5) * temp * (pow( (1+temp),fabs(2*u-1) )-1);			
			
			cnew[p] = c[p] + y * (mmax[p]-mmin[p]);
			
				/*condição de janelamento dos parâmetros*/
				if (cnew[p] >= mmax[p]) {
					
					srand(time(NULL)*q*p*u);
					u=(float)(rand()%1000)/1000;
					cnew[p] = (mmax[p]-mmin[p]) * u + mmin[p];
				} else if (cnew[p] <= mmin[p]) {
					
					srand(time(NULL)*q*p*u);
					u=(float)(rand()%1000)/1000;
					cnew[p] = (mmax[p]-mmin[p]) * u + mmin[p];
					
				}
									
		}
		
		
		/* os parâmetros são atualizados 
		 após serem perturbados e janelados */
		R_N=cnew[0];	
		R_NIP=cnew[1];
		BETA=cnew[2];
		
		/* Restrição: RN deve ser maior que RNIP*/ 
		if (R_N < R_NIP) {
			t_R_NIP=R_N;
			R_N=R_NIP;
			R_NIP=t_R_NIP;
		}
		
	
		/* inicializando o valor do semblance */
		semb=semb_in;
		
		/* Escolhendo a aproximação de tempo de trânsito CRS */
		
		semb=1.;
		
		switch (app)
		{
		   case 1: //aproximação fomel (CRS NÃO hiperbólico)
			app_s="Fomel";
			semb=fomel(t0, m0,  h0, x0, v0,  R_N,  R_NIP,  BETA,  nh,  dh,  nm,  dm, t);
			sf_warning("Aproximação (%s) iteração=%i/25000 Semb=%f",app_s,q,semb);
		   break;

		   case 2: //aproximação jager (CRS hiperbólico)
			app_s="Jager";
			semb=jager(t0, m0,  h0, x0, v0,  R_N,  R_NIP,  BETA,  nh,  dh,  nm,  dm, t);
			sf_warning("Aproximação (%s) iteração=%i/25000 Semb=%f",app_s,q,semb);
		   break;
	   
		   case 3: //aproximação germam-t (CRS quarta ordem - parabólico)
			app_s="Germam t";
			semb=germam_t(t0, m0,  h0, x0, v0,  R_N,  R_NIP,  BETA,  nh,  dh,  nm,  dm, t);
			sf_warning("Aproximação (%s) iteração=%i/25000 Semb=%f",app_s,q,semb);
		   break;
		   
		   case 4: //aproximação germam-t2 (CRS quarta ordem - quadrático)
			app_s="Germam t2";
			semb=germam_t2(t0, m0,  h0, x0, v0,  R_N,  R_NIP,  BETA,  nh,  dh,  nm,  dm, t);
			sf_warning("Aproximação (%s) iteração=%i/25000 Semb=%f",app_s,q,semb);
		   break;
		   
		   case 5: //aproximação germam-tshift (CRS quarta ordem - hipérbole deslocada)
			app_s="Germam tshift";
			semb=germam_tshift(t0, m0,  h0, x0, v0,  R_N,  R_NIP,  BETA,  nh,  dh,  nm,  dm, t);
			sf_warning("Aproximação (%s) iteração=%i/25000 Semb=%f",app_s,q,semb);
		   break;
		   
		   case 6: //aproximação Padé parabólico (CRS Padé parabólico expansão em h)
			app_s="Padé th";
			semb=pade_th(t0, m0,  h0, x0, v0,  R_N,  R_NIP,  BETA,  nh,  dh,  nm,  dm, t);
			sf_warning("Aproximação (%s) iteração=%i/25000 Semb=%f",app_s,q,semb);
		   break;
		   
		   case 7: //aproximação Padé parabólico (CRS Padé parabólico expansão em m)
			app_s="Padé tm";
			semb=pade_tm(t0, m0,  h0, x0, v0,  R_N,  R_NIP,  BETA,  nh,  dh,  nm,  dm, t);
			sf_warning("Aproximação (%s) iteração=%i/25000 Semb=%f",app_s,q,semb);
		   break;
		   
		   case 8: //aproximação Padé hiperbólico (CRS Padé hiperbólico expansão em h)
			app_s="Padé t2h";
			semb=pade_t2h(t0, m0,  h0, x0, v0,  R_N,  R_NIP,  BETA,  nh,  dh,  nm,  dm, t);
			sf_warning("Aproximação (%s) iteração=%i/25000 Semb=%f",app_s,q,semb);
		   break;
		   
		   case 9: //aproximação Padé hiperbólico (CRS Padé hiperbólico expansão em m)
			app_s="Padé t2m";
			semb=pade_t2m(t0, m0,  h0, x0, v0,  R_N,  R_NIP,  BETA,  nh,  dh,  nm,  dm, t);
			sf_warning("Aproximação (%s) iteração=%i/25000 Semb=%f",app_s,q,semb);
		   break;
		   
		   case 10: //aproximação Padé Deslocado (CRS Padé hipérbole deslocada expansão em h)
			app_s="Padé tsh";			 
			semb=pade_tsh(t0, m0,  h0, x0, v0,  R_N,  R_NIP,  BETA,  nh,  dh,  nm,  dm, t);
			sf_warning("Aproximação (%s) iteração=%i/25000 Semb=%f",app_s,q,semb);
		   break;
		   
		   case 11: //aproximação Padé Deslocado (CRS Padé hipérbole deslocada expansão em m)
			app_s="Padé tsm";
			semb=pade_tsm(t0, m0,  h0, x0, v0,  R_N,  R_NIP,  BETA,  nh,  dh,  nm,  dm, t);
			sf_warning("Aproximação (%s) iteração=%i/25000 Semb=%f",app_s,q,semb);
		   break;

		   default:
			 sf_error("Opção app=%i Não disponível", app);
		}

						
		/* condição de convergência dos parâmetros no algoritmo VFSA */		
		if(fabs(semb) > fabs(Em0) ){
			otsemb=semb;
			otrn=R_N;
			otrnip=R_NIP;
			otbeta=BETA;
			
		}
		
		/* condições de atualização dos parâmetros do VFSA */
		deltaE = -semb - Em0;
		
		/* Critério de metrópolis */
		PM = expf(-deltaE/temp);
		
		if (deltaE<=0){
			c[0]=cnew[0];
			c[1]=cnew[1];
			c[2]=cnew[2];
			Em0=-semb;
			
		} else {
			
			u=rand();
			
			if (PM > u){
				
				c[0]=cnew[0];
				c[1]=cnew[1];
				c[2]=cnew[2];
				Em0=-semb;
				
			}
			
		}	
		
	}

	/* Salvar os parâmetros otimizados no arquivo 'param' */
	otm=sf_floatalloc(3);

	otm[0] = otrn;
	otm[1] = otrnip;
	otm[2] = otbeta;

	/* Mostre os parâmetros otimizados ao usuário antes de salvar */
	sf_warning("Parâmetros otimizados:\n RN=%f, RNIP=%f, BETA=%f, SEMB=%f",otrn,otrnip,otbeta,otsemb);

	/* eixo = sf_maxa(n,o,d)*/
	ax = sf_maxa(3, 0, 1);
	ay = sf_maxa(1, 0, 1);
	az = sf_maxa(1, 0, 1);

	/* sf_oaxa(arquivo, eixo, índice do eixo) */
	sf_oaxa(param,ax,1);
	sf_oaxa(param,ay,2);
	sf_oaxa(param,az,3);
	sf_floatwrite(otm,3,param);

	/* Escrever a superfície otimizada no arquivo 'out'*/
	sf_floatwrite(t[0],nm*nh,out);
	
	exit(0);
}
