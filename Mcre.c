/* Versão 1.0 - Gerar a superfície de tempo de trânsito CRE zero offset a partir de uma superfície de tempo de trânsito CRS e dos parâmetros (RN,RNIP, BETA)

Para montar a superfície CRE a partir de uma superfície CRS se utiliza a condição RN=RNIP (condição CRE) nas
aproximações de tempo de trânsito CRS. Por isto este programa é igual ao programa de plotagem da superfície CRS
com a única modificação de igualar os parâmetros RN e RNIP.

A variável app define uma aproximação de tempo de trânsito CRS a ser utilizada.
esta aproximação depende dos parâmetros do CRS zero offset fornecidos no arquivo .rsf
na variável par=, na ordem (RNIP e BETA), que pode ser gerado com sfvfsa.
O arquivo de entrada padrão .rsf é uma superfície de tempo de trânsito crs t(m,h). Onde:
	x1	é o meio afastamento h
	x2	é a coordenada do cmp

A aproximação é dada em função de um CMP central m0 e da velocidade v0 próxima a superfície

Exemplo de uso:

	<in.rsf sfcre app=1 m0=5 v0=1.5 verb=1 > out.rsf
	< out.rsf sfgrey > out.vpl
	sfpen out.vpl

Referências sobre o método CRE:

	-CRUZ, J. et al. The common reflecting element (cre) method revisited. Geophysics, v. 65, p. 979–993, 2000.

Referências sobre as aproximações de tempo de trânsito CRS utilizadas:

	-FOMEL, S.; KAZINNIK, R. Nonhyperbolic common reflection surface. Geophysical Prospecting, v. 61, p. 21–27, 2013.

	-HöCHT, G. Traveltime approximations for 2D and 3D media and kinematic wavefield attributes. Tese (Doutorado) — Faculdade de Física Karlsruhe (TH) genehmigte, 2002.

	-NEVES, R. Aproximações não hiperbólicas do tempo de trânsito utilizando aproximantes de Padé. Dissertação (Mestrado) — Universidade Federal do Pará - UFPa, Belém - PA, 2017.

	-JAGER, R. et al. Common-reflection-surface stack: image and attributes. Geophysics,v. 66, p. 97–109, 2001.

Programador: Rodolfo A. C. Neves (Dirack) 11/03/2019

Email:  rodolfo_profissional@hotmail.com  

Acesse conteúdo exclusivo e tire dúvidas no nosso site:
	http://www.dirackslounge.online

*/
/*
  Copyright (C) 2019 grupo de programação aplicada à geofísica (GPGEOF)
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

#include "f_cre.h" // Biblioteca de funções deste programa

int main(int argc, char* argv[])
{	
	float x0; // Origem da coordenada do CMP dos dados
	int im0; // índice do CMP
	float m0; // CMP central da aproximação
	float dm; // intervalo de amostragem dos CMP's
	int nm; // Número de CMP's
	float h0; // Origem da coordenada do meio afastamento
	float dh; // Intervalo de amostragem entre os meios-afastamentos
	int nh; // Número de meio-afastamentos
	float *otm; // Parâmetros CRS zero offset (RN, RNIP, BETA)
	bool verb; // Chave modo silencioso 0 e modo ativo 1
	int app; // Id da aproximação de tempo de trânsito escolhida
	float **t; // Superfície de tempo de trânsito CRS modelada
	float t0; // tempo de trânsito do raio normal
	float RN, RNIP, BETA; // parâmetros do crs zero offset
	float v0; // velocidade próxima a superfície
	char *app_s; // String armazena o nome da aproximação CRS escolhida

/**************************[ Configuração ]*******************************************************
*
*	Abrir arquivos .rsf de entrada e saída e inicializar os parâmetros obitidos a partir da 
*	linhas de comandos.
*
***************************************************************************************************/
	
	/* Inicializa arquivos rsf*/
	sf_file in, out, param; // Arquivos

	/* Permite receber variáveis pela linhas de comando */
	sf_init(argc,argv); 

	in = sf_input("in"); // Superfície de tempo de trânsito CRS
	param = sf_input("param"); // Parâmetros CRS zero offset (RN, RNIP, BETA)
	out = sf_output("out"); // Superfície CRE aproximada

	/***********************[ Configuração ]****************************/

	if (!sf_getint("app",&app)) app=1;
	/* Aproximação de tempo de trânsito CRS/CRE:  
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

		/*Informe ao usuário qual 
		aproximação CRS/CRE está sendo utilizada */
		f_vfsa_aviso(app); 
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

	if (!sf_getfloat("v0",&v0)) v0=1.5;
	/* Velocidade próxima a superfície (Km/s) */

	/* Ler os parâmetros do CRS de par */
	otm=sf_floatalloc(3);

	sf_floatread(otm,3,param);	
	     
	/* Para uma aproximação CRE RN=RNI (condição CRE) */
	RNIP=otm[1];
	RN=RNIP;
	BETA=otm[2];
    
	/* Ler a superfície CRS modelada */
	t=sf_floatalloc2(nh,nm);

	sf_floatread(t[0],nh*nm,in);

	/* Índice da amostra do CMP central m0 */
	im0=(m0/dm);

	/* t0  tempo de trânsito do raio normal */
	t0=t[im0][0];
	
/************************[ Montar a superfície de tempo de trânsito CRS aproximada ]**********************
*
*	Cada função abaixo corresponde a uma aproximação de tempo de trânsito CRS zero offset, escolhida
*	com a variável app (cada índice app corresponde a uma aproximação a ser utilizada).
*
**********************************************************************************************************/
		
	
		switch (app)
		{
		   case 1: //aproximação fomel (CRS NÃO hiperbólico)
			app_s="Fomel";
			fomel(t0, m0,  h0,  x0,  v0,  RN,  RNIP,  BETA, nh,  dh, nm,  dm,  t);
			sf_warning("Aproximação (%s) RN=%f RNIP=%f BETA=%f",app_s,RN,RNIP,BETA);
		   break;

		   case 2: //aproximação jager (CRS hiperbólico)
			app_s="Jager";			
			jager(t0, m0,  h0,  x0,  v0,  RN,  RNIP,  BETA, nh,  dh, nm,  dm,  t);
			sf_warning("Aproximação (%s) RN=%f RNIP=%f BETA=%f",app_s,RN,RNIP,BETA);
		   break;
		   
		   case 3: //aproximação germam-t (CRS quarta ordem - parabólico)
			app_s="Germam t";	
			germam_t(t0, m0,  h0,  x0,  v0,  RN,  RNIP,  BETA, nh,  dh, nm,  dm,  t);
			sf_warning("Aproximação (%s) RN=%f RNIP=%f BETA=%f",app_s,RN,RNIP,BETA);
		   break;
		   
		   case 4: //aproximação germam-t2 (CRS quarta ordem - quadrático)
			app_s="Germam t2";
			germam_t2(t0, m0,  h0,  x0,  v0,  RN,  RNIP,  BETA, nh,  dh, nm,  dm,  t);
			sf_warning("Aproximação (%s) RN=%f RNIP=%f BETA=%f",app_s,RN,RNIP,BETA);
		   break;
		   
		   case 5: //aproximação germam-tshift (CRS quarta ordem - hipérbole deslocada)
			app_s="Germam tshift";
			germam_tshift(t0, m0,  h0,  x0,  v0,  RN,  RNIP,  BETA, nh,  dh, nm,  dm,  t);
			sf_warning("Aproximação (%s) RN=%f RNIP=%f BETA=%f",app_s,RN,RNIP,BETA);
		   break;
		   
		   case 6: //aproximação Padé parabólico (CRS Padé parabólico expansão em h)
			app_s="Padé th";
			pade_th(t0, m0,  h0,  x0,  v0,  RN,  RNIP,  BETA, nh,  dh, nm,  dm,  t);
			sf_warning("Aproximação (%s) RN=%f RNIP=%f BETA=%f",app_s,RN,RNIP,BETA);
		   break;
		   
		   case 7: //aproximação Padé parabólico (CRS Padé parabólico expansão em m)
			app_s="Padé tm";
			pade_tm(t0, m0,  h0,  x0,  v0,  RN,  RNIP,  BETA, nh,  dh, nm,  dm,  t);
			sf_warning("Aproximação (%s) RN=%f RNIP=%f BETA=%f",app_s,RN,RNIP,BETA);
		   break;
		   
		   case 8: //aproximação Padé hiperbólico (CRS Padé hiperbólico expansão em h)
			app_s="Padé t2h";
			pade_t2h(t0, m0,  h0,  x0,  v0,  RN,  RNIP,  BETA, nh,  dh, nm,  dm,  t);
			sf_warning("Aproximação (%s) RN=%f RNIP=%f BETA=%f",app_s,RN,RNIP,BETA);
		   break;
		   
		   case 9: //aproximação Padé hiperbólico (CRS Padé hiperbólico expansão em m)
			app_s="Padé t2m";
			pade_t2m(t0, m0,  h0,  x0,  v0,  RN,  RNIP,  BETA, nh,  dh, nm,  dm,  t);
			sf_warning("Aproximação (%s) RN=%f RNIP=%f BETA=%f",app_s,RN,RNIP,BETA);
		   break;
		   
		   case 10: //aproximação Padé Deslocado (CRS Padé hipérbole deslocada expansão em h)
			app_s="Padé tsh";
			pade_tsh(t0, m0,  h0,  x0,  v0,  RN,  RNIP,  BETA, nh,  dh, nm,  dm,  t);
			sf_warning("Aproximação (%s) RN=%f RNIP=%f BETA=%f",app_s,RN,RNIP,BETA);
		   break;
		   
		   case 11: //aproximação Padé Deslocado (CRS Padé hipérbole deslocada expansão em m)
			app_s="Padé tsm";
			pade_tsm(t0, m0,  h0,  x0,  v0,  RN,  RNIP,  BETA, nh,  dh, nm,  dm,  t);
			sf_warning("Aproximação (%s) RN=%f RNIP=%f BETA=%f",app_s,RN,RNIP,BETA);
		   break;

		   default:
			 sf_error("Opção app=%i Não disponível",app);
	}

	/*********************************************************************/

	/* Superfície de tempo de trânsito aproximada */
	sf_floatwrite(t[0],nm*nh,out);
	
	exit(0);
}
