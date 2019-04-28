/* Versão 1.0 - Biblioteca de funções de MvfsaSemb.c 

Programador: Rodolfo A. C. Neves 14/04/2019

Email:  rodolfo_profissional@hotmail.com  

Acesse conteúdo exclusivo e tire dúvidas no nosso site:
	http://www.dirackslounge.online
*/

#define hMAX 150
#define mMAX 200
#define PI 3.14159
#include <stdio.h> // biblioteca padrão define operações de entrada/saída
#include <stdlib.h> // biblioteca padrão define alocação de memória, controle de processos
#include <math.h> // biblioteca de funções matemáticas como exp()
#include <stdlib.h> // necessário para gerar número aleatório com: rand() e srand()
#include <time.h> // necessário para gerar semente de numero aleatório: time()
#include <string.h> // necessário para lidar com strings e strcat()
#include <rsf.h> // biblioteca padrão do madagascar
/*^*/

/********************[Funções]**********************************/

float encontrarT0(float ***t, int im0, float dt, float nt){
/*< Dado o índice de m0, encontrar o tempo t0 (será um pico de amplitude no vetor t[im0][0][ti]) >*/

	float amp=0;
	int it;
	float maior=0;
	float t0=0;

	/* Procure por um pico de amplitude, ele marca a primeira chegada */
	for(it=0; it<nt;it++){

		amp=t[im0][0][it];

		//if(amp>0.000)sf_warning("amp(%f)=%f",it*dt,amp);
		
		if(amp>maior){
			maior=amp;
			t0=(it-5)*dt;

			sf_warning("amp(%f)=%f",it*dt,amp);
			
		}

	}

	return t0;

}

void f_vfsa_aviso(int app){
/*< Informar usuário sobre a aproximação de tempo de trânsito CRS escolhida >*/
		
		switch (app)
		{
		   case 1: //aproximação fomel (CRS NÃO hiperbólico)
			 sf_warning("Aproximação Fomel (CRS NÃO hiperbólico)");
		   break;

		   case 2: //aproximação jager (CRS hiperbólico)
			 sf_warning("Aproximação Jager (CRS hiperbólico)");
		   break;
		   
		   case 3: //aproximação germam-t (CRS quarta ordem - parabólico)
			 sf_warning("Aproximação germam-t (CRS quarta ordem - parabólico)");
		   break;
		   
		   case 4: //aproximação germam-t2 (CRS quarta ordem - quadrático)
			 sf_warning("Aproximação germam-t2 (CRS quarta ordem - quadrático)");
		   break;
		   
		   case 5: //aproximação germam-tshift (CRS quarta ordem - hipérbole deslocada)
			 sf_warning("Aproximação germam-tshift (CRS quarta ordem - hipérbole deslocada)");
		   break;
		   
		   case 6: //aproximação Padé parabólico (CRS Padé parabólico expansão em h)
			 sf_warning("Aproximação Padé parabólico (CRS Padé parabólico expansão em h)");
		   break;
		   
		   case 7: //aproximação Padé parabólico (CRS Padé parabólico expansão em m)
			 sf_warning("Aproximação Padé parabólico (CRS Padé parabólico expansão em m)");
		   break;
		   
		   case 8: //aproximação Padé hiperbólico (CRS Padé hiperbólico expansão em h)
			 sf_warning("Aproximação Padé hiperbólico (CRS Padé hiperbólico expansão em h)");
		   break;
		   
		   case 9: //aproximação Padé hiperbólico (CRS Padé hiperbólico expansão em m)
			 sf_warning("Aproximação Padé hiperbólico (CRS Padé hiperbólico expansão em m)");
		   break;
		   
		   case 10: //aproximação Padé Deslocado (CRS Padé hipérbole deslocada expansão em h)
			 sf_warning("Aproximação Padé Deslocado (CRS Padé hipérbole deslocada expansão em h)");
		   break;
		   
		   case 11: //aproximação Padé Deslocado (CRS Padé hipérbole deslocada expansão em m)
			 sf_warning("Aproximação Padé Deslocado (CRS Padé hipérbole deslocada expansão em m)");
		   break;

		   default:
			 sf_error("Opção app=%i não disponível!!!", app);
		}
}


float sinal(float s) { 
/*< função sinal >*/

	if(s > 0){
		
		return s = 1;		
	}
	else if (s < 0){
		
		return s=-1;
		
	}
}


float fomel(float dt, float t0, float m0, float h0, float x0, float v0, float R_N, float R_NIP, float BETA, int nh, float dh, int nm, float dm, float ***t){ 
/*< Semblance da aproximação de tempo de trânsito do CRS não hiperbólico (FOMEL; KAZINNIK, 2013) >*/

		float ampt=0.; // amplitude da amostra
		float amp2t=0.; // amplitude da amostra ao quadrado
		float amp=0.; // soma das amplitudes das amostras
		float amp2=0.; // soma das amplitudes da amostras ao quadrado
		int M=0; // contador de amostras
		float semb=0; // semblance
		int im, ih; // índice da amostra no CMP e no meio afastamento
		int im0=(m0/dm); // índice da amostra do CMP central m0
		float m, h; // coordenada do CMP e do meio afastamento
		float teta; //amostra no tempo
		int tetai;

		/* Parâmetros da aproximação de Fomel */
		float a1=(2*sin(BETA))/(v0);		
		float a2=(2*cos(BETA)*cos(BETA)*t0)/(v0*R_N);
		float b2=(2*cos(BETA)*cos(BETA)*t0)/(v0*R_NIP);
		float c1=2*b2+a1*a1-a2;
		float Fd, Fdmenos, Fdmais;
			
		for (im=im0-mMAX; im < im0+mMAX; im++){
			
			m=im*dm+x0; //coordenada do CMP em relação a origem x0
	
			m=m-m0; // distância em relação ao CMP central m0
			
			for(ih=0;ih<hMAX;ih++){
			
				h=ih*dh+h0; // coordenada do meio afastamento
										
				/* Calcular superfície de tempo de trânsito aproximada */
										
				Fd=(t0+a1*m)*(t0+a1*m)+a2*m*m;				
				Fdmenos=(t0+a1*(m-h))*(t0+a1*(m-h))+a2*(m-h)*(m-h);
				Fdmais=(t0+a1*(m+h))*(t0+a1*(m+h))+a2*(m+h)*(m+h);		
								
				teta=sqrt(Fdmenos*Fdmais);
			
				teta=(Fd+c1*h*h+teta);

				teta=sqrt(teta/2);

				tetai=teta/dt; //índice da amostra

				//sf_warning("rn=%g; rnip=%g; beta=%g",R_N,R_NIP,BETA);
				//sf_warning("indice=%i; teta=%f",tetai,teta);

				if(tetai < 0 || tetai > 250*250*314){
					continue;
				}

				/* Restrição para não permitir tempo negativo 
				if ( teta < 0.){
					sf_warning("m=%g; h=%g; t0=%g",m,h,t0);
					sf_warning("rn=%g; rnip=%g; beta=%g",R_N,R_NIP,BETA);
					 sf_error("a1=%g;a2=%g;b2=%g;c1=%g\n;teta=%g;",a1,a2,b2,c1,teta);
				}*/

				/* Restrição para não permitir tempo negativo */
				//if ( teta < 0.) teta=0.;
				
				/* incremento a quantidade de amostras */
				M++;
				
				/* se o tempo de trânsito da aproximação teta(h,m) está próximo do tempo
				   de trânsito modelado t(h,m) considerar essa amostra no cálculo do 
				   semblance */
				//if (fabs(t[im][ih] - teta) <= 0.2) {	
					
					/* as amplitudes ampt são dadas simulando a fonte por um
					   pulso ricker
					   distância da amostra em relação ao centro do pulso */
					//teta=teta-t[im][ih];
					
					/* amplitude da amostra */
					ampt=t[im][ih][tetai];//(1-teta*teta*0.5)*exp(-8000*(teta*teta)*0.25);

					/*sf_warning("amp=%f",ampt);*/
					
					/* some as amplitude na variável amp */
					amp=amp+ampt;
					
					/* elevar ao quadrado e armazenar em amp2 */
					amp2t=ampt*ampt;
					amp2=amp2+amp2t; 
				//}
				
			}
		
	}		

		//sf_error("fim!");

		/* restrição: garantir que não haja divisão por zero */
		if(amp2==0 || amp==0)		
		       return semb=0;
		else
		  return semb=(amp*amp)/(M*amp2);
		
}


float jager(float t0, float m0, float h0, float x0, float v0, float R_N, float R_NIP, float BETA, int nh, float dh,  int nm,  float dm, float **t) {
/*< Semblance da aproximação do CRS hiperbólico (JAGER et al., 2001) >*/	
		
		float ampt=0.; // amplitude da amostra
		float amp2t=0.; // amplitude da amostra ao quadrado
		float amp=0.; // soma das amplitudes das amostras
		float amp2=0.; // soma das amplitudes da amostras ao quadrado
		int M=0; // contador de amostras
		float semb=0; // semblance
		int im, ih; // índice da amostra no CMP e no meio afastamento
		int im0=m0/dm; // índice da amostra do CMP central m0
		float m, h; // coordenada do CMP e do meio afastamento
		float teta; //amostra no tempo


		/* Parâmetros da aproximação de Jager */
		float a1=(2*sin(BETA))/v0;		
		float a2=(2*cos(BETA)*cos(BETA)*t0)/(v0*R_N);
		float b2=(2*cos(BETA)*cos(BETA)*t0)/(v0*R_NIP);
		float Fd;
				
			
		for (im=im0-mMAX; im < im0+mMAX; im++){
			
			m=im*dm+x0; //coordenada do CMP em relação a origem x0
	
			m=m-m0; // distância em relação ao CMP central m0
			
			for(ih=0;ih<hMAX;ih++){
			
				h=ih*dh+h0; // coordenada do meio afastamento
										
				/* calcular teta (superfície aproximada com a fórmula de jager) */
										
				Fd=(t0+a1*m)*(t0+a1*m)+a2*m*m;
					
				teta=(Fd+b2*h*h);
				teta=sqrt(teta);

				/* Restrição para não permitir tempo negativo */
				if ( teta < 0.){
					sf_warning("m=%g; h=%g; t0=%g",m,h,t0);
					sf_warning("rn=%g; rnip=%g; beta=%g",R_N,R_NIP,BETA);
					 sf_error("a1=%g;a2=%g;b2=%g\n;teta=%g;",a1,a2,b2,teta);
				}
				
				/* Quantidade de amostras */
				M++;
				
				/* se o tempo de trânsito da aproximação teta(h,m) está próximo do tempo
				   de trânsito modelado t(h,m) considero essa amostra no cálculo do 
				   semblance */
				if (fabs(t[im][ih] - teta) <= 0.2) {	
					
					/* as amplitudes ampt são dadas simulando a fonte por um
					   pulso ricker
					   distância da amostra em relação ao centro do pulso */
					teta=teta-t[im][ih];
					
					/* amplitude da amostra */
					ampt=(1-teta*teta*0.5)*exp(-8000*(teta*teta)*0.25);
					
					/* some as amplitudes na variável amp */
					amp=amp+ampt;
					
					/* eleve ao quadrado e armazene em amp2 */
					amp2t=ampt*ampt;
					amp2=amp2+amp2t; 
				}
				
			}
		
	}
		/* restrição: garantir que não haja divisão por zero */
		if(amp2==0 || amp==0)		
		       return semb=0;
		else
		  return semb=(amp*amp)/(M*amp2);
		
}



float germam_t(float t0, float m0, float h0, float x0, float v0, float R_N, float R_NIP, float BETA, int nh, float dh,  int nm,  float dm, float **t) {
/*< Semblance da aproximação do CRS Quarta ordem parabólico (HöCHT, 2002) >*/	

		float ampt=0.; // amplitude da amostra
		float amp2t=0.; // amplitude da amostra ao quadrado
		float amp=0.; // soma das amplitudes das amostras
		float amp2=0.; // soma das amplitudes da amostras ao quadrado
		int M=0; // contador de amostras
		float semb=0; // semblance
		int im, ih; // índice da amostra no CMP e no meio afastamento
		int im0=m0/dm; // índice da amostra do CMP central m0
		float m, h; // coordenada do CMP e do meio afastamento
		float teta; //amostra no tempo

		/* Parâmetros da aproximação de germam_t */
		float s1=(2*sin(BETA)/v0);
		float s2=((cos(BETA)*cos(BETA))/(v0*R_N));
		float s3=(cos(BETA)*cos(BETA)/(v0*R_NIP));	
		float s4=((sin(BETA)*cos(BETA)*cos(BETA))/(v0*R_N*R_N));
		float s5=(sin(BETA)*cos(BETA)*cos(BETA)/(v0*R_NIP*R_NIP*R_N))*(2*R_NIP+R_N);
		float s6=(cos(BETA)*cos(BETA)/(2*v0*R_NIP*R_NIP*R_NIP*R_N*R_N))*(R_NIP*R_NIP*(8*cos(BETA)*cos(BETA)-6)+R_NIP*R_N*(5*cos(BETA)*cos(BETA)-4)-2*R_N*R_N*sin(BETA)*sin(BETA));
		float s7=(cos(BETA)*cos(BETA)*(5*cos(BETA)*cos(BETA)-4))/(4*v0*R_N*R_N*R_N);
		float s8=cos(BETA)*cos(BETA)*(4*R_NIP*sin(BETA)*sin(BETA)-R_N*cos(BETA)*cos(BETA))/(4*v0*R_NIP*R_NIP*R_NIP*R_N);
		float g1, g2, g3, g4,g5;
				
			
		for (im=im0-mMAX; im < im0+mMAX; im++){
			
			m=im*dm+x0; //coordenada do CMP em relação a origem x0
	
			m=m-m0; // distância em relação ao CMP central m0
			
			for(ih=0;ih<hMAX;ih++){
			
				h=ih*dh+h0; // coordenada do meio afastamento
										
				/* APROXIMAÇÃO PARABÓLICA SRC QUARTA ORDEM germam_t */
				g1=t0+m*s1+s2*m*m+h*h*s3-s4*m*m*m;
				g2=-h*h*m*s5;
				g3=-h*h*m*m*s6;
				g4=-m*m*m*m*s7;
				g5=h*h*h*h*s8;
				teta=g1+g2+g3+g4+g5;

				/* Restrição para não permitir tempo negativo 
				if ( teta < 0.){
					sf_warning("m=%g; h=%g; t0=%g",m,h,t0);
					sf_warning("rn=%g; rnip=%g; beta=%g",R_N,R_NIP,BETA);
					 sf_error("g1=%g;g2=%g;g3=%g\n;g4=%g;g5=%g;teta=%g;",g1,g2,g3,g4,g5,teta);
				}*/
				
				/* Incremente a Quantidade de amostras */
				M++;
				
				/* se o tempo de trânsito da aproximação teta(h,m) está próximo do tempo
				   de trânsito modelado t(h,m) considero essa amostra no cálculo do 
				   semblance */
				if (fabs(t[im][ih] - teta) <= 0.2) {	
					
					/* as amplitudes ampt são dadas simulando a fonte por um
					   pulso ricker					
					   distância da amostra em relação ao centro do pulso */
					teta=teta-t[im][ih];
					
					/* amplitude da amostra */
					ampt=(1-teta*teta*0.5)*exp(-8000*(teta*teta)*0.25);
					
					/* some as amplitudes na variável amp */
					amp=amp+ampt;
					
					/* eleve ao quadrado e armazene em amp2 */
					amp2t=ampt*ampt;
					amp2=amp2+amp2t; 
				}
				
			}
		
	}
		/* restrição: garantir que não haja divisão por zero */
		if(amp2==0 || amp==0)		
		       return semb=0;
		else
		  return semb=(amp*amp)/(M*amp2);
		
}


float germam_t2(float t0, float m0, float h0, float x0, float v0, float R_N, float R_NIP, float BETA, int nh, float dh,  int nm,  float dm, float **t) {
/*< Semblance da aproximação do CRS quarta ordem hiperbólico (HöCHT, 2002) >*/	

		float ampt=0.; // amplitude da amostra
		float amp2t=0.; // amplitude da amostra ao quadrado
		float amp=0.; // soma das amplitudes das amostras
		float amp2=0.; // soma das amplitudes da amostras ao quadrado
		int M=0; // contador de amostras
		float semb=0; // semblance
		int im, ih; // índice da amostra no CMP e no meio afastamento
		int im0=m0/dm; // índice da amostra do CMP central m0
		float m, h; // coordenada do CMP e do meio afastamento
		float teta; //amostra no tempo

		/* Parâmetros da aproximação de germam_t2 */
		float s1=t0*t0;
		float s2=((4*t0*sin(BETA))/v0);
		float s3=(2*(v0*t0*cos(BETA)*cos(BETA)+2*R_N*sin(BETA)*sin(BETA))/(v0*v0*R_N));
		float s4=((2*t0*cos(BETA)*cos(BETA))/(v0*R_NIP));
		float s5=((2*sin(BETA)*cos(BETA)*cos(BETA)*(2*R_N-v0*t0))/(v0*v0*R_N*R_N));
		float s6=(2*sin(BETA)*cos(BETA)*cos(BETA)*(2*R_NIP*R_N-2*v0*t0*R_NIP-v0*t0*R_N))/(v0*v0*R_NIP*R_NIP*R_N);
		float s7=(cos(BETA)*cos(BETA)*(R_N*(10*cos(BETA)*cos(BETA)-8)+v0*t0*(4-5*cos(BETA)*cos(BETA))))/(2*v0*v0*R_N*R_N*R_N);
		float s8=(cos(BETA)*cos(BETA)/(v0*v0*R_NIP*R_NIP*R_NIP*R_N*R_N))*(v0*t0*R_NIP*R_NIP*(6-8*cos(BETA)*cos(BETA))+v0*t0*R_NIP*R_N*(4-5*cos(BETA)*cos(BETA))+2*v0*t0*R_N*R_N*sin(BETA)*sin(BETA)-4*R_NIP*R_N*R_N*sin(BETA)*sin(BETA)+R_NIP*R_NIP*R_N*(10*cos(BETA)-8));
		float s9=(cos(BETA)*cos(BETA)*(4*v0*t0*R_NIP*sin(BETA)*sin(BETA)-v0*t0*R_N*cos(BETA)*cos(BETA)+2*R_NIP*R_N*cos(BETA)*cos(BETA)))/(2*v0*v0*R_NIP*R_NIP*R_NIP*R_N);
		float g1, g2, g3, g4, g5, g6;
						
		for (im=im0-mMAX; im < im0+mMAX; im++){
			
			m=im*dm+x0; //coordenada do CMP em relação a origem x0
	
			m=m-m0; // distância em relação ao CMP central m0
			
			for(ih=0;ih<hMAX;ih++){
			
				h=ih*dh+h0; // coordenada do meio afastamento
										
				/* APROXIMAÇÃO QUADRÁTICA SRC QUARTA ORDEM germam t2 */
				g1=s1+s2*m+s3*(m*m)+s4*h*h;
				g2=s5*m*m*m;
				g3=s6*(m)*h*h;
				g4=s7*m*m*m*m;
				g5=s8*(m*m*h*h);
				g6=s9*h*h*h*h;
				teta=g1+g2+g3+g4+g5+g6;
				teta=sqrt(teta);
				
				/* Restrição para não permitir tempo negativo */
				if ( teta < 0.) teta=0.;

				/* Incremente a quantidade de amostras */
				M++;
				
				/* se o tempo de trânsito da aproximação teta(h,m) está próximo do tempo
				   de trânsito modelado t(h,m) considero essa amostra no cálculo do 
				   semblance */
				if (fabs(t[im][ih] - teta) <= 0.2) {	
					
					/* as amplitudes ampt são dadas simulando a fonte por um
					   pulso ricker
					   distância da amostra em relação ao centro do pulso */
					teta=teta-t[im][ih];
					
					/* amplitude da amostra */
					ampt=(1-teta*teta*0.5)*exp(-8000*(teta*teta)*0.25);
					
					/* some as amplitude na variável amp */
					amp=amp+ampt;
					
					/* eleve ao quadrado e armazene em amp2 */
					amp2t=ampt*ampt;
					amp2=amp2+amp2t; 
				}
				
			}
		
	}
	
		/* restrição: garantir que não haja divisão por zero */
		if(amp2==0 || amp==0)		
		       return semb=0;
		else
		  return semb=(amp*amp)/(M*amp2);
		
}


float germam_tshift(float t0, float m0, float h0, float x0, float v0, float R_N, float R_NIP, float BETA, int nh, float dh,  int nm,  float dm, float **t){
/*< Semblance da aproximação do CRS quarta ordem deslocado (HöCHT, 2002) >*/	

		float ampt=0.; // amplitude da amostra
		float amp2t=0.; // amplitude da amostra ao quadrado
		float amp=0.; // soma das amplitudes das amostras
		float amp2=0.; // soma das amplitudes da amostras ao quadrado
		int M=0; // contador de amostras
		float semb=0; // semblance
		int im, ih; // índice da amostra no CMP e no meio afastamento
		int im0=m0/dm; // índice da amostra do CMP central m0
		float m, h; // coordenada do CMP e do meio afastamento
		float teta; //amostra no tempo

		/* Parâmetros da aproximação de germam_tshift */	
		float s1=(2*R_NIP/v0)*(2*R_NIP/v0);
		float s2=(8*R_NIP*sin(BETA)/(v0*v0));
		float s3=4*((R_NIP*cos(BETA)*cos(BETA)+R_N*sin(BETA)*sin(BETA))/(v0*v0*R_N));
		float s4=(4*cos(BETA)*cos(BETA)/(v0*v0));
		float s5=((4*sin(BETA)*cos(BETA)*cos(BETA)*(R_N-R_NIP))/(v0*v0*R_N*R_N));
		float s6=((8*sin(BETA)*cos(BETA)*cos(BETA))/(v0*v0*R_N));
		float s7=(cos(BETA)*cos(BETA)/(v0*v0*R_N*R_N*R_N))*(R_N*(5*cos(BETA)*cos(BETA)-4)+R_NIP*(4-5*cos(BETA)*cos(BETA)));
		float s8=((4*cos(BETA)*cos(BETA)*(3-4*cos(BETA)*cos(BETA)))/(v0*v0*R_N*R_N));
		float s9=(4*cos(BETA)*cos(BETA)*sin(BETA)*sin(BETA))/(v0*v0*R_NIP*R_N);
		float s10=t0-(2*R_NIP/v0);
		float g1, g2, g3, g4,g5;
						
		for (im=im0-mMAX; im < im0+mMAX; im++){
			
			m=im*dm+x0; //coordenada do CMP em relação a origem x0
	
			m=m-m0; // distância em relação ao CMP central m0
			
			for(ih=0;ih<hMAX;ih++){
			
				h=ih*dh+h0; // coordenada do meio afastamento
										
				/* APROXIMAÇÃO DESLOCADA CRS QUARTA ORDEM germam tshift */               
				g1=s1+m*s2+m*m*s3;
				g2=h*h*s4+m*m*m*s5-m*h*h*s6;
				g3=s7*m*m*m*m;
				g4=m*m*h*h*s8;
				g5=h*h*h*h*s9;
				
				teta=g1+g2+g3+g4+g5;
				teta=sqrt(teta);
				teta=teta+s10;

				/* Restrição para não permitir tempo negativo */
				if ( teta < 0.) teta=0.;		
				
				/* Incremente a quantidade de amostras */
				M++;
				
				/* se o tempo de trânsito da aproximação teta(h,m) está próximo do tempo
				   de trânsito modelado t(h,m) considere essa amostra no cálculo do 
				   semblance */
				if (fabs(t[im][ih] - teta) <= 0.2) {	
					
					/* as amplitudes ampt são dadas simulando a fonte por um
					   pulso ricker
					   distância da amostra em relação ao centro do pulso */
					teta=teta-t[im][ih];
					
					/* amplitude da amostra */
					ampt=(1-teta*teta*0.5)*exp(-8000*(teta*teta)*0.25);
					
					/* some as amplitude na variável amp */
					amp=amp+ampt;
					
					/* eleve ao quadrado e armazeno em amp2 */
					amp2t=ampt*ampt;
					amp2=amp2+amp2t; 
				}
				
			}
		
	}

		/* restrição: garantir que não haja divisão por zero */
		if(amp2==0 || amp==0)		
		       return semb=0;
		else
		  return semb=(amp*amp)/(M*amp2);
		
}


float pade_th(float t0, float m0, float h0, float x0, float v0, float R_N, float R_NIP, float BETA, int nh, float dh,  int nm,  float dm, float **t){
/*< Semblance da aproximação do CRS Padé parabólico expansão em h (NEVES, 2017) >*/

		float ampt=0.; // amplitude da amostra
		float amp2t=0.; // amplitude da amostra ao quadrado
		float amp=0.; // soma das amplitudes das amostras
		float amp2=0.; // soma das amplitudes da amostras ao quadrado
		int M=0; // contador de amostras
		float semb=0; // semblance
		int im, ih; // índice da amostra no CMP e no meio afastamento
		int im0=m0/dm; // índice da amostra do CMP central m0
		float m, h; // coordenada do CMP e do meio afastamento
		float teta; //amostra no tempo	

		/* Parâmetros do CRS padé-t expansão em h */
		float s1=(2*sin(BETA)/v0);
		float s2=((cos(BETA)*cos(BETA))/(v0*R_N));
		float s3=(sin(BETA)*cos(BETA)*cos(BETA))/(v0*R_N*R_N);
		float s4=(cos(BETA)*cos(BETA)*(5*cos(BETA)*cos(BETA)-4))/(4*v0*R_N*R_N*R_N);
		float s5=(cos(BETA)*cos(BETA)/(v0*R_NIP));
		float s6=(cos(BETA)*cos(BETA)/(2*v0*R_NIP*R_NIP*R_NIP*R_N*R_N))*(R_NIP*R_NIP*(8*cos(BETA)*cos(BETA)-6)+R_NIP*R_N*(5*cos(BETA)*cos(BETA)-4)-2*R_N*R_N*sin(BETA)*sin(BETA));
		float s7=(sin(BETA)*cos(BETA)*cos(BETA)/(v0*R_NIP*R_NIP*R_N))*(2*R_NIP+R_N);
		float s8=cos(BETA)*cos(BETA)*(4*R_NIP*sin(BETA)*sin(BETA)-R_N*cos(BETA)*cos(BETA))/(4*v0*R_NIP*R_NIP*R_NIP*R_N);
		float CO,C1,C2;
								
		for (im=im0-mMAX; im < im0+mMAX; im++){
			
			m=im*dm+x0; //coordenada do CMP em relação a origem x0
	
			m=m-m0; // distância em relação ao CMP central m0
			
			for(ih=0;ih<hMAX;ih++){
			
				h=ih*dh+h0; // coordenada do meio afastamento
										
				/* APROXIMAÇÃO PARABÓLICA SRC PADÉ EXPANSÃO EM h */
				CO=t0+m*s1+s2*m*m-s3*m*m*m-m*m*m*m*s4;
				C1=s5-m*m*s6-m*s7;
				C2=s8;
				
				teta=(h*h*C1)/(1+h*h*(-C2/C1));
				teta=CO + teta;

				/* Restrição para não permitir tempo negativo */
				//if ( teta < 0. && teta > 10.) teta=0.;
				
				/* Incremente a quantidade de amostras */
				M++;
				
				/* se o tempo de trânsito da aproximação teta(h,m) está próximo do tempo
				   de trânsito modelado t(h,m) considero essa amostra no cálculo do 
				   semblance */
				if (fabs(t[im][ih] - teta) <= 0.2) {	
					
					/* as amplitudes ampt são dadas simulando a fonte por um
					   pulso ricker
					   distância da amostra em relação ao centro do pulso */
					teta=teta-t[im][ih];
					
					/* amplitude da amostra */
					ampt=(1-teta*teta*0.5)*exp(-8000*(teta*teta)*0.25);
					
					/* some as amplitude na variável amp */
					amp=amp+ampt;
					
					/* eleve ao quadrado e armazene em amp2 */
					amp2t=ampt*ampt;
					amp2=amp2+amp2t; 
				}
				
			}
		
	}
		/* restrição: garantir que não haja divisão por zero */
		if(amp2==0 || amp==0)		
		       return semb=0;
		else
		  return semb=(amp*amp)/(M*amp2);
		
}


float pade_tm(float t0, float m0, float h0, float x0, float v0, float R_N, float R_NIP, float BETA, int nh, float dh,  int nm,  float dm, float **t) {
/*< Semblance da aproximação do CRS Padé parabólico expansão em m (NEVES, 2017) >*/	

		float ampt=0.; // amplitude da amostra
		float amp2t=0.; // amplitude da amostra ao quadrado
		float amp=0.; // soma das amplitudes das amostras
		float amp2=0.; // soma das amplitudes da amostras ao quadrado
		int M=0; // contador de amostras
		float semb=0; // semblance
		int im, ih; // índice da amostra no CMP e no meio afastamento
		int im0=m0/dm; // índice da amostra do CMP central m0
		float m, h; // coordenada do CMP e do meio afastamento
		float teta; //amostra no tempo	

		/* Parâmetros do CRS padé-t em m */
		float s1=(cos(BETA)*cos(BETA)/(v0*R_NIP));
		float s2=((cos(BETA)*cos(BETA)*(4*R_NIP*sin(BETA)*sin(BETA)-R_N*cos(BETA)*cos(BETA)))/(4*v0*R_NIP*R_NIP*R_NIP*R_N));
		float s3=(2*sin(BETA)/v0);
		float s4=((sin(BETA)*cos(BETA)*cos(BETA)*(2*R_NIP+R_N))/(v0*R_NIP*R_NIP*R_N));
		float kapa=R_NIP*R_NIP*(8*cos(BETA)*cos(BETA)-6)+R_NIP*R_N*(5*cos(BETA)*cos(BETA)-4)-2*R_N*R_N*sin(BETA)*sin(BETA);
		float s5=(cos(BETA)*cos(BETA)/(v0*R_N));
		float s6=((cos(BETA)*cos(BETA)/(2*v0*R_NIP*R_NIP*R_NIP*R_N*R_N))*kapa);
		float s7=-((sin(BETA)*cos(BETA)*cos(BETA))/(v0*R_N*R_N));
		float s8=-((cos(BETA)*cos(BETA)*(5*cos(BETA)*cos(BETA)-4))/(4*v0*R_N*R_N*R_N));
		float CO,C1,C2,C3,C4,q1,q2,po,p1,p2,pp,qq;
								
		for (im=im0-mMAX; im < im0+mMAX; im++){
			
			m=im*dm+x0; //coordenada do CMP em relação a origem x0
	
			m=m-m0; // distância em relação ao CMP central m0
			
			for(ih=0;ih<hMAX;ih++){
			
				h=ih*dh+h0; // coordenada do meio afastamento
										
				/* APROXIMAÇÃO PARABÓLICA SRC PADE EXPANSÃO EM m */
        			CO=t0+s1*h*h+s2*h*h*h*h;
				C1=s3-h*h*s4;
				C2=s5-h*h*s6;
				C3=s7;
				C4=s8;
				
				q1=(-C3*C2)/(C2*C2-C1*C2+C1*C3);
				q2=-(C4+C3*q1)/C2;
				po=CO;
				p1=CO*q1+C1;
				p2=CO*q2+C2+C1*q1;
				
				pp=po+p1*m+p2*m*m;
				qq=1+q1*m+q2*m*m;
				
				teta=pp/qq;

				/* Restrição para não permitir tempo negativo */
				//if ( teta < 0.) teta=0.;
			
				/* Incremente a quantidade de amostras */
				M++;
				
				/* se o tempo de trânsito da aproximação teta(h,m) está próximo do tempo
				   de trânsito modelado t(h,m) considero essa amostra no cálculo do 
				   semblance */
				if (fabs(t[im][ih] - teta) <= 0.2) {	
					
					/* as amplitudes ampt são dadas simulando a fonte por um
					   pulso ricker
					   distância da amostra em relação ao centro do pulso */
					teta=teta-t[im][ih];
					
					/* amplitude da amostra */
					ampt=(1-teta*teta*0.5)*exp(-8000*(teta*teta)*0.25);
					
					/* some as amplitude na variável amp */
					amp=amp+ampt;
					
					/* eleve ao quadrado e armazene em amp2 */
					amp2t=ampt*ampt;
					amp2=amp2+amp2t; 
				}
				
			}
		
	}
		/* restrição: garantir que não haja divisão por zero */
		if(amp2==0 || amp==0)		
		       return semb=0;
		else
		  return semb=(amp*amp)/(M*amp2);
		
}


float pade_t2h(float t0, float m0, float h0, float x0, float v0, float R_N, float R_NIP, float BETA, int nh, float dh,  int nm,  float dm, float **t) {
/*< Semblance da aproximação do CRS Padé hiperbólico expansão em h (NEVES, 2017) >*/

		float ampt=0.; // amplitude da amostra
		float amp2t=0.; // amplitude da amostra ao quadrado
		float amp=0.; // soma das amplitudes das amostras
		float amp2=0.; // soma das amplitudes da amostras ao quadrado
		int M=0; // contador de amostras
		float semb=0; // semblance
		int im, ih; // índice da amostra no CMP e no meio afastamento
		int im0=m0/dm; // índice da amostra do CMP central m0
		float m, h; // coordenada do CMP e do meio afastamento
		float teta; //amostra no tempo	
	
		/* Parâmetros do CRS padé-t2 em h */
		double s1=t0*t0;
		double s2=((4*t0*sin(BETA))/v0);
		double s3=(2*(v0*t0*cos(BETA)*cos(BETA)+2*R_N*sin(BETA)*sin(BETA))/(v0*v0*R_N));
		double s4=((2*sin(BETA)*cos(BETA)*cos(BETA)*(2*R_N-v0*t0))/(v0*v0*R_N*R_N));
		double s5=(cos(BETA)*cos(BETA)*(R_N*(10*cos(BETA)*cos(BETA)-8)+v0*t0*(4-5*cos(BETA)*cos(BETA))))/(2*v0*v0*R_N*R_N*R_N);
		double s6=((2*t0*cos(BETA)*cos(BETA))/(v0*R_NIP));
		double s7=(2*sin(BETA)*cos(BETA)*cos(BETA)*(2*R_NIP*R_N-2*v0*t0*R_NIP-v0*t0*R_N))/(v0*v0*R_NIP*R_NIP*R_N);
		double s8=(cos(BETA)*cos(BETA)/(v0*v0*R_NIP*R_NIP*R_NIP*R_N*R_N))*(v0*t0*R_NIP*R_NIP*(6-8*cos(BETA)*cos(BETA))+v0*t0*R_NIP*R_N*(4-5*cos(BETA)*cos(BETA))+(2*v0*t0*R_N*R_N*sin(BETA)*sin(BETA)-4*R_NIP*R_N*R_N*sin(BETA)*sin(BETA)+R_NIP*R_NIP*R_N*(10*cos(BETA)*cos(BETA)-8)));
		double s9=((cos(BETA)*cos(BETA))/(2*v0*v0*R_NIP*R_NIP*R_NIP*R_N))*(4*v0*t0*R_NIP*sin(BETA)*sin(BETA)-v0*t0*R_N*cos(BETA)*cos(BETA)+2*R_NIP*R_N*cos(BETA)*cos(BETA));
		
		double k3=-(((4*R_N*v0*t0*R_NIP*sin(BETA)*sin(BETA)-v0*t0*R_N*R_N*cos(BETA)*cos(BETA)+2*R_NIP*R_N*cos(BETA)*cos(BETA)))/2);
		double k1=2*v0*t0*R_NIP*R_NIP*R_N;
		double k2=2*R_N*R_NIP*sin(BETA)*((2*R_NIP*R_N-2*v0*t0*R_NIP-v0*t0*R_N));
		double k4=((v0*t0*R_NIP*R_NIP*(6-8*cos(BETA)*cos(BETA))+v0*t0*R_NIP*R_N*(4-5*cos(BETA)*cos(BETA))));
		double CO,C1,C3;
							
		for (im=im0-mMAX; im < im0+mMAX; im++){
			
			m=im*dm+x0; //coordenada do CMP em relação a origem x0
	
			m=m-m0; // distância em relação ao CMP central m0
			
			for(ih=0;ih<hMAX;ih++){
			
				h=ih*dh+h0; // coordenada do meio afastamento
										
				/* APROXIMAÇÃO QUADRÁTICA SRC-PADÉ EXPANSÃO EM h */
				CO=s1+s2*m+s3*m*m+s4*m*m*m+s5*m*m*m*m;
				C1=s6+s7*m+s8*m*m;
				C3=-s9/C1;//k3/(k1+k2*m+k4*m*m);
									
				teta=(h*h*C1)/(1+h*h*(C3));

				//sf_warning("\ndh=%f;\ndm=%f;\nh=%f;m=%f;\nteta=%f\n",dh,dm,h,m,teta);
				teta=CO+teta;
				//sf_warning("\nh=%f;m=%f;teta=%f\n",h,m,teta);

				/* Restrição para não permitir tempo negativo */
				//if ( teta < 0.){
				
					/*sf_warning("\ndh=%f;\ndm=%f;\nh=%f;m=%f",dh,dm,h,m);
					sf_warning("\ns1=%f;\ns2=%f;\ns3=%f;\ns4=%f;\n",s1,s2,s3,s4);
					sf_warning("\ns5=%f;\ns6=%f;\ns7=%f;\ns8=%f;\n",s5,s6,s7,s8);*/
					//sf_warning("\nCO=%f;\nC1=%f;\nC3=%f;\n",CO,C1,C3);
					//sf_error("\nk1=%f;\nk2=%f\nk3=%f\nk4=%f\n",k1,k2,k3,k4);
					//sf_warning("\nh=%f;m=%f;teta=%f\n",h,m,teta);
				//}

				teta=sqrt(teta);
				//sf_error("\ndh=%f;\ndm=%f;\nh=%f;m=%f;\nteta=%f\n",dh,dm,h,m,teta);

				/*teta=CO+C1*h*h+(C3*h*h*h*h)/(1-C3*h*h/C1);
				teta=sqrt(teta);*/ 
				//t[im][ih]=teta;
				
				/* Incremente a quantidade de amostras */
				M++;
				
				/* se o tempo de trânsito da aproximação teta(h,m) está próximo do tempo
				   de trânsito modelado t(h,m) considero essa amostra no cálculo do 
				   semblance */
				if (fabs(t[im][ih] - teta) <= 0.2) {	
					
					/* as amplitudes ampt são dadas simulando a fonte por um
					   pulso ricker	
					   distância da amostra em relação ao centro do pulso */
					teta=teta-t[im][ih];
					
					/* amplitude da amostra */
					ampt=(1-teta*teta*0.5)*exp(-8000*(teta*teta)*0.25);
					
					/* some as amplitude na variável amp */
					amp=amp+ampt;
					
					/* eleve ao quadrado e armazene em amp2 */
					amp2t=ampt*ampt;
					amp2=amp2+amp2t; 
				}
				
			}
		
	}
		/* restrição: garantir que não haja divisão por zero */
		if(amp2==0 || amp==0)		
		       return semb=0;
		else
		  return semb=(amp*amp)/(M*amp2);
		
}


float pade_t2m(float t0, float m0, float h0, float x0, float v0, float R_N, float R_NIP, float BETA, int nh, float dh,  int nm,  float dm, float **t) {
/*< Semblance da aproximação do CRS Padé hiperbólico expansão em m (NEVES, 2017)>*/
	
		float ampt=0.; // amplitude da amostra
		float amp2t=0.; // amplitude da amostra ao quadrado
		float amp=0.; // soma das amplitudes das amostras
		float amp2=0.; // soma das amplitudes da amostras ao quadrado
		int M=0; // contador de amostras
		float semb=0; // semblance
		int im, ih; // índice da amostra no CMP e no meio afastamento
		int im0=m0/dm; // índice da amostra do CMP central m0
		float m, h; // coordenada do CMP e do meio afastamento
		float teta; //amostra no tempo	

		/* Parâmetros do CRS padé-t2 em h */
		float s1=t0*t0;
		float s2=((2*t0*cos(BETA)*cos(BETA))/(v0*R_NIP));
		float s3=(cos(BETA)*cos(BETA)*(4*v0*t0*R_NIP*sin(BETA)*sin(BETA)-v0*t0*R_N*cos(BETA)*cos(BETA)+2*R_NIP*R_N*cos(BETA)*cos(BETA)))/(2*v0*v0*R_NIP*R_NIP*R_NIP*R_N);
		float s4=(((4*t0*sin(BETA)))/v0)+(2*sin(BETA)*cos(BETA)*cos(BETA)*(2*R_NIP*R_N-2*v0*t0*R_NIP-v0*t0*R_N))/(v0*v0*R_NIP*R_NIP*R_N);
		float s5=((2*(v0*t0*cos(BETA)*cos(BETA)+2*R_N*sin(BETA)*sin(BETA)))/(v0*v0*R_N));
		float s6=(cos(BETA)*cos(BETA)/(v0*v0*R_NIP*R_NIP*R_NIP*R_N*R_N))*(v0*t0*R_NIP*R_NIP*(6-8*cos(BETA)*cos(BETA))+v0*t0*R_NIP*R_N*(4-5*cos(BETA)*cos(BETA))+2*v0*t0*R_N*R_N*sin(BETA)*sin(BETA)-4*R_NIP*R_N*R_N*sin(BETA)*sin(BETA)+R_NIP*R_NIP*R_N*(10*cos(BETA)-8));
		float s7=((2*sin(BETA)*cos(BETA)*cos(BETA)*(2*R_N-v0*t0))/(v0*v0*R_N*R_N));
		float s8=(cos(BETA)*cos(BETA)*(R_N*(10*cos(BETA)*cos(BETA)-8)+v0*t0*(4-5*cos(BETA)*cos(BETA))))/(2*v0*v0*R_N*R_N*R_N);
		float CO,C1,C2,C3,C4,q1,q2,po,p1,p2,pp,qq;
					
		for (im=im0-mMAX; im < im0+mMAX; im++){
			
			m=im*dm+x0; //coordenada do CMP em relação a origem x0
	
			m=m-m0; // distância em relação ao CMP central m0
			
			for(ih=0;ih<hMAX;ih++){
			
				h=ih*dh+h0; // coordenada do meio afastamento
										
				/* APROXIMAÇÃO QUADRÁTICA SRC PADE EXAPANSÃO EM d */
				CO=s1+h*h*s2+h*h*h*h*s3;
				C1=s4*(h*h);
				C2=s5+s6*(h*h);
				C3=s7;
				C4=s8;
				
				q1=(-C3*C2)/(C2*C2-C1*C2+C1*C3);
				q2=-(C4+C3*q1)/C2;
				po=CO;
				p1=CO*q1+C1;
				p2=CO*q2+C2+C1*q1;
				
				pp=po+p1*m+p2*m*m;
				qq=1+q1*m+q2*m*m;
				
				teta=pp/qq;
				teta=sqrt(teta);
			
				/* Restrição para não permitir tempo negativo */
				//if ( teta < 0.) teta=0.;
				
				/* Incremente a quantidade de amostras */
				M++;
				
				/* se o tempo de trânsito da aproximação teta(h,m) está próximo do tempo
				   de trânsito modelado t(h,m) considero essa amostra no cálculo do 
				   semblance */
				if (fabs(t[im][ih] - teta) <= 0.2) {	
					
					/* as amplitudes ampt são dadas simulando a fonte por um
					   pulso ricker
					   distância da amostra em relação ao centro do pulso */
					teta=teta-t[im][ih];
					
					/* amplitude da amostra */
					ampt=(1-teta*teta*0.5)*exp(-8000*(teta*teta)*0.25);
					
					/* some as amplitude na variável amp */
					amp=amp+ampt;
					
					/* eleve ao quadrado e armazene em amp2 */
					amp2t=ampt*ampt;
					amp2=amp2+amp2t; 
				}
				
			}
		
	}
		/* restrição: garantir que não haja divisão por zero */
		if(amp2==0 || amp==0)		
		       return semb=0;
		else
		  return semb=(amp*amp)/(M*amp2);
		
}


float pade_tsh(float t0, float m0, float h0, float x0, float v0, float R_N, float R_NIP, float BETA, int nh, float dh,  int nm,  float dm, float **t){
/*< Semblance da aproximação do CRS Padé deslocada expansão em h (NEVES, 2017) >*/

		float ampt=0.; // amplitude da amostra
		float amp2t=0.; // amplitude da amostra ao quadrado
		float amp=0.; // soma das amplitudes das amostras
		float amp2=0.; // soma das amplitudes da amostras ao quadrado
		int M=0; // contador de amostras
		float semb=0; // semblance
		int im, ih; // índice da amostra no CMP e no meio afastamento
		int im0=m0/dm; // índice da amostra do CMP central m0
		float m, h; // coordenada do CMP e do meio afastamento
		float teta; //amostra no tempo	
	
		/* Parâmetros do CRS padé-tshift em h */
		float s1=(2*R_NIP/v0)*(2*R_NIP/v0);	
		float s2=((8*R_NIP*sin(BETA))/(v0*v0));
		float s3=((R_NIP*cos(BETA)*cos(BETA)+R_N*sin(BETA)*sin(BETA))/(v0*v0*R_N));
		float s4=((4*sin(BETA)*cos(BETA)*cos(BETA)*(R_N-R_NIP))/(v0*v0*R_N*R_N));
		float s5=(cos(BETA)*cos(BETA)/(v0*v0*R_N*R_N*R_N))*(R_N*(5*cos(BETA)*cos(BETA)-4)+R_NIP*(4-5*cos(BETA)*cos(BETA)));
		float s6=((4*cos(BETA)*cos(BETA))/(v0*v0));
		float s7=((8*sin(BETA)*cos(BETA)*cos(BETA))/(v0*v0*R_N));
		float s8=((4*cos(BETA)*cos(BETA)*(3-4*cos(BETA)*cos(BETA)))/(v0*v0*R_N*R_N));
		float s9=(4*cos(BETA)*cos(BETA)*sin(BETA)*sin(BETA))/(v0*v0*R_NIP*R_N);
		float s10=t0-(2*R_NIP/v0);
		float CO,C1,C2;
			
		for (im=im0-mMAX; im < im0+mMAX; im++){
			
			m=im*dm+x0; //coordenada do CMP em relação a origem x0
	
			m=m-m0; // distância em relação ao CMP central m0
			
			for(ih=0;ih<hMAX;ih++){
			
				h=ih*dh+h0; // coordenada do meio afastamento
										
			/* APROXIMAÇÃO DESLOCADA SRC-PADÉ EXPANSÃO EM h */
				CO=s1+m*s2+4*m*m*s3+m*m*m*s4+s5*m*m*m*m;
				C1=s6-m*s7+m*m*s8;
				C2=s9;
				
				
				teta=(h*h*C1)/(1+h*h*(-C2/C1));
				teta=CO + teta;
				teta=sqrt(teta);
				teta=teta+s10;
			
				/* Restrição para não permitir tempo negativo */
				//if ( teta < 0.) teta=0.;
				
				/* Incremente a quantidade de amostras */
				M++;
				
				/* se o tempo de trânsito da aproximação teta(h,m) está próximo do tempo
				   de trânsito modelado t(h,m) considero essa amostra no cálculo do 
				   semblance */
				if (fabs(t[im][ih] - teta) <= 0.2) {	
					
					/* as amplitudes ampt são dadas simulando a fonte por um
					   pulso ricker	
					   distância da amostra em relação ao centro do pulso */
					teta=teta-t[im][ih];
					
					/* amplitude da amostra */
					ampt=(1-teta*teta*0.5)*exp(-8000*(teta*teta)*0.25);
					
					/* some as amplitudes na variável amp */
					amp=amp+ampt;
					
					/* elevo ao quadrado e armazeno em amp2 */
					amp2t=ampt*ampt;
					amp2=amp2+amp2t; 
				}
				
			}
		
	}
		/* restrição: garantir que não haja divisão por zero */
		if(amp2==0 || amp==0)		
		       return semb=0;
		else
		  return semb=(amp*amp)/(M*amp2);
		
}


float pade_tsm(float t0, float m0, float h0, float x0, float v0, float R_N, float R_NIP, float BETA, int nh, float dh,  int nm,  float dm, float **t){
/*< Semblance da aproximação do CRS Padé deslocada expansão em m (NEVES, 2017) >*/

		float ampt=0.; // amplitude da amostra
		float amp2t=0.; // amplitude da amostra ao quadrado
		float amp=0.; // soma das amplitudes das amostras
		float amp2=0.; // soma das amplitudes da amostras ao quadrado
		int M=0; // contador de amostras
		float semb=0; // semblance
		int im, ih; // índice da amostra no CMP e no meio afastamento
		int im0=m0/dm; // índice da amostra do CMP central m0
		float m, h; // coordenada do CMP e do meio afastamento
		float teta; //amostra no tempo	
	
		/* Parâmetros do CRS padé-tshift em m */
		float s1=(2*R_NIP/v0)*(2*R_NIP/v0);
		float s2=((4*cos(BETA)*cos(BETA))/(v0*v0));
		float s3=(4*cos(BETA)*cos(BETA)*sin(BETA)*sin(BETA))/(v0*v0*R_NIP*R_N);
		float s4=((8*R_NIP*sin(BETA))/(v0*v0));
		float s5=((8*sin(BETA)*cos(BETA)*cos(BETA))/(v0*v0*R_N));
		float s6=4*((R_NIP*cos(BETA)*cos(BETA)+R_N*sin(BETA)*sin(BETA))/(v0*v0*R_N));
		float s7=((4*cos(BETA)*cos(BETA)*(3-4*cos(BETA)*cos(BETA)))/(v0*v0*R_N*R_N));
		float s8=((4*sin(BETA)*cos(BETA)*cos(BETA)*(R_N-R_NIP))/(v0*v0*R_N*R_N));
		float s9=(cos(BETA)*cos(BETA))*(R_N*(5*cos(BETA)*cos(BETA)-4)+R_NIP*(4-5*cos(BETA)*cos(BETA)));
		float s10=(v0*v0*R_N*R_N*R_N);
		float s11=t0-(2*R_NIP/v0);
		float CO,C1,C2,C3,C4,q1,q2,po,p1,p2,pp,qq;
					
		for (im=im0-mMAX; im < im0+mMAX; im++){
			
			m=im*dm+x0; //coordenada do CMP em relação a origem x0
	
			m=m-m0; // distância em relação ao CMP central m0
			
			for(ih=0;ih<hMAX;ih++){
			
				h=ih*dh+h0; // coordenada do meio afastamento
										
			/* APROXIMAÇÃO DESLOCADA SRC PADE EXPANSÃO EM d */
				CO=s1+h*h*s2+h*h*h*h*s3;
				C1=s4-h*h*s5;
				C2=s6+h*h*s7;
				C3=s8;
				C4=s9;
				C4=C4/s10;
				
				q1=(-C3*C2)/(C2*C2-C1*C2+C1*C3);
				q2=-(C4+C3*q1)/C2;
				po=CO;
				p1=CO*q1+C1;
				p2=CO*q2+C2+C1*q1;
				
				pp=po+p1*m+p2*m*m;
				qq=1+q1*m+q2*m*m;
				
				teta=pp/qq;
				teta=sqrt(teta);
				teta=teta+s11;
			
				/* Restrição para não permitir tempo negativo */
				//if ( teta < 0.) teta=0.;
				
				/* Incremente a quantidade de amostras */
				M++;
				
				/* se o tempo de trânsito da aproximação teta(h,m) está próximo do tempo
				   de trânsito modelado t(h,m) considero essa amostra no cálculo do 
				   semblance */
				if (fabs(t[im][ih] - teta) <= 0.2) {	
					
					/* as amplitudes ampt são dadas simulando a fonte por um
					   pulso ricker
					   distância da amostra em relação ao centro do pulso */
					teta=teta-t[im][ih];
					
					/* amplitude da amostra */
					ampt=(1-teta*teta*0.5)*exp(-8000*(teta*teta)*0.25);
					
					/* some as amplitudes na variável amp */
					amp=amp+ampt;
					
					/* eleve ao quadrado e armazene em amp2 */
					amp2t=ampt*ampt;
					amp2=amp2+amp2t; 
				}
				
			}
		
	}
		/* restrição: garantir que não haja divisão por zero */
		if(amp2==0 || amp==0)		
		       return semb=0;
		else
		  return semb=(amp*amp)/(M*amp2);
		
}
