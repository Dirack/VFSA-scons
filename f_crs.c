/* Versão 1.0 - Biblioteca de funções de Mcrs.c 

Programador: Rodolfo A. C. Neves 17/11/2018

Email:  rodolfo_profissional@hotmail.com  

Acesse conteúdo exclusivo e tire dúvidas no nosso site:
	http://www.dirackslounge.online
*/

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
			 sf_warning("Aproximação Padé parabólico (CRS Padé parabólico expansão em m)");
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


void fomel(float t0, float m0, float h0, float x0, float v0, float RN, float RNIP, float BETA, int nh, float dh, int nm, float dm, float **t){ 
/*< Aproximação de tempo de trânsito do CRS não hiperbólico (FOMEL; KAZINNIK, 2013) >*/

		int im, ih; // contadores de laço
		float m; //coordenadas do CMP
		float h; // coordenadas do meio-afastamento
		float teta; //amostra no tempo

		/* Parâmetros da aproximação de Fomel */
		float a1=(2*sin(BETA))/v0;		
		float a2=(2*cos(BETA)*cos(BETA)*t0)/(v0*RN);
		float b2=(2*cos(BETA)*cos(BETA)*t0)/(v0*RNIP);
		float c1=2*b2+a1*a1-a2;	
		float Fd, Fdmenos, Fdmais;	
			
		for (im=0; im < nm; im++){
			
			m=im*dm+x0; //coordenada do cmp

			m=m-m0; // distância em relação ao CMP central
				
			for(ih=0;ih<nh;ih++){
			
				h=ih*dh+h0; // coordenada do meio afastamento
										
				/* Calcular superfície de tempo de trânsito aproximada teta(h,m) */
										
				Fd=(t0+a1*m)*(t0+a1*m)+a2*m*m;				
				Fdmenos=(t0+a1*(m-h))*(t0+a1*(m-h))+a2*(m-h)*(m-h);
				Fdmais=(t0+a1*(m+h))*(t0+a1*(m+h))+a2*(m+h)*(m+h);		
								
				teta=sqrt(Fdmenos*Fdmais);			
				teta=(Fd+c1*h*h+teta);
				teta=sqrt(teta/2);

				t[im][ih]=teta;	
				
			}
		
	}
				
}


void jager(float t0, float m0, float h0, float x0, float v0, float RN, float RNIP, float BETA, int nh, float dh, int nm, float dm, float **t) {
/*< Semblance da aproximação do CRS hiperbólico (JAGER et al., 2001) >*/

		int im, ih; // contadores de laço
		float m; //coordenadas do CMP
		float h; // coordenadas do meio-afastamento
		float teta; //amostra no tempo
	
		/* Parâmetros da aproximação de Jager */
		float a1=(2*sin(BETA))/v0;		
		float a2=(2*cos(BETA)*cos(BETA)*t0)/(v0*RN);
		float b2=(2*cos(BETA)*cos(BETA)*t0)/(v0*RNIP);
		float Fd;	
			
		for (im=0; im < nm; im++){
			
			m=im*dm+x0; //coordenada do cmp

			m=m-m0; // distância em relação ao CMP central
				
			for(ih=0;ih<nh;ih++){
			
				h=ih*dh+h0; // coordenada do meio afastamento
										
				/* calcular teta (superfície aproximada com a fórmula de jager) */
									
				Fd=(t0+a1*m)*(t0+a1*m)+a2*m*m;
				teta=(Fd+b2*h*h);
				teta=sqrt(teta);
				t[im][ih]=teta;				
				
			}
		
	}

		
}



void germam_t(float t0, float m0, float h0, float x0, float v0, float RN, float RNIP, float BETA, int nh, float dh, int nm, float dm, float **t) {
/*< Semblance da aproximação do CRS Quarta ordem parabólico (HöCHT, 2002) >*/	

		int im, ih; // contadores de laço
		float m; //coordenadas do CMP
		float h; // coordenadas do meio-afastamento
		float teta; //amostra no tempo

		/* Parâmetros da aproximação de germam_t */
		
		float s1=(2*sin(BETA)/v0);
		float s2=((cos(BETA)*cos(BETA))/(v0*RN));
		float s3=(cos(BETA)*cos(BETA)/(v0*RNIP));	
		float s4=((sin(BETA)*cos(BETA)*cos(BETA))/(v0*RN*RN));
		float s5=(sin(BETA)*cos(BETA)*cos(BETA)/(v0*RNIP*RNIP*RN))*(2*RNIP+RN);
		float s6=(cos(BETA)*cos(BETA)/(2*v0*RNIP*RNIP*RNIP*RN*RN))*(RNIP*RNIP*(8*cos(BETA)*cos(BETA)-6)+RNIP*RN*(5*cos(BETA)*cos(BETA)-4)-2*RN*RN*sin(BETA)*sin(BETA));
		float s7=(cos(BETA)*cos(BETA)*(5*cos(BETA)*cos(BETA)-4))/(4*v0*RN*RN*RN);
		float s8=cos(BETA)*cos(BETA)*(4*RNIP*sin(BETA)*sin(BETA)-RN*cos(BETA)*cos(BETA))/(4*v0*RNIP*RNIP*RNIP*RN);
		float g1, g2, g3, g4, g5;		
			
		for (im=0; im < nm; im++){
			
			m=im*dm+x0; //coordenada do cmp

			m=m-m0; // distância em relação ao CMP central
				
			for(ih=0;ih<nh;ih++){
			
				h=ih*dh+h0; // coordenada do meio afastamento
										
				/* APROXIMAÇÃO PARABÓLICA SRC QUARTA ORDEM germam_t */
                
				g1=t0+m*s1+s2*m*m+h*h*s3-s4*m*m*m;
				g2=-h*h*m*s5;
				g3=-h*h*m*m*s6;
				g4=-m*m*m*m*s7;
				g5=h*h*h*h*s8;
				teta=g1+g2+g3+g4+g5;
				t[im][ih]=teta;
				
			}
		
	}
		
}


void germam_t2(float t0, float m0, float h0, float x0, float v0, float RN, float RNIP, float BETA, int nh, float dh, int nm, float dm, float **t){
/*< Semblance da aproximação do CRS quarta ordem hiperbólico (HöCHT, 2002) >*/	

		int im, ih; // contadores de laço
		float m; //coordenadas do CMP
		float h; // coordenadas do meio-afastamento
		float teta; //amostra no tempo

		/* Parâmetros da aproximação de germam_t2 */
		float s1=t0*t0;
		float s2=((4*t0*sin(BETA))/v0);
		float s3=(2*(v0*t0*cos(BETA)*cos(BETA)+2*RN*sin(BETA)*sin(BETA))/(v0*v0*RN));
		float s4=((2*t0*cos(BETA)*cos(BETA))/(v0*RNIP));
		float s5=((2*sin(BETA)*cos(BETA)*cos(BETA)*(2*RN-v0*t0))/(v0*v0*RN*RN));
		float s6=(2*sin(BETA)*cos(BETA)*cos(BETA)*(2*RNIP*RN-2*v0*t0*RNIP-v0*t0*RN))/(v0*v0*RNIP*RNIP*RN);
		float s7=(cos(BETA)*cos(BETA)*(RN*(10*cos(BETA)*cos(BETA)-8)+v0*t0*(4-5*cos(BETA)*cos(BETA))))/(2*v0*v0*RN*RN*RN);
		float s8=(cos(BETA)*cos(BETA)/(v0*v0*RNIP*RNIP*RNIP*RN*RN))*(v0*t0*RNIP*RNIP*(6-8*cos(BETA)*cos(BETA))+v0*t0*RNIP*RN*(4-5*cos(BETA)*cos(BETA))+2*v0*t0*RN*RN*sin(BETA)*sin(BETA)-4*RNIP*RN*RN*sin(BETA)*sin(BETA)+RNIP*RNIP*RN*(10*cos(BETA)-8));
		float s9=(cos(BETA)*cos(BETA)*(4*v0*t0*RNIP*sin(BETA)*sin(BETA)-v0*t0*RN*cos(BETA)*cos(BETA)+2*RNIP*RN*cos(BETA)*cos(BETA)))/(2*v0*v0*RNIP*RNIP*RNIP*RN);
		float g1, g2, g3, g4, g5, g6;
								
		for (im=0; im < nm; im++){
			
			m=im*dm+x0; //coordenada do cmp

			m=m-m0; // distância em relação ao CMP central
				
			for(ih=0;ih<nh;ih++){
			
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
				t[im][ih]=teta;
				
			}
		
	}
		
}


void germam_tshift(float t0, float m0, float h0, float x0, float v0, float RN, float RNIP, float BETA, int nh, float dh, int nm, float dm, float **t) {
/*< Semblance da aproximação do CRS quarta ordem deslocado (HöCHT, 2002) >*/

		int im, ih; // contadores de laço
		float m; //coordenadas do CMP
		float h; // coordenadas do meio-afastamento
		float teta; //amostra no tempo

		/* Parâmetros da aproximação de germam_tshift */	
		float s1=(2*RNIP/v0)*(2*RNIP/v0);
		float s2=(8*RNIP*sin(BETA)/(v0*v0));
		float s3=4*((RNIP*cos(BETA)*cos(BETA)+RN*sin(BETA)*sin(BETA))/(v0*v0*RN));
		float s4=(4*cos(BETA)*cos(BETA)/(v0*v0));
		float s5=((4*sin(BETA)*cos(BETA)*cos(BETA)*(RN-RNIP))/(v0*v0*RN*RN));
		float s6=((8*sin(BETA)*cos(BETA)*cos(BETA))/(v0*v0*RN));
		float s7=(cos(BETA)*cos(BETA)/(v0*v0*RN*RN*RN))*(RN*(5*cos(BETA)*cos(BETA)-4)+RNIP*(4-5*cos(BETA)*cos(BETA)));
		float s8=((4*cos(BETA)*cos(BETA)*(3-4*cos(BETA)*cos(BETA)))/(v0*v0*RN*RN));
		float s9=(4*cos(BETA)*cos(BETA)*sin(BETA)*sin(BETA))/(v0*v0*RNIP*RN);
		float s10=t0-(2*RNIP/v0);
		float g1, g2, g3, g4, g5;		
			
		for (im=0; im < nm; im++){
			
			m=im*dm+x0; //coordenada do cmp

			m=m-m0; // distância em relação ao CMP central
				
			for(ih=0;ih<nh;ih++){
			
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
				t[im][ih]=teta;						

			}
		
	}
		
}


void pade_th(float t0, float m0, float h0, float x0, float v0, float RN, float RNIP, float BETA, int nh, float dh, int nm, float dm, float **t) {
/*< Semblance da aproximação do CRS Padé parabólico expansão em h (NEVES, 2017) >*/

		int im, ih; // contadores de laço
		float m; //coordenadas do CMP
		float h; // coordenadas do meio-afastamento
		float teta; //amostra no tempo	

		/* Parâmetros do CRS padé-t expansão em h */
		float s1=(2*sin(BETA)/v0);
		float s2=((cos(BETA)*cos(BETA))/(v0*RN));
		float s3=(sin(BETA)*cos(BETA)*cos(BETA))/(v0*RN*RN);
		float s4=(cos(BETA)*cos(BETA)*(5*cos(BETA)*cos(BETA)-4))/(4*v0*RN*RN*RN);
		float s5=(cos(BETA)*cos(BETA)/(v0*RNIP));
		float s6=(cos(BETA)*cos(BETA)/(2*v0*RNIP*RNIP*RNIP*RN*RN))*(RNIP*RNIP*(8*cos(BETA)*cos(BETA)-6)+RNIP*RN*(5*cos(BETA)*cos(BETA)-4)-2*RN*RN*sin(BETA)*sin(BETA));
		float s7=(sin(BETA)*cos(BETA)*cos(BETA)/(v0*RNIP*RNIP*RN))*(2*RNIP+RN);
		float s8=cos(BETA)*cos(BETA)*(4*RNIP*sin(BETA)*sin(BETA)-RN*cos(BETA)*cos(BETA))/(4*v0*RNIP*RNIP*RNIP*RN);
		float CO,C1,C2;		
			
		for (im=0; im < nm; im++){
			
			m=im*dm+x0; //coordenada do cmp

			m=m-m0; // distância em relação ao CMP central
				
			for(ih=0;ih<nh;ih++){
			
				h=ih*dh+h0; // coordenada do meio afastamento
										
				/* APROXIMAÇÃO PARABÓLICA SRC PADÉ EXPANSÃO EM h */
				CO=t0+m*s1+s2*m*m-s3*m*m*m-m*m*m*m*s4;
				C1=s5-m*m*s6-m*s7;
				C2=s8;
				
				
				teta=(h*h*C1)/(1+h*h*(-C2/C1));
				teta=CO + teta;
				
				t[im][ih]=teta; 
				
			}
		
	}
		
}


void pade_tm(float t0, float m0, float h0, float x0, float v0, float RN, float RNIP, float BETA, int nh, float dh, int nm, float dm, float **t){
/*< Semblance da aproximação do CRS Padé parabólico expansão em m (NEVES, 2017) >*/

		int im, ih; // contadores de laço
		float m; //coordenadas do CMP
		float h; // coordenadas do meio-afastamento
		float teta; //amostra no tempo	
	
		/* Parâmetros do CRS padé-t em m */
		float s1=(cos(BETA)*cos(BETA)/(v0*RNIP));
		float s2=((cos(BETA)*cos(BETA)*(4*RNIP*sin(BETA)*sin(BETA)-RN*cos(BETA)*cos(BETA)))/(4*v0*RNIP*RNIP*RNIP*RN));
		float s3=(2*sin(BETA)/v0);
		float s4=((sin(BETA)*cos(BETA)*cos(BETA)*(2*RNIP+RN))/(v0*RNIP*RNIP*RN));
		float kapa=RNIP*RNIP*(8*cos(BETA)*cos(BETA)-6)+RNIP*RN*(5*cos(BETA)*cos(BETA)-4)-2*RN*RN*sin(BETA)*sin(BETA);
		float s5=(cos(BETA)*cos(BETA)/(v0*RN));
		float s6=((cos(BETA)*cos(BETA)/(2*v0*RNIP*RNIP*RNIP*RN*RN))*kapa);
		float s7=-((sin(BETA)*cos(BETA)*cos(BETA))/(v0*RN*RN));
		float s8=-((cos(BETA)*cos(BETA)*(5*cos(BETA)*cos(BETA)-4))/(4*v0*RN*RN*RN));
		float CO,C1,C2,C3,C4,q1,q2,po,p1,p2,pp,qq;			
			
		for (im=0; im < nm; im++){
			
			m=im*dm+x0; //coordenada do cmp

			m=m-m0; // distância em relação ao CMP central
				
			for(ih=0;ih<nh;ih++){
			
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
				t[im][ih]=teta;
				
			}
		
	}
		
}


void pade_t2h(float t0, float m0, float h0, float x0, float v0, float RN, float RNIP, float BETA, int nh, float dh, int nm, float dm, float **t) {
/*< Semblance da aproximação do CRS Padé hiperbólico expansão em h (NEVES, 2017) >*/

		int im, ih; // contadores de laço
		float m; //coordenadas do CMP
		float h; // coordenadas do meio-afastamento
		float teta; //amostra no tempo	
	
		/* Parâmetros do CRS padé-t2 em h */
		double s1=t0*t0;
		double s2=((4*t0*sin(BETA))/v0);
		double s3=(2*(v0*t0*cos(BETA)*cos(BETA)+2*RN*sin(BETA)*sin(BETA))/(v0*v0*RN));
		double s4=((2*sin(BETA)*cos(BETA)*cos(BETA)*(2*RN-v0*t0))/(v0*v0*RN*RN));
		double s5=(cos(BETA)*cos(BETA)*(RN*(10*cos(BETA)*cos(BETA)-8)+v0*t0*(4-5*cos(BETA)*cos(BETA))))/(2*v0*v0*RN*RN*RN);
		double s6=((2*t0*cos(BETA)*cos(BETA))/(v0*RNIP));
		double s7=(2*sin(BETA)*cos(BETA)*cos(BETA)*(2*RNIP*RN-2*v0*t0*RNIP-v0*t0*RN))/(v0*v0*RNIP*RNIP*RN);
		double s8=(cos(BETA)*cos(BETA)/(v0*v0*RNIP*RNIP*RNIP*RN*RN))*(v0*t0*RNIP*RNIP*(6-8*cos(BETA)*cos(BETA))+v0*t0*RNIP*RN*(4-5*cos(BETA)*cos(BETA))+(2*v0*t0*RN*RN*sin(BETA)*sin(BETA)-4*RNIP*RN*RN*sin(BETA)*sin(BETA)+RNIP*RNIP*RN*(10*cos(BETA)*cos(BETA)-8)));
		double s9=((cos(BETA)*cos(BETA))/(2*v0*v0*RNIP*RNIP*RNIP*RN))*(4*v0*t0*RNIP*sin(BETA)*sin(BETA)-v0*t0*RN*cos(BETA)*cos(BETA)+2*RNIP*RN*cos(BETA)*cos(BETA));
		double CO,C1,C3;
			
		for (im=0; im < nm; im++){
			
			m=im*dm+x0; //coordenada do cmp

			m=m-m0; // distância em relação ao CMP central
				
			for(ih=0;ih<nh;ih++){
			
				h=ih*dh+h0; // coordenada do meio afastamento
										
				/* APROXIMAÇÃO QUADRÁTICA SRC-PADÉ EXPANSÃO EM h */
				CO=s1+s2*m+s3*m*m+s4*m*m*m+s5*m*m*m*m;
				C1=s6+s7*m+s8*m*m;
				C3=-s9/C1;
									
				teta=(h*h*C1)/(1+h*h*(C3));
				teta=CO+teta;
				teta=sqrt(teta); 
				t[im][ih]=teta;

				/*teta=CO+C1*h*h+(C3*h*h*h*h)/(1-C3*h*h/C1);
				teta=sqrt(teta); 
				t[im][ih]=teta;*/
				
			}
		
	}
		
}


void pade_t2m(float t0, float m0, float h0, float x0, float v0, float RN, float RNIP, float BETA, int nh, float dh, int nm, float dm, float **t) {
/*< Semblance da aproximação do CRS Padé hiperbólico expansão em m (NEVES, 2017) >*/

		int im, ih; // contadores de laço
		float m; //coordenadas do CMP
		float h; // coordenadas do meio-afastamento
		float teta; //amostra no tempo	
	
		/* Parâmetros do CRS padé-t2 em h */
		float s1=t0*t0;
		float s2=((2*t0*cos(BETA)*cos(BETA))/(v0*RNIP));
		float s3=(cos(BETA)*cos(BETA)*(4*v0*t0*RNIP*sin(BETA)*sin(BETA)-v0*t0*RN*cos(BETA)*cos(BETA)+2*RNIP*RN*cos(BETA)*cos(BETA)))/(2*v0*v0*RNIP*RNIP*RNIP*RN);
		float s4=(((4*t0*sin(BETA)))/v0)+(2*sin(BETA)*cos(BETA)*cos(BETA)*(2*RNIP*RN-2*v0*t0*RNIP-v0*t0*RN))/(v0*v0*RNIP*RNIP*RN);
		float s5=((2*(v0*t0*cos(BETA)*cos(BETA)+2*RN*sin(BETA)*sin(BETA)))/(v0*v0*RN));
		float s6=(cos(BETA)*cos(BETA)/(v0*v0*RNIP*RNIP*RNIP*RN*RN))*(v0*t0*RNIP*RNIP*(6-8*cos(BETA)*cos(BETA))+v0*t0*RNIP*RN*(4-5*cos(BETA)*cos(BETA))+2*v0*t0*RN*RN*sin(BETA)*sin(BETA)-4*RNIP*RN*RN*sin(BETA)*sin(BETA)+RNIP*RNIP*RN*(10*cos(BETA)-8));
		float s7=((2*sin(BETA)*cos(BETA)*cos(BETA)*(2*RN-v0*t0))/(v0*v0*RN*RN));
		float s8=(cos(BETA)*cos(BETA)*(RN*(10*cos(BETA)*cos(BETA)-8)+v0*t0*(4-5*cos(BETA)*cos(BETA))))/(2*v0*v0*RN*RN*RN);
		float CO,C1,C2,C3,C4,q1,q2,po,p1,p2,pp,qq;		

		for (im=0; im < nm; im++){
			
			m=im*dm+x0; //coordenada do cmp

			m=m-m0; // distância em relação ao CMP central
				
			for(ih=0;ih<nh;ih++){
			
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
				t[im][ih]=teta;
				
			}
		
	}
		
}


void pade_tsh(float t0, float m0, float h0, float x0, float v0, float RN, float RNIP, float BETA, int nh, float dh, int nm, float dm, float **t) {
/*< Semblance da aproximação do CRS Padé deslocada expansão em h (NEVES, 2017) >*/

		int im, ih; // contadores de laço
		float m; //coordenadas do CMP
		float h; // coordenadas do meio-afastamento
		float teta; //amostra no tempo	
	
		/* Parâmetros do CRS padé-tshift em h */
		float s1=(2*RNIP/v0)*(2*RNIP/v0);	
		float s2=((8*RNIP*sin(BETA))/(v0*v0));
		float s3=((RNIP*cos(BETA)*cos(BETA)+RN*sin(BETA)*sin(BETA))/(v0*v0*RN));
		float s4=((4*sin(BETA)*cos(BETA)*cos(BETA)*(RN-RNIP))/(v0*v0*RN*RN));
		float s5=(cos(BETA)*cos(BETA)/(v0*v0*RN*RN*RN))*(RN*(5*cos(BETA)*cos(BETA)-4)+RNIP*(4-5*cos(BETA)*cos(BETA)));
		float s6=((4*cos(BETA)*cos(BETA))/(v0*v0));
		float s7=((8*sin(BETA)*cos(BETA)*cos(BETA))/(v0*v0*RN));
		float s8=((4*cos(BETA)*cos(BETA)*(3-4*cos(BETA)*cos(BETA)))/(v0*v0*RN*RN));
		float s9=(4*cos(BETA)*cos(BETA)*sin(BETA)*sin(BETA))/(v0*v0*RNIP*RN);
		float s10=t0-(2*RNIP/v0);
		float CO,C1,C2;
			
		for (im=0; im < nm; im++){
			
			m=im*dm+x0; //coordenada do cmp

			m=m-m0; // distância em relação ao CMP central
				
			for(ih=0;ih<nh;ih++){
			
				h=ih*dh+h0; // coordenada do meio afastamento
										
			/* APROXIMAÇÃO DESLOCADA SRC-PADÉ EXPANSÃO EM h */
				CO=s1+m*s2+4*m*m*s3+m*m*m*s4+s5*m*m*m*m;
				C1=s6-m*s7+m*m*s8;
				C2=s9;
				
				
				teta=(h*h*C1)/(1+h*h*(-C2/C1));
				teta=CO + teta;
				teta=sqrt(teta);
				teta=teta+s10;
				t[im][ih]=teta;
				
			}
		
	}
		
}


void pade_tsm(float t0, float m0, float h0, float x0, float v0, float RN, float RNIP, float BETA, int nh, float dh, int nm, float dm, float **t) {
/*< Semblance da aproximação do CRS Padé deslocada expansão em m (NEVES, 2017) >*/

		int im, ih; // contadores de laço
		float m; //coordenadas do CMP
		float h; // coordenadas do meio-afastamento
		float teta; //amostra no tempo	
	
		/* Parâmetros do CRS padé-tshift em m */
		float s1=(2*RNIP/v0)*(2*RNIP/v0);
		float s2=((4*cos(BETA)*cos(BETA))/(v0*v0));
		float s3=(4*cos(BETA)*cos(BETA)*sin(BETA)*sin(BETA))/(v0*v0*RNIP*RN);
		float s4=((8*RNIP*sin(BETA))/(v0*v0));
		float s5=((8*sin(BETA)*cos(BETA)*cos(BETA))/(v0*v0*RN));
		float s6=4*((RNIP*cos(BETA)*cos(BETA)+RN*sin(BETA)*sin(BETA))/(v0*v0*RN));
		float s7=((4*cos(BETA)*cos(BETA)*(3-4*cos(BETA)*cos(BETA)))/(v0*v0*RN*RN));
		float s8=((4*sin(BETA)*cos(BETA)*cos(BETA)*(RN-RNIP))/(v0*v0*RN*RN));
		float s9=(cos(BETA)*cos(BETA))*(RN*(5*cos(BETA)*cos(BETA)-4)+RNIP*(4-5*cos(BETA)*cos(BETA)));
		float s10=(v0*v0*RN*RN*RN);
		float s11=t0-(2*RNIP/v0);
		float CO,C1,C2,C3,C4,q1,q2,po,p1,p2,pp,qq;			
			
		for (im=0; im < nm; im++){
			
			m=im*dm+x0; //coordenada do cmp

			m=m-m0; // distância em relação ao CMP central
				
			for(ih=0;ih<nh;ih++){
			
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
				t[im][ih]=teta;

			}
		
	}
	
}
