/* This file is automatically generated. DO NOT EDIT! */

#ifndef _f_crs_h
#define _f_crs_h


#define PI 3.14159
#include <stdio.h> // biblioteca padrão define operações de entrada/saída
#include <stdlib.h> // biblioteca padrão define alocação de memória, controle de processos
#include <math.h> // biblioteca de funções matemáticas como exp()
#include <stdlib.h> // necessário para gerar número aleatório com: rand() e srand()
#include <time.h> // necessário para gerar semente de numero aleatório: time()
#include <string.h> // necessário para lidar com strings e strcat()
#include <rsf.h> // biblioteca padrão do madagascar


void f_vfsa_aviso(int app);
/*< Informar usuário sobre a aproximação de tempo de trânsito CRS escolhida >*/


float sinal(float s);
/*< função sinal >*/


void fomel(float t0, float m0, float h0, float x0, float v0, float RN, float RNIP, float BETA, int nh, float dh, int nm, float dm, float **t);
/*< Aproximação de tempo de trânsito do CRS não hiperbólico (FOMEL; KAZINNIK, 2013) >*/


void jager(float t0, float m0, float h0, float x0, float v0, float RN, float RNIP, float BETA, int nh, float dh, int nm, float dm, float **t);
/*< Semblance da aproximação do CRS hiperbólico (JAGER et al., 2001) >*/


void germam_t(float t0, float m0, float h0, float x0, float v0, float RN, float RNIP, float BETA, int nh, float dh, int nm, float dm, float **t);
/*< Semblance da aproximação do CRS Quarta ordem parabólico (HöCHT, 2002) >*/


void germam_t2(float t0, float m0, float h0, float x0, float v0, float RN, float RNIP, float BETA, int nh, float dh, int nm, float dm, float **t);
/*< Semblance da aproximação do CRS quarta ordem hiperbólico (HöCHT, 2002) >*/


void germam_tshift(float t0, float m0, float h0, float x0, float v0, float RN, float RNIP, float BETA, int nh, float dh, int nm, float dm, float **t);
/*< Semblance da aproximação do CRS quarta ordem deslocado (HöCHT, 2002) >*/


void pade_th(float t0, float m0, float h0, float x0, float v0, float RN, float RNIP, float BETA, int nh, float dh, int nm, float dm, float **t);
/*< Semblance da aproximação do CRS Padé parabólico expansão em h (NEVES, 2017) >*/


void pade_tm(float t0, float m0, float h0, float x0, float v0, float RN, float RNIP, float BETA, int nh, float dh, int nm, float dm, float **t);
/*< Semblance da aproximação do CRS Padé parabólico expansão em m (NEVES, 2017) >*/


void pade_t2h(float t0, float m0, float h0, float x0, float v0, float RN, float RNIP, float BETA, int nh, float dh, int nm, float dm, float **t);
/*< Semblance da aproximação do CRS Padé hiperbólico expansão em h (NEVES, 2017) >*/


void pade_t2m(float t0, float m0, float h0, float x0, float v0, float RN, float RNIP, float BETA, int nh, float dh, int nm, float dm, float **t);
/*< Semblance da aproximação do CRS Padé hiperbólico expansão em m (NEVES, 2017) >*/


void pade_tsh(float t0, float m0, float h0, float x0, float v0, float RN, float RNIP, float BETA, int nh, float dh, int nm, float dm, float **t);
/*< Semblance da aproximação do CRS Padé deslocada expansão em h (NEVES, 2017) >*/


void pade_tsm(float t0, float m0, float h0, float x0, float v0, float RN, float RNIP, float BETA, int nh, float dh, int nm, float dm, float **t);
/*< Semblance da aproximação do CRS Padé deslocada expansão em m (NEVES, 2017) >*/

#endif
