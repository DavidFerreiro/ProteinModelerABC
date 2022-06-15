// Mutation parameters
//extern char coded_aa[64], *codon[64]; 
#define MUTPAR 10
extern float mut_par[MUTPAR];
extern int NPAR;
// contact definition
//extern char cont_type;
extern float cont_thr, cont_thr2;
// contact statistics
extern int DTAIL_MAX;
extern int LEN_MAX;
extern int LEN;
extern float *Cont_L;
extern float *Cont_1L;
extern float *Cont_2L;
extern float *Cont_3L;
extern float *C2_L;
extern float *C3_L;
extern float *Cnc1_L;
extern float **nc_nc_L;
extern float **nc_nc_Nc_L;
extern char FILE_STR[200];
// Contact energy function
extern char AA_code[21];
//extern float Econt[21][21];
extern float **Econt_T, T_ratio;
extern float **E_loc_over_T;
extern float hydro[21];
extern float SEC_STR;
//extern char SEC_EL[16];
extern float TEMP;
extern float sC1, sC0, SC, s0; //, Conf_entropy, K2Thr
// general
extern int Verbose;
// AA distr
extern float Entr_ave, Entr_reg;
extern int *i_sec;

