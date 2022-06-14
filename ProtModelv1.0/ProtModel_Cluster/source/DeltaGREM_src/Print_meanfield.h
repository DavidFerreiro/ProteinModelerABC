int Print_profiles(float **P_MF_ia, char *tag, double DG_ave, float Lambda,
		   float *P_mut_a, short *aa_seq, int L,
		   char *name_file, int PRINT, float *wi, FILE *out);
int Print_exchange(float **P_MF_ia, char *TAG,
		   struct res_short *res, int L,
		   float *P_mut_a, float *P_cod, float **Q_cod,
		   float tt_ratio, float TWONUC,
		   char *nameout, int FORMAT, char EXCHANGE,
		   char *MATRIX, int PRINT_EXCH_ALL, float *wi, FILE *out);
/*int Print_subst_rate(float **P_MF_ia, float *P_mut_a,
		     float *P_cod, float **Q_cod,
		     float tt_ratio, short *aa_seq,
		     int **Cnat, char *sec_str, int L,
		     char *name_file, float TWONUC);*/
int Print_profile_evo(char *name, double **P_ia, short *aa_seq, int L,
		      double DG_ave, long it_sum);
float Entropy(float *P, int n);
