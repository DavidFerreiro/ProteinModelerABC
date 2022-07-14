/* Program DeltaGREM
   Computes contact free energy (native-unfolding-misfolding)
   of a set of sequences with respect to the best of a list of
   target structures. Sequences can also be given as mutations.
*/

float LMIN=0.5; // Minimum length for selecting alignments
int PRINT_MUT=1; // Print details of mutations?
float SEC_STR=0;
char FILE_STR[200];

#define FILE_CODE_DEF "gen_code_ATGC.in"

#include "REM.h"
#include "coord.h"
#include "mut_del.h"
#include "alignments.h"
#include "gen_code.h"
#include "allocate.h"
#include "protein3.h"
#include "read_pdb.h"
#include "random3.h"           /* Generating random numbers */
#include "mutation.h"
#include "codes.h"
#include "input.h"
#include "get_para_DeltaGREM.h"
#include "Sec_str_all.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>

// External variables
int DTAIL_MAX=20;
float cont_thr=4.5;
int LEN; // For REM calculations
char AA_code[21];
float hydro[21];
int Verbose=0;

int ini_print;


/* 
   Give as input a PDB file and either an alignment file or a file with
   mutations. Output is the DeltaG of wild type and mutant sequences.
*/

#define N_CHAR 300          // Max. length of file names

struct result{
  float DeltaG_1;   // Mean of DeltaG over native models
  float DeltaG_2;   // Meean square of DeltaG over native models
  float DeltaG_min; // Minimum of DeltaG over native models
  float G_nat;      // Minimum of Gnat over native models
  int pdbbest;      // Best model where minimum is attained
  float seqid_opt;  // Largest seqid over possible models
  float seqid_mod;  // Seqid for best energy model
  float hydro;      // Average hydrophobicity
  int n, L;         // Number of models and length of optimal model
  int nmut;         // Number of mutations
};

// Input codes
int Remove_gaps(char **MSA, int n_seq, int L_ali, char **name,
		char *node, int c, char gap);

// Output
struct result *Allocate_res(int n);
void Record(struct result *r, struct REM E, int i_pdb, float seq_id, int L);
char *Get_name(char *file_name);
void Print_result(struct result *r, double *DDG_1, double *DDG_2, int *n,
		  char *name, char **name_pdb, float DG_wt, float G_nat_wt,
		  float hydro_wt, float fit_wt, FILE *file_out);
void Print_mut(short *aa_seq, short *aa_seq0, int **C_nat, int len_amm,
	       FILE *file_mut_out);
void Print_DG(struct REM E, int *i_sec, char *sec_str,
	      char *name_pdb, int ini, char *name);
void Empty_prot(struct protein *prot);

// Input parameters
// A: Input files defining the protein
static char dir_out[N_CHAR];
// B: Thermodynamic parameters
float sC1=0.065, sC0=0, sU1=0.140; // Configuration entropy
int REM=2;   // Use 1st (1), 2nd (2) and 3rd (3) moment of misfolded energy 
// C2: Mean-field model
int MEANFIELD=1;  // Compute site-specific mean-field distributions
int OPT_LAMBDA=1; // 1=Optimize Lambda by maximum likelihood
float LAMBDA=0.9; // Lambda parameter for meanfield if OPT_LAMBDA=0
char MODEL[N_CHAR]="ALL"; // Type of mean-field model (ALL, NAT, DG)
float DG_OPT=-1;  // DG target of the optimization if MODEL="DG"
// D1: Mutation model, P_mut
#define MUTPAR 10
float mut_par[MUTPAR], tt_ratio=1, kCpG=1;
// D1: Mutation model, exchangeability
char MATRIX[40]="WAG"; // Empirical exchangeability matrix
char EXCHANGE='F'; // exchangeability model. M=MUT F=FLUX Q=RATE E=EXCH
float TWONUC=0; // Rate of two to one nucleotide substitutions  if MUT
// E: Output
int FORMAT=1;   // PAML format

// Derived data
unsigned long iran;
static int len_amm;
float DG_thr;

float Econt[21][21];
float **Econt_T=NULL;
float A_LOC=1; // Coefficient for local interactions

int main(int argc, char **argv){


  /***********************
          INPUT
  ************************/
  // Input files
  //char Input_dir[N_CHAR];
  char **file_pdb, chain[]="\0\0\0";
  char *file_ali=NULL, *file_mut=NULL;

  /***********************
         SEQUENCES
  ************************/
  //int nuc_mut, res_mut=0, aa_new;

  /***********************
         Wild Type
  ************************/
  //float fitness_wt;

  /***********************
         MUTANT
  ************************/
  //float fitness_mut;

  /***********************
          DUMMIES
  ************************/
  int i, j;

  /***********************
          OUTPUT
  ************************/

  /******************** Input operations   ************************/
  float TEMP=0.5; // default
  int N_pdb=Get_para(argc, argv, &file_pdb, chain,
		     FILE_STR, &TEMP, &sU1, &sC0, &sC1, &REM, &SEC_STR,
		     &file_ali, &file_mut, dir_out);
  /*printf("%d PDB files to read\n", N_pdb);*/
  /*printf("Temperature= %.3f\n", TEMP);*/

  // alignment
  int n_seq=0, L_ali=0, *selected=NULL;
  char **name_seq;
  char **MSA=Read_MSA(&n_seq, &L_ali, &name_seq, &selected, file_ali, LMIN);
  if(MSA==NULL){printf("No multiple alignment given, checking mutations\n");}
  else{Remove_gaps(MSA, n_seq, L_ali, name_seq, "node", 4, 'L');}
 
  // Mutations-deletions
  int Nmut=0; struct mutation *mutations=NULL;
  int Ndel=0; struct deletion *deletions=NULL;
  if(file_mut){
    Nmut=Read_mut(&mutations, file_mut, &deletions, &Ndel);
    printf("%d mutations, %d deletions\n", Nmut, Ndel);
  }
  if( n_seq &&((Nmut)||(Ndel))){
      printf("ERROR, it is not allowed to input MSA (%d seq.) and ", n_seq);
      printf("mutations (%d a.a. changes and %d indels)\n", Nmut, Ndel);
      exit(8);
  }

  /**************************  Allocate  ****************************/
  char **name_pdb=malloc(N_pdb*sizeof(char *));
  float *seq_id=NULL; int *seq_L=NULL, i_pdb;
  if(n_seq){
    seq_id=malloc(n_seq*sizeof(float));
    seq_L=malloc(n_seq*sizeof(int));
  }
  struct result *result_seq=Allocate_res(n_seq);
  struct result *result_mut=Allocate_res(Nmut);
  struct result *result_del=Allocate_res(Ndel);
  struct result *result_wt=Allocate_res(1);

  int n_pdb=0, lchain=1;
  if(chain[0]=='\0')lchain=0; // chain is not defined
  for(i_pdb=0; i_pdb<N_pdb; i_pdb++){
    
    int N_mod=Count_models_PDB(file_pdb[i_pdb]), imod;
    if(N_mod<0){
      printf("WARNING, no PDB file found: %s\n", file_pdb[i_pdb]);
      continue;
    }
    /*printf("PDB file %s, %d models\n", file_pdb[i_pdb], N_mod);*/
    name_pdb[i_pdb]=Get_name(file_pdb[i_pdb]);
    if(lchain==0)chain[0]='\0'; // Read the first chain

    for(imod=0; imod<N_mod; imod++){
      if(imod==0)n_pdb++;

      struct residue *res=NULL; struct res_short *res_short; int *res_index;
      len_amm=Get_pdb(&target, &res, &res_short, file_pdb[i_pdb],
		      chain, &res_index);
      //int kmod; if(N_mod==1){kmod=-1;}else{kmod=imod;}
      //len_amm=Get_pdb(&target, &res, file_pdb[i_pdb], chain, kmod);
      if(len_amm<=0){
	printf("WARNING, no chain found: %s chain %s\n",
	       file_pdb[i_pdb], chain); continue;
      }

      int *i_sec=target.i_sec;
      short *aa_seq0=target.aa_seq;
      int ali_seq[len_amm]; 
      if(MSA){ // Align PDB with MSA
	char seq_PDB[len_amm];
	for(i=0; i<len_amm; i++)seq_PDB[i]=Amin_code(aa_seq0[i]);
	int k_pdb=Find_seq(ali_seq,seq_L,seq_id, seq_PDB,len_amm,
			   MSA,n_seq,L_ali);
	if(k_pdb<0){
	  printf("WARNING, PDB sequence %s chain %s not found in MSA %s\n",
		 file_pdb[i_pdb], chain, file_ali);
	  free(res); Empty_prot(&target); continue;
	}
      }

      /******************** Folding stability **********************/
      int **C_nat=Fill_C_nat(len_amm, target.contact);
      float S_C=sC0+len_amm*sC1, S_U=len_amm*sU1;
      struct REM E_wt;
      Initialize_E_REM(&E_wt, len_amm, REM, TEMP, S_C, S_U, FILE_STR, 0);
      short aa_seq[len_amm];

      if(MSA){
	for(i=0; i<n_seq; i++){
	  for(j=0; j<len_amm; j++){
	    if(ali_seq[j]>=0){aa_seq[j]=Code_AA(MSA[i][ali_seq[j]]);}
	    else{aa_seq[j]=20;} // gap
	  }	  
	  //E_wt.DeltaG=Compute_DG_overT_contfreq(&E_wt, aa_seq, C_nat, i_sec);
	  E_wt.DeltaG=Compute_DG_overT_threading(&E_wt, aa_seq, C_nat, i_sec);
	  E_wt.DeltaG+=(sC1+sU1)*(seq_L[i]-len_amm);
	  if(imod==0)result_seq[i].hydro=Compute_hydro(aa_seq, len_amm);
	  Record(result_seq+i, E_wt, i_pdb, seq_id[i], seq_L[i]);
	}

      }else{

	// Wild type
	//E_wt.DeltaG=Compute_DG_overT_contfreq(&E_wt, aa_seq0, C_nat, i_sec);
	E_wt.DeltaG=Compute_DG_overT_threading(&E_wt, aa_seq0, C_nat, i_sec);
	if(imod==0)result_wt[0].hydro=Compute_hydro(aa_seq0, len_amm);
	Record(result_wt, E_wt, i_pdb, 1.000, len_amm);
	Print_DG(E_wt, i_sec, target.sec_str, name_pdb[i_pdb], 1, "WT");
	Test_contfreq(&E_wt, aa_seq0, C_nat, i_sec, name_pdb[i_pdb]);
	//Test_contfreq(&E_wt, aa_seq0, name_pdb[i_pdb]);

	struct REM E_mut;
	Initialize_E_REM(&E_mut, len_amm, REM, TEMP, S_C, S_U, FILE_STR, 0);

	// Mutation
	FILE *file_mut_out=NULL;
	if(PRINT_MUT && Nmut){
	  char name_mut[900];
	  sprintf(name_mut, "%s_mutations.dat", file_mut);
	  file_mut_out=fopen(name_mut, "w");
	  //printf("Writing %s\n", name_mut);
	}
	for(i=0; i<Nmut; i++){
	  for(j=0; j<len_amm; j++)aa_seq[j]=aa_seq0[j];
	  if(Construct_mut(aa_seq, mutations+i, res, len_amm)<0)continue;
	  Print_mut(aa_seq, aa_seq0, C_nat, len_amm, file_mut_out);
	  printf("Mutation %s\n", mutations[i].name);
	  E_mut.DeltaG=
	    //Compute_DG_overT_contfreq(&E_mut, aa_seq, C_nat, i_sec);
	    Compute_DG_overT_threading(&E_mut, aa_seq, C_nat, i_sec);
	  if(imod==0)result_mut[i].hydro=Compute_hydro(aa_seq, len_amm);
	  // Statistics
	  float id=1.-mutations[i].nmut/(float)len_amm;
	  Record(result_mut+i, E_mut, i_pdb, id, len_amm);
	  result_mut[i].nmut=mutations[i].nmut;
	  Print_DG(E_mut, NULL, NULL, name_pdb[i_pdb], 0, mutations[i].name);
	}
      
	// Deletion
	for(i=0; i<Ndel; i++){
	  for(j=0; j<len_amm; j++)aa_seq[j]=aa_seq0[j];
	  if(Construct_del(aa_seq, deletions+i, res, len_amm)<0)continue;
	  int del=deletions[i].Ldel;
	  printf("Deletion %s %d\n", deletions[i].name, del);
	  E_mut.DeltaG=
	    //Compute_DG_overT_contfreq(&E_mut, aa_seq, C_nat, i_sec);
	    Compute_DG_overT_threading(&E_mut, aa_seq, C_nat, i_sec);
	  E_mut.DeltaG-=(sC1+sU1)*del;
	  if(imod==0)result_del[i].hydro=Compute_hydro(aa_seq, len_amm);
	  // Statistics
	  float id=1.-del/(float)len_amm;
	  Record(result_del+i, E_mut, i_pdb, id, len_amm-del);
	}
      }
      
      for(i=0; i<len_amm; i++){free(C_nat[i]);} free(C_nat);
      free(res);
      Empty_prot(&target);
    }
  }
  if(n_pdb==0){
    printf("ERROR, no pdb file found\n"); exit(8);
  }
  
  /************************* Output files ***************************/
  char nameout[300], para[80];
  sprintf(para, "T%.2f_SU%.3f_SC%.3f", TEMP, sU1, sC1);
  if(n_seq>0){
    // MSA
    sprintf(nameout, "%s%s", dir_out, file_ali);
  }else if((Nmut) || (Ndel)){
    sprintf(nameout, "%s%s", dir_out, file_mut);
  }else{
    sprintf(nameout, "%s%s", dir_out, name_pdb[0]);
  }
  strcat(nameout, "_DeltaG.dat");

  char header[2000];
  sprintf(header,
	  "# T= %.2f sU1= %.3f sC1= %.3f sC0= %.3f REM= %d\n"
	  "# %d PDB files used\n", TEMP, sU1, sC1, sC0, REM, n_pdb);
  if((Nmut==0)&&(Ndel==0)){
    strcat(header, "#L\tDelta_G/L(best)\tsigma(Delta_G)/L\tG_nat/N");
    if(n_pdb>1)strcat(header, "\tseq.id(best)");
    strcat(header, "\tseq.id(model)\thydro\tlog(f)\ttemplate\tsequence");
  }else{
    strcat(header, "#L\tDDG\tsigma(DDG)\tDG_nat");
    if(n_pdb>1)strcat(header,"\tseq.id(best)");
    strcat(header,"\tseq.id(model)\thydro\tlog(f)\tDeltaG/nmut\ttemplate\tmut");
  }

  FILE *file_out=fopen(nameout, "w");
  //Mi primer codigo C :'(
  char *numero;
  char *checker = NULL;
  int val;
  int len;
  numero = nameout;
  len = strlen(numero);
  checker = strstr(numero, "sequences");
  if (checker == numero) {
      numero[len-10] = '\0';
      numero += 9;
      val = atoi(numero);
      /*if (val == 1){
          printf(" Calculating simulation # %d summary statistics\n", val);
          }*/
      if (val % 50 == 0){
          printf(" Calculating simulation # %d summary statistics\n", val);
          }
  }
  
  //printf(" Writing %s\n", nameout);
  fprintf(file_out, "%s\n", header);
  double DDG_1=0, DDG_2=0; int n=0;
  double Dh_1=0, Dh_2=0;
  float DG_wt=result_wt[0].DeltaG_min;
  float G_nat_wt=result_wt[0].G_nat;
  float hydro_wt=result_wt[0].hydro, h;
  float fit_wt=-log(1+exp(DG_wt));

  if(n_seq){
    for(i=0; i<n_seq; i++){
      Print_result(result_seq+i, &DDG_1, &DDG_2, &n, name_seq[i],
		   name_pdb, 0, 0, 0, 0, file_out);
    }
    h=result_seq[i].hydro; Dh_1+=h; Dh_2+=h*h;
  }else{
    result_wt[0].nmut=1;
    Print_result(result_wt, &DDG_1, &DDG_2, &n, "WT",
		 name_pdb, 0, 0, 0, 0, file_out);

    DDG_1=0; DDG_2=0; n=0; Dh_1=0; Dh_2=0;
    for(i=0; i<Nmut; i++){
      Print_result(result_mut+i, &DDG_1, &DDG_2, &n, mutations[i].name,
		   name_pdb, DG_wt, G_nat_wt, hydro_wt, fit_wt, file_out);
      h=result_mut[i].hydro; Dh_1+=h; Dh_2+=h*h;
    }
    for(i=0; i<Ndel; i++){
      result_del[i].nmut=1;
      Print_result(result_del+i, &DDG_1, &DDG_2, &n, deletions[i].name,
		   name_pdb, DG_wt, G_nat_wt, hydro_wt, fit_wt, file_out);
      h=result_del[i].hydro; Dh_1+=h; Dh_2+=h*h;
    }
  }

  if(n>1){
    DDG_1/=n;
    DDG_2=sqrt((DDG_2/n-DDG_1*DDG_1)/(n-1));
    DDG_1-=DG_wt;
    fprintf(file_out, "# Mean (D)DG= %.3f Error= %.3f\n", DDG_1, DDG_2);
    /*printf("%s Mean (D)DG= %.3f Error= %.3f\n", file_mut, DDG_1, DDG_2);*/
    Dh_1/=n;
    Dh_2=sqrt((Dh_2/n-Dh_1*Dh_1)/(n-1));
    Dh_1-=hydro_wt;
    fprintf(file_out, "# Mean (D)hydro= %.3g Error= %.2g\n", Dh_1, Dh_2);
  }
  fclose(file_out);

  return(0);
}
 
void Empty_prot(struct protein *prot){
  if(prot->length==0)return;
  for(int i=0; i<prot->length; i++){
    free(prot->contact[i]);
  }
  free(prot->contact);
  free(prot->cont_list);
  free(prot->aa_seq);
  free(prot->sec_str);
  free(prot->i_sec);
}


char *Get_name(char *file_name){
  char *name=malloc(100*sizeof(char)), *n=name;
  char *s=file_name, *s1=s;
  while(*s!='\0'){if(*s=='/')s1=s; s++;} s1++;
  while(*s1!='\0'){if(*s1=='.')break; *n=*s1; n++; s1++;}
  *n='\0';
  return(name);
}

int Remove_gaps(char **MSA, int n_seq, int L_ali, char **name,
		char *node, int c, char cgap){
  int iseq, i, Lgap=0, inode=0, ini=1, undecided=0, skip=0, mseq=0;
  for(i=0; i<L_ali; i++){
    int n_gap=0, n_c=0;
    for(iseq=0; iseq<n_seq; iseq++){
      if(strncmp(name[iseq], node, c)==0)continue;
      if(MSA[iseq][i]=='-'){n_gap++;}
      else if(MSA[iseq][i]==cgap){n_c++;}
      if(i==0)mseq++;
    }
    if(n_gap==0)continue;    // No gap present in real sequence
    if(n_c)undecided++;      // Not clear decision
    if(n_c > n_gap){skip++; continue;} // More characters than gaps in column
    for(iseq=0; iseq<n_seq; iseq++){
      if(strncmp(name[iseq], node, c)!=0)continue;
      if(MSA[iseq][i]==cgap){MSA[iseq][i]='-'; Lgap++;}
      if(ini)inode++;
    }
    if(ini)ini=0;
  }
  /*printf("%d undecided columns over %d, %d of them treated as gap\n",
	 undecided, L_ali, undecided-skip);
  printf("%d gap characters introduced in %d internal nodes\n",
	 Lgap, inode);*/
  //exit(8);

  char nameout[50]; strcpy(nameout, "newalignment.fasta");
  FILE *file_out=fopen(nameout, "w");
  for(iseq=0; iseq<n_seq; iseq++){
    fprintf(file_out, ">%s\n", name[iseq]);
    for(i=0; i<L_ali; i++)fprintf(file_out, "%c", MSA[iseq][i]);
    fprintf(file_out, "\n");
  }
  fclose(file_out);
  /*printf("Writing %s\n", nameout);*/
  return(0);
}

struct result *Allocate_res(int n)
{
  if(n==0)return(NULL);
  struct result *res=malloc(n*sizeof(struct result)), *r=res;
  int i;
  for(i=0; i<n; i++){
    r->DeltaG_1=0;
    r->DeltaG_2=0;
    r->DeltaG_min=0;
    r->G_nat=0;
    r->pdbbest=-1;
    r->seqid_opt=0;
    r->seqid_mod=0;
    r->n=0;
    r->nmut=0;
    r++;
  }
  return(res);
}

void Record(struct result *r, struct REM E, int i_pdb, float seq_id, int L)
{
  r->DeltaG_1+=E.DeltaG;
  r->DeltaG_2+=E.DeltaG*E.DeltaG;
  if((r->n==0)||(E.DeltaG<r->DeltaG_min)){
    r->DeltaG_min=E.DeltaG; r->G_nat=E.E_nat;
    r->pdbbest=i_pdb; r->seqid_mod=seq_id; r->L=L;
  }
  if((r->n==0)||(seq_id > r->seqid_opt))r->seqid_opt=seq_id;
  (r->n)++;
}

void Print_result(struct result *r, double *DDG_1, double *DDG_2, int *n,
		  char *name, char **name_pdb, float DG_wt, float G_nat_wt,
		  float hydro_wt, float fit_wt, FILE *file_out)
{
  if(r->n==0)return;
  *DDG_1+=r->DeltaG_min;
  *DDG_2+=r->DeltaG_min*r->DeltaG_min;
  (*n)++;
  r->DeltaG_1/=r->n;
  float s=r->DeltaG_2-r->n*r->DeltaG_1*r->DeltaG_1;
  if(r->n>1)s=sqrt(s/(r->n-1));
  float fit=-log(1+exp(r->DeltaG_min));


  fprintf(file_out, "%3d", r->L);
  if(DG_wt==0){ // wild type
    fprintf(file_out, "\t%.4f", r->DeltaG_min/r->L);
    fprintf(file_out, "\t%.2f", s/r->L);
    fprintf(file_out, "\t%.4f", r->G_nat/r->L);
    if(r->n>1)fprintf(file_out, "\t%.3f", r->seqid_opt);
    fprintf(file_out, "\t%.3f", r->seqid_mod);
    fprintf(file_out, "\t%.4f", r->hydro);
    fprintf(file_out, "\t%.3g", fit);
    fprintf(file_out, "\t%.3f", 0.00);
  }else{
    float DDG=r->DeltaG_min-DG_wt;
    fprintf(file_out, "\t%6.3f", DDG);
    fprintf(file_out, "\t%.2f", s);
    fprintf(file_out, "\t%6.3f", r->G_nat-G_nat_wt);
    if(r->n>1)fprintf(file_out, "\t%.3f", r->seqid_opt);
    fprintf(file_out, "\t%.3f", r->seqid_mod);
    fprintf(file_out, "\t%.4f", r->hydro);
    fprintf(file_out, "\t%.3g", fit);
    if(r->nmut)DDG/=r->nmut;
    fprintf(file_out, "\t%.3f", DDG);
  }
  fprintf(file_out, "\t%s", name_pdb[r->pdbbest]);
  fprintf(file_out, "\t%s", name);
  fprintf(file_out, "\n");
}

void Print_mut(short *aa_seq, short *aa_seq0, int **C_nat, int len_amm,
	       FILE *file_out)
{
  if(file_out==NULL)return;
  int i, j;
  for(i=0; i<len_amm; i++){
    if(aa_seq[i]==aa_seq0[i])continue;
    float *U=Econt[aa_seq[i]];
    float *U_old=Econt[aa_seq0[i]];
    for(j=0; j<len_amm; j++){
      if(C_nat[i][j]==0)continue;
      fprintf(file_out, "%c%c-%c%c-%d ",
	      Amin_code(aa_seq0[i]), Amin_code(aa_seq0[j]),
	      Amin_code(aa_seq[i]), Amin_code(aa_seq[j]), j);
      fprintf(file_out, " %.3f ", U[aa_seq[j]]-U_old[aa_seq0[j]]);
    }
  }
  fprintf(file_out, "\n");
}

void Print_DG(struct REM E, int *i_sec, char *sec_str,
	      char *name_pdb, int ini, char *name)
{
  char nameout[200];
  sprintf(nameout, "%s_DG.dat", name_pdb);
  char mode[4];
  if(ini){strcpy(mode, "w");}else{strcpy(mode, "a");}
  FILE *file_out=fopen(nameout, mode);
  if(ini){
    /*printf("Writing %s\n", nameout);*/
    if(SEC_STR){
      fprintf(file_out, "# sec_str: ");
      for(int i=0; i<E.L; i++)fprintf(file_out, "%c", sec_str[i]);
      fprintf(file_out, "\n");
      fprintf(file_out, "# sec_str: ");
      for(int i=0; i<E.L; i++)fprintf(file_out, "%c", SEC_EL[i_sec[i]]);
      fprintf(file_out, "\n");
    }
    fprintf(file_out, "#E_nat\tE_loc\tE1\tE2\tE3\tL\tDG\tmut\n");
  }
  fprintf(file_out, "%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%d\t%.2f\t%s\n",
	  E.E_nat, E.E_loc, E.E1, -E.E2, E.E3, E.L, E.DeltaG, name);
  fclose(file_out);
}
