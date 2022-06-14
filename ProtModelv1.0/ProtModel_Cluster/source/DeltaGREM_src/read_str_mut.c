#include "coord.h"
#include "protein3.h"
#include "read_str_mut.h"
#include "codes.h"
#include "allocate.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

static int Get_rank(int *code_AA, char *string);

float **Str_mut_matrix(char *STR_MUT_TYPE, char *file_str_mut,
		       struct protein target, int *res_index)
{
  if(file_str_mut[0]=='\0')return(NULL);
  FILE *file_in=fopen(file_str_mut, "r");
  if(file_in==NULL){
    printf("WARNING, file %s with structural mutations does not exist\n",
	   file_str_mut); return(NULL);
  }
  int wrong_thr=10;
  int L=target.L_PDB; short *aseq=target.aa_seq;
  char string[1000]; int k=0, wrong=0;
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#')continue;
    if(k>=L){
      printf("WARNING, line %s in file %s > number of residues %d\n",
	     string, file_str_mut, L); wrong++; k++; continue;
    }
    int i=res_index[k];
    if(i<0 || i>= target.length){
      printf("WARNING, sequence index of res %d = %d not in [0,%d]\n",
	     k, i, target.length-1); wrong++;
    }else if(string[0]!=AMIN_CODE[aseq[i]]){
      printf("WARNING, amino acid %c%d in file %s different from %c ",
	     string[0], k+1, file_str_mut, AMIN_CODE[aseq[i]]);
      printf("in target sequence\n"); wrong++;
    }
    k++;
  }
  fclose(file_in);

  if(k!=L){
    printf("WARNING, file %s does contains %d lines but %d residues in prot\n",
	   file_str_mut, k, L); wrong+=abs(k-L);
  }
  if(wrong){
    printf("WARNING, there were %d errors in file %s\n", wrong, file_str_mut);
    if(wrong>wrong_thr){printf("str mut not read\n"); return(NULL);}
  }

  printf("Reading predicted structural effects of mutations in file %s\n",
	 file_str_mut);
  file_in=fopen(file_str_mut, "r");
  int len_amm=target.length;
  float **Str_mut=Allocate_mat2_f(len_amm, 20, "Str_mut");
  int num[len_amm]; for(k=0; k<len_amm; k++)num[k]=0;
  k=0; char dumm[10]; int code_AA[20], ncont, a;
  for(a=0; a<20; a++)code_AA[a]=a;
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#'){
      if(strncmp(string, "#Mut", 4)==0){
	if(Get_rank(code_AA, string+10)<0)
	  for(a=0; a<20; a++)code_AA[a]=a;
      }
      continue;
    }
    if(k>=L)goto next;
    int i=res_index[k];
    if(i<0 || i>= len_amm)goto next;
    float s[20], *m=Str_mut[i]; num[i]++;
    sscanf(string,"%s%d %f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f",
	   dumm, &ncont, s, s+1, s+2, s+3, s+4, s+5, s+6, s+7, s+8, s+9,
	   s+10, s+11, s+12, s+13, s+14, s+15, s+16, s+17, s+18, s+19);
    for(i=0; i<20; i++)m[code_AA[i]]+=s[i];
  next:
    k++;
  }
  fclose(file_in);
  for(int i=0; i<len_amm; i++){
    if(num[i]>1)for(a=0; a<20; a++)Str_mut[i][a]/=num[i];
  }

  // Determine type of input
  char *s=file_str_mut; // is file name mut_DE or mut_RMSD?
  while(*s!='\0'){
    if((*s=='m')&&(*(s+1)=='u')&&(*(s+2)=='t')){
      if(*(s+4)=='D'){strcpy(STR_MUT_TYPE, "DE");}
      else if(*(s+4)=='R'){strcpy(STR_MUT_TYPE, "RMSD");}
      else{strcpy(STR_MUT_TYPE, "UNK");
	printf("WARNING, unknown file name %s\n", file_str_mut);
      }
      break;
    }
    s++;
  }
  printf("Type of mutational effects: %s\n", STR_MUT_TYPE);
  return(Str_mut);
}

int Get_rank(int *code_AA, char *string)
{
  char *AA[20]; int a; for(a=0; a<20; a++)AA[a]=malloc(8*sizeof(char));
  sscanf(string,"%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s",
	 AA[0], AA[1], AA[2], AA[3], AA[4], AA[5], AA[6], AA[7], AA[8], AA[9],
	 AA[10], AA[11], AA[12], AA[13], AA[14],
	 AA[15], AA[16], AA[17], AA[18], AA[19]);
  printf("Amino acid order: ");
  for(a=0; a<20; a++)printf("%s ",AA[a]);
  printf("\n");
  for(a=0; a<20; a++){
    int k; for(k=0; k<20; k++)if(AA[a][0]==AMIN_CODE[k])break;
    if(k==20){
      printf("WARNING, amino acid %d = %s not found in string %s\n",
	     a, AA[a], string); return(-1);
    }
    code_AA[a]=k;
  }
  return(0);
}
