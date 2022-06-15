#include "codes.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "read_ali.h"

static int Align_nogap(int *ali_PDB, char *aa_seq, int L_PDB,
		       char *ali, int L_ali, int nalimin);

int Read_ali(float **freq_ia, float *Seq_id, short ***ali_seq, int *N_ali,
	     char *file_ali, short *i_seq, int L_seq_PDB, int L_str_PDB)
{
  if(file_ali==NULL){
    printf("WARNING, there is no name of FASTA file with alignment\n");
    return(0);
  }
  FILE *file_in=fopen(file_ali, "r");
  if(file_in==NULL){
    printf("WARNING, FASTA file with alignment %s does not exist\n",
	   file_ali); return(0);
  }

  int n=0, l=0; char string[10000], *s;
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='>'){n++; continue;}
    else if(n==1){
      s=string; while((*s!='\n')&&(*s!='\r')&&(*s!='\0')){l++; s++;}
    }
  }
  fclose(file_in);
  printf("%d sequences of length %d read in alignment %s\n",n,l,file_ali);

  // Allocate
  int N_seq=n, L_ali=l;
  char *ali[N_seq], *ali_n;
  for(n=0; n<N_seq; n++)ali[n]=malloc(L_ali*sizeof(char));

  file_in=fopen(file_ali, "r"); n=-1;
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='>'){
      if((n>=0)&&(l!=L_ali)){
	printf("ERROR, variable number of columns in alignment\n");
	printf("Sequence %d has %d symbols instead of %d\n",n,l,L_ali);
	printf("Just read: %s\n", string);
	for(n=0; n<N_seq; n++)free(ali[n]);
	return(0);
      }
      n++; ali_n=ali[n]; l=0; continue;
    }
    s=string; 
    while((*s!='\n')&&(*s!='\r')&&(*s!='\0')){
      if(l>=L_ali){
	printf("ERROR reading alignment, sequence %d has %d > %d symbols\n",
	       n, l, L_ali); printf("Just read: %s\n", string);
	for(n=0; n<N_seq; n++)free(ali[n]);
	return(0);
      }
      *ali_n=*s; s++; l++; ali_n++;
    }
  }
  fclose(file_in);
  printf("Read %d sequences with length %d\n", N_seq,L_ali);

  // Align to PDB sequence
  int ali_PDB[L_seq_PDB], gapmax=10;
  char aa_seq[L_seq_PDB]; int nalimax=0, nalimin=L_str_PDB-gapmax;
  for(int i=0; i<L_seq_PDB; i++)aa_seq[i]=AMIN_CODE[i_seq[i]];
  for(n=0; n<N_seq; n++){
    int nali=
      Align_nogap(ali_PDB, aa_seq, L_seq_PDB, ali[n], L_ali, nalimin);
    if(nali > nalimax)nalimax=nali;
    if(nali >= nalimin){break;}
  }
  if(n==N_seq){
    printf("WARNING, the PDB sequence could not be aligned to any seq\n");
    printf("Maximum alignment length: %d\n", nalimax);
    for(n=0; n<N_seq; n++)free(ali[n]);
    //exit(8);
    return(0);
  }
  int npdb=n;
  printf("PDB sequence aligned to sequence %d nali= %d\n",n,nalimax);
  printf("Seq PDB: ");
  for(l=0; l<L_seq_PDB; l++)printf("%c",aa_seq[l]);
  printf("\n");
  printf("Seq ali: ");
  for(l=0; l<L_seq_PDB; l++){
    if(ali_PDB[l]<0){printf("-");}
    else{printf("%c", ali[npdb][ali_PDB[l]]);}
  }
  printf("\n");

  *ali_seq=malloc(N_seq*sizeof(short *));
  for(n=0; n<N_seq; n++){
    (*ali_seq)[n]=malloc(L_seq_PDB*sizeof(short));
  }
  int a;
  for(int i=0; i<L_seq_PDB; i++){
    for(a=0; a<20; a++)freq_ia[i][a]=0;
    int j=ali_PDB[i];
    for(n=0; n<N_seq; n++){
      if((j>=0)&&(ali[n][j]!='-')&&(ali[n][j]!='.')){
	a=Code_AA(ali[n][j]); 
	if((a<0)||(a>19)){printf("ERROR wrong AA %c\n", ali[n][j]); a=-1;}
	else{freq_ia[i][a]++;}
	(*ali_seq)[n][i]=a;
      }else{
	(*ali_seq)[n][i]=-1;
      }
    }
    float sum=0; for(a=0; a<20; a++)sum+=freq_ia[i][a];
    if(sum==0){
      printf("WARNING, column %d corresponding to site %d is empty\n",j,i); 
    }
    // Mean sequence identity
    *Seq_id=0;
    for(n=0; n<N_seq; n++){
      if(n==npdb)continue;
      int s=0;
      for(int j=0; j<L_ali; j++){
	if(ali[n][j]==ali[npdb][j] && ali[n][j]!='-')s++;
      } 
      /*for(int i=0; i<L_seq_PDB; i++){
	int j=ali_PDB[i]; if(j<0)continue;
	if(ali[n][j]==aa_seq[i])s++;
	}*/
      *Seq_id+=s;
    }
    *Seq_id/=(L_seq_PDB*(N_seq-1));
  }
  *N_ali=L_ali;
  return(N_seq);

}

int Align_nogap(int *ali_PDB, char *aa_seq, int L_PDB,
		char *ali, int L_ali, int nalimin)
{
  char *a1, *a2; int L1, i1, L2, i2;
  int DBG=0;
  int nali=0, gap=0, start=0;
  a1=aa_seq; L1=L_PDB;
  a2=ali; L2=L_ali; i2=0;
  if(DBG){
    for(i1=0; i1<L_ali; i1++)printf("%c", ali[i1]);
    printf("\n");
  }
  for(i1=0; i1<L_PDB; i1++)ali_PDB[i1]=-1;
  for(i1=0; i1<L1; i1++){
    while(a2[i2]=='-' && i2<L2)i2++;
    if(i2==L2)break;
    if(a1[i1]!=a2[i2]){
      if(start==0)continue; // Initial part of PDB sequence, it may be 
                            // disordered and some programs do not read it
      if(a2[i2]==a1[i1+1] && gap==0){gap=1; continue;}
      if(DBG && nali>100)printf("%c %c %d %d\n", a1[i1],a2[i2],i1,nali);
      while((a2[i2]!=a1[i1])&&(i2<L2))i2++;
      if(i2==L2)break;
    }
    if(start==0)start=1;
    if(a1[i1]==a2[i2]){
      ali_PDB[i1]=i2; i2++; nali++; if(gap)gap=0;
    }else if(gap){
      break;
    }
  }
  if(DBG)printf("nali= %d\n", nali);
  if(nali>=nalimin)return(nali);
  
  int nali2=0; a1=ali; L1=L_ali; gap=0;
  a2=aa_seq; L2=L_PDB; i2=0; nali2=0;
  for(i1=0; i1<L_PDB; i1++)ali_PDB[i1]=-1;
  for(i1=0; i1<L1; i1++){
    if(a1[i1]!=a2[i2]){
      if(a1[i1]=='-')continue;
      if(a2[i2]==a1[i1+1] && gap==0){gap=1; continue;}
      if(DBG && nali>100)printf("%c %c %d %d\n",a2[i2],a1[i1],i2,nali2);
      while((a2[i2]!=a1[i1])&&(i2<L2))i2++;
      if(i2==L2)break;
    }
    if(a1[i1]==a2[i2]){
      ali_PDB[i2]=i1; gap=0;
      int j1=i1-1;
      if(i2 && ali_PDB[i2-1]<j1 && a2[i2-1]==a1[j1]){
	ali_PDB[i2-1]=j1;
	for(int j2=i2-2; j2>=0; j2--)
	  if(ali_PDB[j2]==j1){ali_PDB[j2]=-1; break;}
      }
      i2++; nali2++;
    }else if(gap){
      break;
    }
  }
  if(DBG)printf("nali2= %d\n", nali2);

  if(nali>nali2){return(nali);}
  else{return(nali2);}
}
