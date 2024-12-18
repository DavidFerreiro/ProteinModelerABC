#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define N_CHAR 400
#include "codes.h"
#include "get_para_DeltaGREM.h"

static void help(char *);
static int
Read_parameters(char *FILE_IN, char **file_pdb, char chain[], char *file_str,
		float *TEMP, float *sU1, float *sC0, float *sC1,
		int *REM, float *SEC_STR,
		char **file_ali, char **file_mut, char *dir_out);
static int Find_string(char *flag, char **argv, int n_arg,int m,char *string);


int Get_para(int argc, char **argv,
	     char ***file_pdb, char chain[], char *file_str,
	     float *TEMP, float *sU1, float *sC0, float *sC1,
	     int *REM, float *SEC_STR,
	     char **file_ali, char **file_mut, char *dir_out)
{
  int PDBMAX=500;
  char FILE_IN[N_CHAR], string[N_CHAR];
  int N_pdb, i, j; float p;

  if(Find_string("-h", argv, argc, 2, string))help(argv[0]);
  if(argc < 2){
    printf("ERROR, either input file or options must be specified\n");
    help(argv[0]);
  }

  strcpy(FILE_IN, argv[1]);
  *file_pdb=malloc(PDBMAX*sizeof(char *));
  N_pdb=Read_parameters(FILE_IN, *file_pdb, chain, file_str,
			TEMP, sU1, sC0, sC1, REM, SEC_STR,
			file_ali, file_mut, dir_out);

  /*********************************************/
  if(argc <= 2){
    if(N_pdb==-1){
      printf("ERROR, input file %s does not exist\n", argv[1]);
      help(argv[0]);
    }else{
      goto check;
    }
  }
  printf("! Changing parameters with command line:\n");

  for(j=1; j<argc; j++){
    if(strncmp(argv[j], "-pdb", 4)==0){
      j++; sscanf(argv[j], "%s", (*file_pdb)[0]);
    }else if(strncmp(argv[j], "-chain", 3)==0){
      j++; sscanf(argv[j], "%s", chain);
    }else if(strncmp(argv[j], "-temp", 5)==0){
      sscanf(argv[j+1], "%f", &p);
      if(p>0){
	*TEMP=p; /*printf("Temperature: %3f\n", p);*/ j++;
      }
    }else if(strncmp(argv[j], "-sU1", 3)==0){
      sscanf(argv[j+1], "%f", &p);
      if(p>0){
	*sU1=p; printf("Unfolded entropy per residue: %3f\n", *sU1); j++;
      }
    }else if(strncmp(argv[j], "-ali", 4)==0){
      *file_ali=malloc(200*sizeof(char));
      sscanf(argv[j+1], "%s", *file_ali);
    }else if(strncmp(argv[j], "-d", 2)==0){
      strcpy(dir_out,argv[j+1]); j++;
      printf("Directory in is: %s\n",dir_out);
    }else{
      printf("WARNING, unrecognized option %s\n", argv[j]);
    }
  }

  printf("! End changing parameters\n\n");

  /*********************************************/

  // Checks
 check:
  if(N_pdb<=0){
      printf("ERROR, PDB file not specified\n");
      help(argv[0]);
  }
  for(i=0; i<N_pdb; i++){
    if((*file_pdb)[i][0]=='\0'){
      printf("ERROR, PDB file not specified\n");
      help(argv[0]);
    }
  }
  return(N_pdb);
}

int Read_parameters(char *FILE_IN, char **file_pdb, char *chain, char *file_str,
		    float *TEMP, float *sU1, float *sC0, float *sC1, 
		    int *REM, float *SEC_STR,
		    char **file_ali, char **file_mut, char *dir_out)
{
  FILE *file_in=fopen(FILE_IN, "r");
  char string[1000], dumm[80]; int N_pdb=0; float x;

  if(file_in==NULL){
    printf("WARNING, input file %s does not exist\n", FILE_IN);
    printf("Revise also the output directory %s\n", dir_out); // Miguel (removing warning)
    return(-1);
  }
  /*printf("Reading parameters in %s\n", FILE_IN);*/

  char READ[100];
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#')continue;
    // Protein files
    if(strncmp(string, "PDB=", 4)==0){
      file_pdb[N_pdb]=malloc(100*sizeof(char));
      sscanf(string+4,"%s", file_pdb[N_pdb]);
      N_pdb++;
    }else if(strncmp(string, "CHAIN=", 6)==0){
      sscanf(string+6,"%s", dumm); *chain=dumm[0];
      /*printf("chain= %c\n", *chain);*/
    }else if(strncmp(string, "ALI=", 4)==0){
      sscanf(string+4,"%s", READ);
      if(strncmp(READ, "NULL", 4)==0)continue;
      *file_ali=malloc(100*sizeof(char));
      strcpy(*file_ali, READ);
    }else if(strncmp(string, "MUT=", 4)==0){
      sscanf(string+4,"%s", READ);
      if(strncmp(READ, "NULL", 4)==0)continue;
      *file_mut=malloc(100*sizeof(char));
      strcpy(*file_mut, READ);
      // Thermodynamic parameters:
    }else if(strncmp(string, "TEMP=", 5)==0){
      sscanf(string+5,"%f", TEMP);
    }else if(strncmp(string, "SU1=", 3)==0){
      sscanf(string+4,"%f", sU1); 
    }else if(strncmp(string, "SC1=", 4)==0){
      sscanf(string+4,"%f", sC1); 
    }else if(strncmp(string, "SC0=", 4)==0){
      sscanf(string+4,"%f", sC0); 
    }else if(strncmp(string, "REM=", 4)==0){
      sscanf(string+4,"%d", REM);
      if((*REM!=1)&&(*REM!=2)&&(*REM!=3)){
	printf("ERROR, the variable REM can only be set to 1, 2 or 3\n");
	exit(8);
      }
    }else if(strncmp(string, "SEC_STR=", 8)==0){
      sscanf(string+8,"%f", &x);
      if(x>=0)*SEC_STR=x;
    }else if(strncmp(string, "FILE_STR=", 9)==0){
      sscanf(string+9,"%s", file_str);
    }else{
      /*printf("WARNING, uninterpreted line:\n%s", string); */
    }
  }
  fclose(file_in);
  /*printf("Parameters: T= %.2f SC/L= %.3f SU/L= %.3f REM= %d\n",
	 *TEMP, *sC1, *sU1, *REM);*/
  return(N_pdb);
}

int Find_string(char *flag, char **argv, int n_arg, int m, char *string){
  int i;
  for(i=0; i<n_arg; i++){
    if(strncmp(argv[i], flag, m)==0){
      i++; if(i<n_arg)strcpy(string, argv[i]);
      return(1);
    }
  }
  return(0);
}

void help(char *prog){
  printf("USAGE %s <input file> ", prog);
  printf("(ex. Prot_evol.in) with default parameters\n");
  printf("ARGUMENTS:\n");
  printf("Protein:\n");
  printf("-pdb <file.pdb>     # Mandatory, unless in input file\n");
  printf("-chain <chain>      # Chain identifier (ex. A)\n");
  printf("-seq  <file_seq>    # file with DNA sequece (optional)\n");
  printf("-ali <file>         # file with protein seq in FASTA format (opt)\n");

  printf("\n\nFORMAT of input file:\n");
  printf("#================================================================\n");
  printf("# A) Input files\n");
  printf("PDB=/data/ortizg/databases/pdb/1tre.pdb  # file_pdb\n");
  printf("CHAIN=  A\n");
  printf("# (one or many pdb files may be provided)\n"); 
  printf("FILE_STR=/home/ubastolla/RESEARCH/PROT_EVOL/INPUT/structures.in\n");
  printf("# List of alternative contact matrices for misfolding computations\n");
  printf("ALI=	family.aln	# FASTA file with proteins to examine (optional)\n");
  printf("MUT=	prot.mut	# list of mutations (optional)\n");
  printf("#================================================================\n");
  printf("# B) Thermodynamic model\n");
  printf("TEMP=	0.5		# Temperature\n");
  printf("SU1=	0.065		# configurational entropy per res (unfold)\n");
  printf("SC1=  0.065		# configurational entropy per res (misfold)\n");
  printf("SC0=  0.0		# configurational entropy offset (misfold)\n");
  printf("REM=   2		# Use up to 1,2,3r moments of misfolding energy?\n");
  printf("\n");
  printf("FORMAT of mutation file:\n");
  printf("All mutations that occur together go in the same line\n");
  printf("If more than one chain is present, the chain must be specified");
  printf(" after the mutation.\n");
  printf("If single chain the chain must NOT be specified. Example:\n");
  printf("T103A A T103A B (polychain)\n");
  printf("T103A (single chain)\n");
  printf("FORMAT of deletion:\n");
  printf("DEL T103-W109 A T103A-W109 B (polychain)\n");
  printf("DEL T103-W109 (single chain)\n");

  printf("\n");
  exit(8);
}
