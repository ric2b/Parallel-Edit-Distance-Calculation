#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

int id=0, p=0;
int *chunkN=NULL, *chunkM=NULL;
MPI_Status status;

typedef struct t_data{ /* Estrutura necessária para a computação da LCS */
	int **matrix1;
	int **matrix2; /*Matriz complementar dos slaves*/
	char *rows;
	char *columns;
	int N;
	int M;
}data;

int **createMatrix(int N, int M){ /* Cria matriz */
	int i=0;
	int **newM=NULL;
	
	newM = (int**)calloc((N+1), sizeof(int*));
	if(newM == NULL){
		printf("\n[ERROR] Creating columns for matrix\n\n");
		exit(-1);
	}
	
	for(i=0; i<(N+1); i++){
		newM[i] = (int*)calloc((M+1), sizeof(int));
		if(newM[i] == NULL){
			printf("\n[ERROR] Creating columns for matrix\n\n");
			MPI_Finalize();
			exit(-1);
		}
	}
	
	return newM;
}

char *createVector(int SIZE){ /* Cria vector para armazenamento das sequências de caractéres */
	char *newV=NULL;
	
	newV = (char*)calloc(SIZE, sizeof(char));
	if(newV == NULL){
		putchar('\n');
		printf("ERROR => Creating Vector!!!\n");
		putchar('\n');
		exit(-1);
	}
	
	return newV;
}

int *defineChunk(int DIM){
	int *aux=NULL;
	int i=0;
	
	aux=(int*)calloc(p,sizeof(int));
	if(aux==NULL){
		printf("\n[ERROR] Allocating memory for aux |Function: defineChunk(int, int)|\n\n");
		MPI_Finalize();
		exit(-1);
	}
	
	if((DIM % p) == 0){
		for(i=0; i<p; i++) aux[i] = DIM/p;
	}else{
		for(i=1; i<p; i++) aux[i] = (DIM/p)+1;
		aux[0] = (DIM - (p-1)*((DIM/p)+1));
	}
	
	return aux;
}

FILE *openFile(char *fname){ /* Abre ficheiro de input */
	FILE *fp=NULL;
	
	fp = fopen(fname, "r");
	if(fp == NULL){
		printf("\n[ERROR] Opening File\n\n");
		MPI_Finalize();
		exit(-1);
	}
	
	return fp;
}

data readFile(char *fname){ /* Lê ficheiro de input e procede à sua descodificação */
	data newData;
	FILE *pFile=NULL;
	int j=0;
	
	if(!id){
		
		pFile = openFile(fname);
		
		if(fscanf(pFile, "%d %d", &newData.N, &newData.M) != 2){
    	printf("\n[ERROR] Reading variables from the input file\n\n");
    	fclose(pFile);
    	MPI_Finalize();
    	exit(-1); 
  	}
  	if(fscanf(pFile, "\n") != 0){
  		printf("\n[ERROR] Reading EOF from the input file\n\n");
   	 	fclose(pFile);
    	MPI_Finalize();
    	exit(-1);
  	}
 	}	
 	MPI_Bcast(&(newData.N), 1, MPI_INT, 0, MPI_COMM_WORLD);
 	MPI_Bcast(&(newData.M), 1, MPI_INT, 0, MPI_COMM_WORLD);
 	
 	/*Calculate Chunks*/
 	chunkN = defineChunk(newData.N);
 	chunkM = defineChunk(newData.M);
  
  if(!id){
  	newData.rows = createVector(newData.N);
  	newData.columns = createVector(newData.M);
  }else{ /*No seu próprio Universo, os slaves apenas necessitam conhecer aquilo com o qual eles vão trabalhar*/
  	newData.rows = createVector(newData.N);
  	newData.M = chunkM[id];
  	newData.columns = createVector(newData.M);
  }
  
  if(!id){
  	int aux=0;
  
  	if(fread(newData.rows, 1, newData.N, pFile) != newData.N){
    	printf("\n[ERROR] Reading information from file\n\n");
    	free(newData.rows);
    	free(newData.columns);
    	fclose(pFile);
    	MPI_Finalize();
    	exit(-1); 
  	}
  	if(fscanf(pFile, "\n") != 0){
			printf("\n[ERROR] Reading EOF from file\n\n");
		  free(newData.rows);
		  free(newData.columns);
		  fclose(pFile);
		  MPI_Finalize();
		  exit(-1);
		}
	
		if(fread(newData.columns, 1, newData.M, pFile) != newData.M){
		  printf("\n[ERROR] Reading information from file\n\n");
		  free(newData.rows);
		  free(newData.columns);
		  fclose(pFile);
		  MPI_Finalize();
		  exit(-1); 
		}
		for(j=1; j<p; j++){
			MPI_Send(&(newData.rows[0]), newData.N, MPI_CHAR, j, j, MPI_COMM_WORLD);
			aux += chunkM[j-1]; 
			MPI_Send(&(newData.columns[aux]), chunkM[j], MPI_CHAR, j, j, MPI_COMM_WORLD);
		}
  }else{
  	MPI_Recv(&(newData.rows[0]), newData.N, MPI_CHAR, 0, id, MPI_COMM_WORLD, &status);
  	MPI_Recv(&(newData.columns[0]), newData.M, MPI_CHAR, 0, id, MPI_COMM_WORLD, &status);
  }
  
  /*Caso (N%p) != 0 || (M%p) != 0, cada slave irá ter duas matrizes de dimensões diferentes para realizar a computação*/
  /*O master, no entanto, só precisa de uma matriz: a matriz completa de dimensões NxM*/
  if(!id){
  	newData.matrix1 = createMatrix(newData.N, newData.M);
  	fclose(pFile);
  }else{
  	newData.N = chunkN[id];
  	newData.matrix1 = createMatrix(newData.N, newData.M);
  	newData.matrix2 = createMatrix(chunkN[0], newData.M);
  }
  
	return newData;
}

int main(int argc, char *argv[]){
	int tid=0;
	double start=0, end=0;
	char *fname=NULL;
	data lcs;
	
	start = MPI_Wtime();
	MPI_Init(&argc, &argv);
	
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	
	if(argc != 2){
		if(!id) printf("\n[ERROR] Missing initial arguments\n\n");
		MPI_Finalize();
		exit(0);
	}
	
	/*Leitura do ficheiro / Cálculo Chunk*/
	if(!id) fname = argv[1];
	
	lcs = readFile(fname);
	
	
	/*Critical section for printing information after reading the file*/
	/*tid=0;
	while(tid < p){
		if(id == tid){
			printf("(id=%d): \n", id);
			printf("N=%d, M=%d\n", lcs.N, lcs.M);
			printf("Rows: %s\n", lcs.rows);
			printf("Columns: %s\n", lcs.columns);
		}
		tid++;
		MPI_Barrier(MPI_COMM_WORLD);
	}*/
	
	MPI_Finalize();
	end = MPI_Wtime();
	
	if(!id) printf("\nTime: %.5g\n\n", (end-start));
	
	exit(0);
}
