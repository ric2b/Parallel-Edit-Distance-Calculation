#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#define max(a,b) (((a) > (b)) ? (a) : (b)) /*Retorna o máximo entre A e B*/
#define min(a,b) (((a) < (b)) ? (a) : (b)) /*Retorna o mínimo entre A e B*/

int id=0, p=0;
int *chunkN=NULL, *chunkM=NULL;
MPI_Status status;
int gN=0, gM=0;

typedef struct t_data{ /*Estrutura necessária para a computação da LCS*/
	int **matrix;
	char *rows;
	char *columns;
}data;

int **createMatrix(int N, int M){ /*Cria matriz como um vector contínuo de elementos*/
	int i=0;
	int *aux=NULL;
	int **newM=NULL;
	
	aux = (int*)calloc((N+1)*(M+1), sizeof(int));
	if(aux == NULL){
		printf("\n[ERROR] Creatingmatrix\n\n");
		MPI_Finalize();
		exit(-1);
	}
	
	newM = (int**)calloc((N+1), sizeof(int*));
	if(newM == NULL){
		printf("\n[ERROR] Creating matrix\n\n");
		MPI_Finalize();
		exit(-1);
	}
	
	for(i=0; i<(N+1); i++){
		newM[i] = &(aux[(M+1)*i]);
	}
	
	return newM;
}

char *createVector(int SIZE){ /* Cria vector para armazenamento das sequências de caractéres */
	char *newV=NULL;
	
	newV = (char*)calloc(SIZE, sizeof(char));
	if(newV == NULL){
		printf("\n[ERROR] Creating Vector\n\n");
		exit(-1);
	}
	
	return newV;
}

int *defineChunk(int DIM){
	int *aux=NULL;
	int i=0;
	
	aux=(int*)calloc(p,sizeof(int));
	if(aux==NULL){
		printf("\n[ERROR] Allocating memory for aux\n\n");
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
		
		if(fscanf(pFile, "%d %d", &gN, &gM) != 2){
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
 	MPI_Bcast(&gN, 1, MPI_INT, 0, MPI_COMM_WORLD);
 	MPI_Bcast(&gM, 1, MPI_INT, 0, MPI_COMM_WORLD);
 	
 	if((gN/p) <= 1.5 || (gM/p) <= 1.5){
 		p=min(gN,gM)/1.5;
 		if(!id) printf("\n[WARNING] The program will only run with %d machines\n\n", p);
 	}
 	
 	if(id < p){
	 	/*Calculate Chunks*/
	 	chunkN = defineChunk(gN);
	 	chunkM = defineChunk(gM);
		
		if(!id){
			int aux=0;
			char *rows=NULL;
			
			rows = createVector(gN);
		
			if(fread(rows, 1, gN, pFile) != gN){
		  	printf("\n[ERROR] Reading information from file\n\n");
		  	free(rows);	
		  	fclose(pFile);
		  	MPI_Finalize();
		  	exit(-1); 
			}
			if(fscanf(pFile, "\n") != 0){
				printf("\n[ERROR] Reading EOF from file\n\n");
				free(rows);
				fclose(pFile);
				MPI_Finalize();
				exit(-1);
			}
	
			newData.columns = createVector(gM);
		
			if(fread(newData.columns, 1, gM, pFile) != gM){
				printf("\n[ERROR] Reading information from file\n\n");
				free(rows);
				free(newData.columns);
				fclose(pFile);
				MPI_Finalize();
				exit(-1); 
			}
			for(j=0; j<p; j++){ /*Cada máquina apenas necessita parte da string correspondente às linhas - incluíndo o "Master"*/
				aux+=chunkN[j-1];
				MPI_Send(&(rows[aux]), chunkN[j], MPI_CHAR, j, j, MPI_COMM_WORLD);	
			}
			free(rows);
		}
		
		newData.rows = createVector(chunkN[id]);
		if(id){
			newData.columns = createVector(gM);
		}
		
		MPI_Recv(&(newData.rows[0]), chunkN[id], MPI_CHAR, 0, id, MPI_COMM_WORLD, &status);
		/*Todas as máquinas precisam de toda a string correspondente às colunas*/
		MPI_Bcast(&(newData.columns[0]), gM, MPI_CHAR, 0, MPI_COMM_WORLD);
		
		/*Cada máquina terá a sua própria matriz de dimensões (chunkN[id]+1)x(M+1)*/
		newData.matrix = createMatrix(chunkN[id], gM);
  
  	if(!id) fclose(pFile);
  }
  
	return newData;
}

short cost(int x){ /* Função cost fornecida pelos professores */
	int i=0, n_iter=20;
	double dcost=0;
	
	for(i=0; i<n_iter; i++){
		dcost += pow(sin((double) x), 2) + pow(cos((double) x), 2);
	}	
	
	return (short) (dcost / n_iter + 0.1);	
}

int **computeMatrix(data lcs, int N, int M, int mj){ /* Realiza a computação do bloco da matriz */
	int i=0, jj=0, j=0;
	
	for(i=1; i<(N+1); i++){
		for(jj=1; jj<(M+1); jj++){
			j=jj+mj;
			if(lcs.rows[i-1] == lcs.columns[j-1]){
				lcs.matrix[i][j] = lcs.matrix[i-1][j-1] + cost(i); /*funcão cost*/
			}else{
				lcs.matrix[i][j] = max(lcs.matrix[i][j-1], lcs.matrix[i-1][j]);
			}
		}
	}
	
	return lcs.matrix;
}

data distributeMatrix(data lcs){
	int k=0, nj=0, idN=0, idM=0, aux=0;
	int i=0;
	
	for(k=0; k < (2*p)-1; k++){
		if(k < p){ /*Aumento da carga de trabalho*/ 
			if(!id){ /*O processo 0 nunca recebe informação de outro processo*/
				idN=chunkN[id];
				idM=chunkM[k];
				
				lcs.matrix = computeMatrix(lcs, idN, idM, nj); /*Computação do bloco correspondente*/
				
				MPI_Send(&(lcs.matrix[idN][nj+1]), idM, MPI_INT, id+1, id, MPI_COMM_WORLD); /*Envia informação para o processo id+1*/
		
				nj += idM; /*Indica a primeira coluna do bloco seguinte*/
			}else{
				if(k >= id){ /*O primeiro bloco a ser processado, por cada processo, poderá ter dimensões diferentes dos próximos blocos*/
					idN=chunkN[id];
					idM=chunkM[k-id];
					
					/*Recebe informação necessária do processo id-1*/
					MPI_Recv(&(lcs.matrix[0][nj+1]), idM, MPI_INT, id-1, id-1, MPI_COMM_WORLD, &status);
					
					lcs.matrix = computeMatrix(lcs, idN, idM, nj); /*Computação do bloco correspondente*/
					
					/*Envia informação para o processo id+1 - o processo id=p-1 nunca envia informação para outro processo*/
					if(id != p-1) MPI_Send(&(lcs.matrix[idN][nj+1]), idM, MPI_INT, id+1, id, MPI_COMM_WORLD);
					
					nj += idM;
				}
			}
		}else{ /*Diminuição da carga de trabalho*/
			if(id > k-p){
				idN=chunkN[id];
				idM=chunkM[k-id];
				
				/*Recebe informação necessária do processo id-1*/
				MPI_Recv(&(lcs.matrix[0][nj+1]), idM, MPI_INT, id-1, id-1, MPI_COMM_WORLD, &status);
				
				lcs.matrix = computeMatrix(lcs, idN, idM, nj); /*Computação do bloco correspondente*/
				
				/*Envia informação para o processo id+1 - o processo id=p-1 nunca envia informação para outro processo*/
				if(id != p-1) MPI_Send(&(lcs.matrix[idN][nj+1]), idM, MPI_INT, id+1, id, MPI_COMM_WORLD);
					
				nj += idM;
			}
		}
	}
	
	return lcs;
}

void freeMem(data lcs){
	
	/*Libertart memória alocada para a matriz de cada máquina*/
	free(lcs.matrix[0]);
	free(lcs.matrix);
	
	/*Libertar memória alocada para os vectores de caractéres de cada máquina*/
	free(lcs.rows);
	free(lcs.columns);
	
	return;
}

int main(int argc, char *argv[]){
	int tid=0;
	/*int i=0;*/
	int i=0, j=0;
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
	
	/*Leitura do ficheiro / Cálculo do Chunk relativo a cada Máquina*/
	if(!id) fname = argv[1];
	
	lcs = readFile(fname);
	
	/*Computação da matriz*/
	if(id < p) lcs = distributeMatrix(lcs);
	
	/*Critical Section for printing the whole Matrix*/
	tid=0;
	while(tid < p){
		if(id < p){
			if(id == tid){
				if(!id){
					for(j=0; j<gM+1; j++){
						printf("%d ", lcs.matrix[0][j]);			
					}
					putchar('\n');
				}
				for(i=1; i<chunkN[id]+1; i++){
					for(j=0; j<gM+1; j++){
						printf("%d ", lcs.matrix[i][j]);			
					}
					putchar('\n');
				}
			}
		}
		tid++;
		MPI_Barrier(MPI_COMM_WORLD);
	}
	
	if(id < p) freeMem(lcs);
	
	MPI_Finalize();
	end = MPI_Wtime();
	
	if(!id) printf("\nTime: %.5g\n\n", (end-start));
	
	exit(0);
}
