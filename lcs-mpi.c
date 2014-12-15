/***********************************************
* Project: Longest Common Subsequence          *
*                                              *
* Course: Parallel Distributed Computing       *
*	                                             *
* Authors: Daniel Arrais, nº 69675             *
*          Miguel Costa, nº 73359              *
*          Ricardo Amendoeira, nº 73373			   *
*                                              *
* Version: MPI                                 *
************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define max(a,b) (((a) > (b)) ? (a) : (b)) /*Retorna o máximo entre A e B*/

int id=0, p=0;
int gN=0, gM=0, sizeLCS=0;
int *chunkN=NULL;
MPI_Status status;

/*Estrutura definida como pilha - armazena a LCS*/
typedef struct t_stack{ 
	char *st;
	int top;
}stack;

/*Aloca memória necessária para a estrutura pilha*/
stack createStack(int sizeStack){
	stack new;

	new.st = (char*)calloc(sizeStack, sizeof(char));
	if(new.st == NULL){
		printf("\n[ERROR] Creating Stack\n\n");
		exit(-1);
	}
	
	new.top = -1;
	
	return new;
}

/*Realoca memória necessária para a estrutura pilha para cada processo, excepto o 0*/
stack reallocStack(stack seq, int count){

	seq.st = (char*)realloc(seq.st, count*sizeof(char));
	if(seq.st == NULL){
		printf("\n[ERROR] Reallocating memory for stack\n\n");
		MPI_Finalize();
		exit(-1);
	}
	
	return seq;
}

/*Insere novo elemento na pilha*/
stack push(stack seq, char item, int x){
  seq.top=x;
  
  seq.st[seq.top] = item;
  
	return seq;
}

/*Imprime os caractéres armazenados na pilha*/
void display(stack seq){ 
  int i=0;
   
	for(i=seq.top; i>=0; i--){
  	printf("%c", seq.st[i]);
	}
	putchar('\n');
  
  return;
}

/*Estrutura necessária para a computação da matriz*/
typedef struct t_data{
	short **matrix;
	char *rows;
	char *columns;
}data;

/*Cria uma matriz como um vector contínuo de elementos*/
short **createMatrix(int N, int M){
	int i=0;
	short *aux=NULL;
	short **newM=NULL;
	
	aux = (short*)calloc((N+1)*(M+1), sizeof(short));
	if(aux == NULL){
		printf("\n[ERROR] Allocating memory in function createMatrix(int,int)\n\n");
		MPI_Finalize();
		exit(-1);
	}
	
	newM = (short**)calloc((N+1), sizeof(short*));
	if(newM == NULL){
		printf("\n[ERROR] Allocating memory in function createMatrix(int,int)\n\n");
		MPI_Finalize();
		exit(-1);
	}
	
	for(i=0; i<(N+1); i++){
		newM[i] = &(aux[(M+1)*i]);
	}
	
	return newM;
}

/*Cria vector para armazenamento das sequências de caractéres*/
char *createVector(int SIZE){ 
	char *newV=NULL;
	
	newV = (char*)calloc(SIZE, sizeof(char));
	if(newV == NULL){
		printf("\n[ERROR] Allocating memory in function createVector(int)\n\n");
		MPI_Finalize();
		exit(-1);
	}
	
	return newV;
}

/*Calcula o tamanho do bloco correspondente a cada processo*/
int *defineChunk(int DIM){
	int *aux=NULL;
	int i=0;
	
	aux=(int*)calloc(p,sizeof(int));
	if(aux==NULL){
		printf("\n[ERROR] Allocating memory in function defineChunk(int)\n\n");
		MPI_Finalize();
		exit(-1);
	}
	
	if((DIM % p) == 0){ /*Se DIM/p for uma divisão inteira, cada processo fica com o mesmo chunk*/
		for(i=0; i<p; i++) aux[i] = DIM/p; 
	}else{ /*Caso contrário, cada processo fica com (DIM/p)+1 excepto o processo p-1 que fica com o restante*/
		for(i=0; i<(p-1); i++) aux[i] = (DIM/p)+1;
		aux[p-1] = (DIM - (p-1)*((DIM/p)+1));
	}
	
	return aux;
}

/*Abre ficheiro de input*/
FILE *openFile(char *fname){
	FILE *fp=NULL;
	
	fp = fopen(fname, "r");
	if(fp == NULL){
		printf("\n[ERROR] Opening File\n\n");
		MPI_Finalize();
		exit(-1);
	}
	
	return fp;
}

/*Leitura do ficheiro de input e consequente descodificação*/
data readFile(char *fname){
	FILE *pFile=NULL;
	data newData;
	
	if(!id){ /*Apenas o processo 0 tem acesso ao ficheiro*/
		pFile = openFile(fname);
	
		if(fscanf(pFile, "%d %d", &gN, &gM) != 2){ /*Leitura do tamanho das strings (N,M)*/
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
	/*Todos os processos têm conhecimento do tamanho total das strings de caractéres*/
	MPI_Bcast(&gN, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&gM, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	chunkN=defineChunk(gN); /*Todos os processos têm conhecimento dos chunks designados a cada processo*/
	
	if(!id){
		int aux=0, j=0;
		char *rows=NULL;
		
		rows = createVector(gN);
	
		if(fread(rows, 1, gN, pFile) != gN){ /*Leitura da string de tamanho N*/
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
	
		if(fread(newData.columns, 1, gM, pFile) != gM){ /*Leitura da string de tamanho M*/
			printf("\n[ERROR] Reading information from file\n\n");
			free(rows);
			free(newData.columns);
			fclose(pFile);
			MPI_Finalize();
			exit(-1); 
		}
		for(j=0; j<p; j++){ /*Cada processo apenas necessita de parte da string de tamanho N*/
			MPI_Send(&(rows[aux]), chunkN[j], MPI_CHAR, j, j, MPI_COMM_WORLD); /*Envio da informação correspondente a cada*/
			aux+=chunkN[j];
		}
		free(rows);
	}
	
	newData.rows = createVector(chunkN[id]); /*Cria vector de tamanho chunkN[id] - armazena parte da string de tamanho N*/
	if(id) newData.columns = createVector(gM); /*Cria vector de tamanho M - armazena toda a string de tamanho M*/
	
	/*Cada processo só precisa de parte da string correspondente às linhas*/
	MPI_Recv(&(newData.rows[0]), chunkN[id], MPI_CHAR, 0, id, MPI_COMM_WORLD, &status);
	/*Todos os processos precisam da totalidade da string correspondente às colunas*/
	MPI_Bcast(&(newData.columns[0]), gM, MPI_CHAR, 0, MPI_COMM_WORLD);
	
	/*Cada processo terá a sua própria matriz de dimensões (chunkN[id]+1)x(M+1)*/
	newData.matrix = createMatrix(chunkN[id], gM);
	
	if(!id) fclose(pFile);
	
	return newData;
}

/*Função cost fornecida*/
short cost(int x){
	int i=0, n_iter=20;
	double dcost=0;
	
	for(i=0; i<n_iter; i++){
		dcost += pow(sin((double) x), 2) + pow(cos((double) x), 2);
	}	
	
	return (short) (dcost / n_iter + 0.1);	
}

/*Computação da matriz*/
data computeMatrix(data lcs){
	int i=0, j=0;
	MPI_Request *request=NULL;
	
	request = (MPI_Request*)calloc(gM, sizeof(MPI_Request)); /*Cria um vector de MPI_Requests - usado no MPI_Isend*/
	for(j=1; j<gM+1; j++){ /*Cada processo percorre a sua matriz ao longo das colunas*/
		if(id){ /*O processo 0 nunca recebe informação de outro processo*/
			/*Recepção de informação necessária para a computação da próxima coluna*/
			MPI_Recv(&(lcs.matrix[0][j]), 1, MPI_SHORT, id-1, id, MPI_COMM_WORLD, &status);
		}
		for(i=1; i<chunkN[id]+1; i++){
			if(lcs.rows[i-1] == lcs.columns[j-1]){
				lcs.matrix[i][j] = lcs.matrix[i-1][j-1] + cost(i);
			}else{
				lcs.matrix[i][j] = max(lcs.matrix[i][j-1], lcs.matrix[i-1][j]);
			}
		}
		if(id < p-1){ /*O processo p-1 nunca envia informação para outro processo*/
			/*Após a computação de cada coluna, cada processo id envia para o processo id+1 o último elemento da coluna computada*/
			MPI_Isend(&(lcs.matrix[i-1][j]), 1, MPI_SHORT, id+1, id+1, MPI_COMM_WORLD, &request[j-1]);
		}
	}
	
	return lcs;
}

/*Retorna a LCS, realizando um traceback na matriz*/
stack computeLCS(data lcs, stack seq, int *flag){
	MPI_Request request;
	int i=0, j=0, k=0, count=0, t=0;
	int *ij=NULL;
	
	if(id == p-1){
		i=chunkN[id];
		j=gM;
		sizeLCS = lcs.matrix[i][j];
		/*O processo p-1 envia para o processo 0 o actual tamanho da LCS*/
		MPI_Isend(&sizeLCS, 1, MPI_INT, 0, id, MPI_COMM_WORLD, &request);
	}
	
	if(!id){
		MPI_Recv(&sizeLCS, 1, MPI_INT, p-1, p-1, MPI_COMM_WORLD, &status);
		seq = createStack(sizeLCS); /*O processo 0 aloca memória necessária para a LCS*/
	}
	
	ij=(int*)calloc(2, sizeof(int));
	if(ij == NULL){
		printf("\n[ERROR] Allocating memory in function computeLCS(data,stack,int)\n\n");
		MPI_Finalize(); 
		exit(-1);
	}
	ij[0]=gN;
	ij[1]=gM;
	
	k=p;
	
	while(ij[0] != 0 && ij[1] != 0){ /*Enquanto nenhum processo atingir a coluna 0 ou a linha 0*/
		
		k--;
		count=0;
		
		if(k == id){ /*Cada processo vai sendo chamado decrescentemente: o processo p-1 é o primeiro, depois o p-2, e assim sucessivamente*/
			i=chunkN[id];
			j=ij[1];
			while(i != 0 && j != 0){ /*Enquanto não atingir a linha i=0 ou a coluna j=0 da sua matriz*/
				if(lcs.rows[i-1] == lcs.columns[j-1]){ /*Caso haja um match*/
					count++;
					(*flag)=1;
					
					if(id){ /*À excepção do processo 0, os restantes vão alocando memória para a pilha à medida que necessitam*/
						if(count==1) seq = createStack(count);
						if(count>1)  seq = reallocStack(seq, count);
					}
					
					if(!id && count==1) count+=t;
					
					seq=push(seq, lcs.rows[i-1], count-1); /*Insere novo elemento na pilha*/
					
					i=i-1;
					j=j-1;
				}else{
					if(lcs.matrix[i][j-1] >= lcs.matrix[i-1][j]){
						j=j-1;
					}else{
						i=i-1;
					}
				}
			}
			if(!id && i==0) ij[0]=0; /*Apenas o processo 0 pode obter i=0*/
			ij[1]=j; /*Ponto de partida para o próximo processo*/
			if(id){ /*Cada processo, após calcular a parte da LCS correspondente ao seu bloco, envia esta para o processo 0*/
				MPI_Send(&count, 1, MPI_INT, 0, k, MPI_COMM_WORLD);
				MPI_Send(&(seq.st[0]), count, MPI_CHAR, 0, k, MPI_COMM_WORLD);
			}
		}
		
		if(!id){
			if(ij[0] != 0 && ij[1] != 0){ /*O processo 0 vai concatenando as LCS calculadas por cada um dos processos*/
				MPI_Recv(&count, 1, MPI_INT, k, k, MPI_COMM_WORLD, &status);
				MPI_Recv(&(seq.st[t]), count, MPI_CHAR, k, k, MPI_COMM_WORLD, &status);
			}
			t+=count;
		}
		
		/*Cada processo precisa de se situar relativamente à linha e coluna do total da matriz*/
		MPI_Bcast(&ij[0], 2, MPI_INT, k, MPI_COMM_WORLD);
	}
	
	free(ij);
			
	return seq;
}

/*Liberta memória alocada*/
void freeMem(data lcs, stack seq, int flag){
	
	if(flag == 1) free(seq.st);
	
	/*Libertar memória alocada para a matriz de cada máquina*/
	free(lcs.matrix[0]);
	free(lcs.matrix);
	
	/*Libertar memória alocada para os vectores de caractéres de cada máquina*/
	free(lcs.rows);
	free(lcs.columns);
	
	return;
}

/*Imprime o resultado final*/
void print(stack seq){

	printf("%d\n", sizeLCS);
	display(seq);
	
	return;
}

/*MAIN*/
int main(int argc, char *argv[]){
	char *fname=NULL;
	data lcs;
	stack seq;
	int flag=0;
	
	MPI_Init(&argc, &argv);
	
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);

	if(argc != 2){
		if(!id) printf("\n[ERROR] Missing initial arguments\n\n");
		MPI_Finalize();
		exit(-1);		
	}
	
	if(!id) fname = argv[1];

	/*Leitura do ficheiro de input*/
	lcs = readFile(fname);

	/*Computação da Matriz*/
	lcs = computeMatrix(lcs);

	MPI_Barrier(MPI_COMM_WORLD);
	
	/*Obtenção da maior subsequência comum*/
	seq = computeLCS(lcs, seq, &flag);
	
	/*Impressão do resultado*/
	if(!id) print(seq);

	/*Libertação da memória alocada*/
	freeMem(lcs, seq, flag);
	
	MPI_Finalize();
	
	exit(0);
}

