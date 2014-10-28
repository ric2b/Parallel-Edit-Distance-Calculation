#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#define max(a,b) (((a) > (b)) ? (a) : (b))

typedef struct t_stack{
	char *st;
	int top;
}stack;

stack createStack(int sizeStack){
	stack new;

	new.st = (char*)calloc(sizeStack, sizeof(char));

	if(new.st == NULL){
		putchar('\n');
		printf("ERROR => Creating Stack!!!\n");
		putchar('\n');
		exit(-1);
	}
	
	new.top = -1;
	
	return new;
}

stack push(stack seq, char item){
  seq.top++;
  seq.st[seq.top] = item;
  
	return seq;
}

int stempty(stack seq){
  if(seq.top == -1) return 1;
  else return 0;
}

void display(stack seq){
  int i=0;
   
	if(stempty(seq)){ 
  	printf("No Longest Common Subsequence Found!\n");
  }else{
 		for(i=seq.top; i>=0; i--){
    	printf("%c", seq.st[i]);
  	}
  	putchar('\n');
  }
  
  return;
}

typedef struct t_data{
	int **matrix;
	char *rows;
	char *columns;
	int N;
	int M;
}data;

data lcs;

int **createMatrix(int N, int M){
	int i=0;
	int **newM=NULL;
	
	newM = (int**)calloc((N+1), sizeof(int*));
	if(newM == NULL){
		putchar('\n');
		printf("ERROR => Creating Matrix!!!\n");
		putchar('\n');
		exit(-1);
	}
	
	for(i=0; i<(N+1); i++){
		newM[i] = (int*)calloc((M+1), sizeof(int));
		if(newM[i] == NULL){
			putchar('\n');
			printf("ERROR => Creating Matrix!!!\n");
			putchar('\n');
			exit(-1);
		}
	}
	
	return newM;
}

char *createVector(int SIZE){
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

FILE *openFile(char *fname){
	FILE *fp=NULL;
	
	
	fp = fopen(fname, "r");
	if(fp == NULL){
		putchar('\n');
		printf("ERROR => Opening File!!!\n");
		putchar('\n');
		exit(-1);
	}
	
	return fp;
}

void readFile(char *fname){
	int N=0, M=0;
	FILE *pFile=NULL;
	data newData;
	
	pFile = openFile(fname);
	if(fscanf(pFile, "%d %d", &N, &M) != 2){
    printf("ERROR => Reading variables from file!!!\n");
    fclose(pFile);
    exit(-1); 
  }
  if(fscanf(pFile, "\n") != 0){
  	printf("ERROR => Reading EOF from file!!!\n");
    fclose(pFile);
    exit(-1);
  }
  
  newData.rows = createVector(N);
  newData.columns = createVector(M);
  
  if(fread(newData.rows, 1, N, pFile) != N){
    printf("ERROR => Reading information from file!!!\n");
    fclose(pFile);
    free(newData.rows);
    free(newData.columns);
    exit(-1); 
  }
  if(fscanf(pFile, "\n") != 0){
  	printf("ERROR => Reading EOF from file!!!\n");
    free(newData.rows);
    free(newData.columns);
    fclose(pFile);
    exit(-1);
  }
	
	if(fread(newData.columns, 1, M, pFile) != M){
    printf("ERROR => Reading information from file!!!\n");
    fclose(pFile);
    free(newData.rows);
    free(newData.columns);
    exit(-1); 
  }
	
	newData.matrix = createMatrix(N,M);
	newData.N = N;
	newData.M = M;
	
	fclose(pFile);
	
	lcs = newData;
	
	return;
}

short cost(int x){
	int i=0, n_iter=20;
	double dcost=0;
	
	for(i=0; i<n_iter; i++){
		dcost += pow(sin((double) x), 2) + pow(cos((double) x), 2);
	}	
	
	return (short) (dcost / n_iter + 0.1);	
}

void ComputeCell(int k, int t){
	if(k==0 || t==0){
		lcs.matrix[k][t] = 0;		
	}else if(lcs.rows[k-1] == lcs.columns[t-1]){
		lcs.matrix[k][t] = lcs.matrix[k-1][t-1] + cost(k);
	}else{
		lcs.matrix[k][t] = max(lcs.matrix[k][t-1], lcs.matrix[k-1][t]);
	}
}


void computematrix(int i, int j, int N, int M, int n, int m, int caso){
	int f=0,g=0,ver=0;
		
	switch(caso){
		case 1:
			for(f=i-n/2; f<n ;f++){
				for(g=j-n/2; g<m ;g++){
					ComputeCell(f, g);
				}
			}
			if(f==N && g==M){
				return;
			}else if((2*n)<=N && (2*n)<=M){
				computematrix(f,0,N,M,2*n,n,2);
			}else{
					for(; f<=N; f++){
						for(g=1 ; g<n; g++){
							ComputeCell(f, g);
						}
					}
					ver++;
					for(; g<=M; g++){
						for(f=1; f<n; f++){
							ComputeCell(f, g);
						}
					}
					ver++;
					for(f=n ;f<=N; f++){
						for(g=n; g<=M; g++){
							ComputeCell(f, g);
						}	
					}		
					return;
			}			
			break;
	
		case 2:
			for(f=i; f<n; f++){
				for(g=j; g<m; g++){
					ComputeCell(f, g);
				}
			}			
			for(g=i; g<n; g++){
				for(f=j; f<m ; f++){
					ComputeCell(f, g);
				}
			}
			computematrix(g,g, N,M, n,n,1);
			break;
			
		default:
			break;
	}
	
	return;
}

stack computeLCS(stack seq){
	int i, j, sizeStack=0;
	
	i = lcs.N;
	j = lcs.M;
	
	sizeStack = lcs.matrix[i][j];
	seq = createStack(sizeStack);
	
	while(i != 0 && j != 0){
		if(lcs.rows[i-1] == lcs.columns[j-1]){
			seq = push(seq, lcs.rows[i-1]);
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

	return seq;
}

void print(stack seq){
	
	printf("%d\n", lcs.matrix[lcs.N][lcs.M]);
	display(seq);
	
	return;
}

void freeMem(stack seq){
	int i=0;

	free(seq.st);
	
	for(i=0; i<lcs.N+1; i++){
		free(lcs.matrix[i]);
	}
	free(lcs.matrix);
	free(lcs.rows);
	free(lcs.columns);

	return;
}

int main(int argc, char *argv[]){
	char *fname=NULL;
	stack seq;
	double start=0, end=0;
	/* int i=0, j=0; */
	
	start = omp_get_wtime();
	
	if(argc != 2){
		putchar('\n');
		printf("ERROR => Arguments!!!\n");
		putchar('\n');
		exit(-1);
	}
	
	fname = argv[1];
	
	readFile(fname);
	
	computematrix(0,0, lcs.N, lcs.M, 1,1,1);
	
	seq = computeLCS(seq);
	
	print(seq);
	
	/* Impressão da Matriz */
	/* for(i=0; i<lcs.N+1; i++){
		for(j=0; j<lcs.M+1; j++){
			printf("%d ", lcs.matrix[i][j]);
		}
		putchar('\n');
	} */
	
	freeMem(seq);
	
	end = omp_get_wtime();
	
	putchar('\n');
	printf("Time: %.5g\n", (end-start));
	
	exit(0);
}
