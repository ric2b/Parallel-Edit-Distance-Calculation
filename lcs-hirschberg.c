#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#define max(a,b) (((a) > (b)) ? (a) : (b))

typedef struct t_data{	
	char *rows;
	char *columns;
	int N;
	int M;
}data;

char *createVector(int SIZE){
	char *newV=NULL;
	
	newV = (char*)calloc(SIZE, sizeof(char));
	if(newV == NULL){
		putchar('\n');
		printf("ERROR => Creating newV Vector!!!\n");
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

data readFile(char *fname){
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
  
  newData.rows = createVector(N+1);
  newData.columns = createVector(M+1);
  
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
	
	newData.N = N;
	newData.M = M;
	
	fclose(pFile);
	
	return newData;
}

short cost(int x){
	int i=0, n_iter=20;
	double dcost=0;
	
	for(i=0; i<n_iter; i++){
		dcost += pow(sin((double) x), 2) + pow(cos((double) x), 2);
	}	
	
	return (short) (dcost / n_iter + 0.1);	
}

static int *ALG_B(int m, int n, const char *A, const char *B){
	int i=0,j=0;
	int *k0=NULL;
	int *k1=NULL;
	
	k0 = (int*)calloc(n+1, sizeof(int));
	if(k0 == NULL){
		putchar('\n');
		printf("ERROR => Creating k0 Vector!!!\n");
		putchar('\n');
		exit(-1);
	}
	
	k1 = (int*)calloc(n+1, sizeof(int));
	if(k1 == NULL){
		putchar('\n');
		printf("ERROR => Creating k1 Vector!!!\n");
		putchar('\n');
		exit(-1);
	}
	
	for(i=1; i<=m; i++){
		for(j=0; j<=n; j++)
			k0[j] = k1[j];
		for(j=1; j<=n; j++){
			if (A[i-1] == B[j-1]){
				k1[j] = k0[j-1]+cost(i);   /*cost */
			}else{
				k1[j] = max(k1[j-1], k0[j]);
			}
		}
	}
	free(k0);
	
	return k1;
}


int findK(int *L1, int *L2, int n){
	int W = 0;
	int k = 0;
	int j=0;
	
	for(j=0; j<=n; j++){
		if(W < (L1[j] + L2[n-j])){
			W = L1[j] + L2[n-j];
			k = j;
		}
	}
	
	return k;
}

char *reverseString(char *string, int n){
	char *aux=NULL;
	
	aux = (char*)calloc((n + 1), sizeof(char));
	if(aux == NULL){
		putchar('\n');
		printf("ERROR => Creating aux Vector!!!\n");
		putchar('\n');
		exit(-1);
	}
	
	aux =aux+n;
	*aux = 0;
	while(n--){
		*--aux = *string++;
	}
	
	return aux;
}

char *ALG_C(int m, int n, char *A, char *B) {
	int j=0, i=0, k=0;
	char *C=NULL, *C1=NULL, *C2=NULL, *A1=NULL, *B1=NULL; 
	int *L1=NULL, *L2=NULL;
	
	
	if(n == 0){
		C = (char*)calloc(1,sizeof(char));
		if(C == NULL){
			putchar('\n');
			printf("ERROR => Creating C Vector!!!\n");
			putchar('\n');
			exit(-1);
		}
		return C;
	}else if(m == 1){
		for (j = 0; j < n; j++) {
			if (A[0] == B[j]) {
				C = (char*)calloc(2, sizeof(char));
				if(C == NULL){
					putchar('\n');
					printf("ERROR => Creating C Vector!!!\n");
					putchar('\n');
					exit(-1);
				}
				C[0] = A[0];
				return C;
			}
		}
		C = (char*)calloc(1, sizeof(char));
		if(C == NULL){
			putchar('\n');
			printf("ERROR => Creating C Vector!!!\n");
			putchar('\n');
			exit(-1);
		}
		return C;
	}else{
		i = m / 2;
	
		L1 = ALG_B(i, n, A, B);

		A1 = reverseString(A + i, m - i);
		B1 = reverseString(B, n);

		L2 = ALG_B(m - i, n, A1, B1);

		k = findK(L1, L2, n);

		C1 = ALG_C(i, k, A, B);
		C2 = ALG_C(m - i, n - k, A + i, B + k);

		C = (char*)calloc((strlen(C1) + strlen(C2) + 1), sizeof(char));
		if(C == NULL){
			putchar('\n');
			printf("ERROR => Creating C Vector!!!\n");
			putchar('\n');
			exit(-1);
		}
	
		strcpy(C, C1);
		strcat(C, C2);
	
		/* Libertar Memória */
		free(C1);
		free(C2);
		free(A1);
		free(B1);
		free(L1);
		free(L2);
	
		return C;
	}
}

char *LCS_Hirschberg(char *A, char *B) {
	int n=0, m=0;
	
	n = strlen(A);
	m = strlen(B);
	
	return ALG_C(n, m, A, B);
}

void freeMem(data lcs, char *C){
	
	free(C);
	free(lcs.rows);
	free(lcs.columns);

	return;
}

int main(int argc, char *argv[]){
	char *fname=NULL;
	data lcs;
	double start=0, end=0;
	char *C=NULL;
	
	start = omp_get_wtime();
	
	if(argc != 2){
		putchar('\n');
		printf("ERROR => Arguments!!!\n");
		putchar('\n');
		exit(-1);
	}
	
	fname = argv[1];
	
	lcs = readFile(fname);
	
	C = LCS_Hirschberg(lcs.rows, lcs.columns);
	
	printf("%d\n%s\n", (int)strlen(C), C);
	
	freeMem(lcs, C);
	
	end = omp_get_wtime();
	
	putchar('\n');
	printf("Time: %.5g\n", (end-start));
	
	exit(0);
}
