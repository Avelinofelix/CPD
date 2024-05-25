#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>


#define ALEATORIO ((double)random() / (double)RAND_MAX)
#define max(x,y) ( (x) > (y) ? (x) : (y) )
#define POS(i,j,columns) ((i)*(columns)+(j))
#define BLOCK_LOW(id,p,n) ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n) - 1)
#define BLOCK_SIZE(id,p,n) (BLOCK_HIGH(id,p,n) - BLOCK_LOW(id,p,n) + 1)
#define BLOCK_OWNER(index,p,n) (((p)*((index)+1)-1)/(n))

int **criar_matriz_compacta(int n);
double **criar_matriz_double(int r, int c);
void liberar_matriz_int(int **matriz, int linhas);
void liberar_matriz_double(double **matriz, int linhas);
void imprimir_matriz_int(int **M, int r, int c);
void preenche_aleatorio_LR(int nU, int nI, int nF);
void ler_entrada(const char *filename);
void criar_estruturas_matrizes();
void liberar_estruturas_matrizes();
void atualizar();
void loop();
void imprimir_recomendacoes();
void imprimir_matriz_double(double **M, int r, int c);
void multiplicar_nao_zeros(double **L, double **RT, double **B, int **A, int nNaoZero, int nCaracteristicas, int nItens, int block_low, int block_size);
void copiar_matriz(double **m1, double **m2, int n, int m);
void multiplicar_matriz(double **X, double **Y, double **Z, int n, int m, int p);
double **transpor_matriz(double **M, int n, int m);

int iteracoes, nCaracteristicas, nUsuarios, nItens, nNaoZero, numThreads, max, block_size, id, nproc, nElements = 0;
int **A, *recomendacoes;
double alpha;
double ***privLsum, ***privRTsum;
double **L, **R, **RT, **B, *Lsum, *RTsum;

void multiplicar_nao_zeros(double **L, double **RT, double **B, int **A, int nNaoZero, int nCaracteristicas, int nItens, int block_low, int block_size) {
    for (int n = 0; n < nNaoZero; n++) {
        int i = A[n][0];
        int j = A[n][1];
        int i2 = i - block_low;

        if (i2 >= 0 && i2 < block_size) {
            double soma = 0;
            for (int k = 0; k < nCaracteristicas; k++) {
                soma += L[i2][k] * RT[j][k];
            }
            B[i2][j] = soma;
        }
    }

     MPI_Allreduce(MPI_IN_PLACE, &B[0][0], nItens * nCaracteristicas, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}
void multiplicar_matriz(double **X, double **Y, double **Z, int n, int m, int p) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            double soma = 0;
            for (int k = 0; k < p; k++)
                soma += X[i][k] * Y[j][k];
            Z[i][j] = soma;
        }
    }
}

double **transpor_matriz(double **M, int n, int m) {
    double **MT = criar_matriz_double(m, n);

    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            MT[j][i] = M[i][j];

    return MT;
}

int **criar_matriz_compacta(int n) {
    int **M = (int **)malloc(n * sizeof(int *));

    for (int i = 0; i < n; i++)
        M[i] = (int *)calloc(3, sizeof(int));
    return M;
}

double **criar_matriz_double(int r, int c) {
    double **M = (double **)malloc(r * sizeof(double *));

    for (int i = 0; i < r; i++)
        M[i] = (double *)calloc(c, sizeof(double));

    return M;
}

void liberar_matriz_int(int **matriz, int linhas) {
    for (int i = 0; i < linhas; i++) {
        free(matriz[i]);
    }
    free(matriz);
}
void liberar_matriz_double(double **matriz, int linhas) {
    for (int i = 0; i < linhas; i++) {
        free(matriz[i]);
    }
    free(matriz);
}

void imprimir_matriz_int(int **M, int r, int c) {
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++)
            printf("%d ", M[i][j]);
        printf("\n");
    }
    printf("\n");
}

void imprimir_matriz_double(double **M, int r, int c) {
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++)
            printf("%f ", M[i][j]);
        printf("\n");
    }
    printf("\n");
}
void preenche_aleatorio_LR(int nU, int nI, int nF) {
    srandom(0);
    for (int i = 0; i < nU; i++) {
        for (int j = 0; j < nF; j++) {
            L[i][j] = ALEATORIO / (double)nF;
        }
    }

    for (int i = 0; i < nF; i++) {
        for (int j = 0; j < nI; j++) {
            R[i][j] = ALEATORIO / (double)nF;
        }
    }
}

void ler_entrada(const char *filename) {
    FILE *ponteiro_arquivo;
    ponteiro_arquivo = fopen(filename, "r");
    if (ponteiro_arquivo == NULL) {
        printf("Erro ao abrir o arquivo %s.\n", filename);
        exit(EXIT_FAILURE);
    }
    fscanf(ponteiro_arquivo, "%d", &iteracoes);
    fscanf(ponteiro_arquivo, "%lf", &alpha);
    fscanf(ponteiro_arquivo, "%d", &nCaracteristicas);
    fscanf(ponteiro_arquivo, "%d", &nUsuarios);
    fscanf(ponteiro_arquivo, "%d", &nItens);
    fscanf(ponteiro_arquivo, "%d", &nNaoZero);

    max = max(nUsuarios, nItens);

    criar_estruturas_matrizes();

    for (int i = 0; i < nNaoZero; i++) {
        int n, m;
        double v;

        fscanf(ponteiro_arquivo, "%d", &n);
        fscanf(ponteiro_arquivo, "%d", &m);
        fscanf(ponteiro_arquivo, "%lf", &v);

        A[i][0] = n;
        A[i][1] = m;
        A[i][2] = v;
    }

    fclose(ponteiro_arquivo);
}

void criar_estruturas_matrizes() {
    A = criar_matriz_compacta(nNaoZero);
    B = criar_matriz_double(nUsuarios, nItens);
    L = criar_matriz_double(nUsuarios, nCaracteristicas);
    R = criar_matriz_double(nCaracteristicas, nItens);
    RT = criar_matriz_double(nItens, nCaracteristicas);
    Lsum = (double *)calloc(nUsuarios * nCaracteristicas, sizeof(double)); 
    RTsum = (double *)calloc(nItens * nCaracteristicas, sizeof(double)); 
}
void liberar_estruturas_matrizes() {
    liberar_matriz_int(A, nNaoZero);
    liberar_matriz_double(B, nUsuarios);
    liberar_matriz_double(L, nUsuarios);
    liberar_matriz_double(R, nCaracteristicas);
    liberar_matriz_double(RT, nItens);

}
void atualizar() {
    int i, j, n, k;

 
    for (n = 0; n < nNaoZero; n++) {
        i = A[n][0];
        j = A[n][1];

        for (k = 0; k < nCaracteristicas; k++) {
            Lsum[i * nCaracteristicas + k] += alpha * (2 * (A[n][2] - B[i][j]) * (-RT[j][k]));
            RTsum[j * nCaracteristicas + k] += alpha * (2 * (A[n][2] - B[i][j]) * (-L[i][k]));
        }
    }

      MPI_Allreduce(MPI_IN_PLACE, Lsum, nUsuarios * nCaracteristicas, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, RTsum, nItens * nCaracteristicas, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    
    for (i = 0; i < nUsuarios; i++) {
        for (k = 0; k < nCaracteristicas; k++) {
            L[i][k] -= Lsum[i * nCaracteristicas + k];
        }
    }

    for (j = 0; j < nItens; j++) {
        for (k = 0; k < nCaracteristicas; k++) {
            RT[j][k] -= RTsum[j * nCaracteristicas + k];
        }
    }
}
void copiar_matriz(double **origem, double **destino, int linhas, int colunas) {
    for (int i = 0; i < linhas; i++) {
        for (int j = 0; j < colunas; j++) {
            destino[i][j] = origem[i][j];
        }
    }
}

void imprimir_recomendacoes() {
    int indice = 0;
    for (int usuario = 0; usuario < nUsuarios; usuario++) {
        double maximo = -1;
        int recomendacao;

        for (int item = 0; item < nItens; item++) {
            if (indice < nNaoZero && A[indice][0] == usuario && A[indice][1] == item) {
                indice++;
                continue;
            }

            if (B[usuario][item] > maximo) {
                maximo = B[usuario][item];
                recomendacao = item;
            }
        }
        printf("%d\n", recomendacao);
    }
}

void loop() {
    int block_low = BLOCK_LOW(id, nproc, nUsuarios);
    int block_size = BLOCK_SIZE(id, nproc, nUsuarios);

    for (int i = 0; i < iteracoes; i++) {
        multiplicar_nao_zeros(L, RT, B, A, nNaoZero, nCaracteristicas, nItens, block_low, block_size);
        atualizar();
    }

     MPI_Allreduce(MPI_IN_PLACE, &B[0][0], nUsuarios * nItens, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}
int main(int argc, char **argv) {
    double start, end;

    if (argc != 2) {
        printf("Número inválido de argumentos fornecidos.\n");
        printf("Uso: %s <nome_do_arquivo_de_entrada>\n", argv[0]);
        return -1;
    }

    MPI_Init(&argc, &argv);
    MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    if (id == 0) {
        ler_entrada(argv[1]);
        preenche_aleatorio_LR(nUsuarios, nItens, nCaracteristicas);
        RT = transpor_matriz(R, nCaracteristicas, nItens);
    }

     MPI_Bcast(&iteracoes, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&alpha, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nCaracteristicas, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nUsuarios, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nItens, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nNaoZero, 1, MPI_INT, 0, MPI_COMM_WORLD);

     criar_estruturas_matrizes();

    MPI_Bcast(&A[0][0], nNaoZero * 3, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&L[0][0], nUsuarios * nCaracteristicas, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&RT[0][0], nItens * nCaracteristicas, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    int block_low = BLOCK_LOW(id, nproc, nUsuarios);
    int block_size = BLOCK_SIZE(id, nproc, nUsuarios);

    
    for (int i = 0; i < iteracoes; i++) {
         multiplicar_nao_zeros(L, RT, B, A, nNaoZero, nCaracteristicas, nItens, block_low, block_size);

        atualizar();

         MPI_Allreduce(MPI_IN_PLACE, &B[0][0], nUsuarios * nItens, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }

     if (id == 0) {
        imprimir_recomendacoes();
    }

   liberar_estruturas_matrizes();
    MPI_Finalize();

    return 0;
}