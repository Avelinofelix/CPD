#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

// Macros
#define RAND01 ((double) random() / (double) RAND_MAX)
#define max(x,y) ( (x) > (y) ? (x) : (y) )
#define POS(i,j,columns) ((i)*(columns)+(j))
#define BLOCK_LOW(id,p,n) ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n) - 1)
#define BLOCK_SIZE(id,p,n) \
 (BLOCK_HIGH(id,p,n) - BLOCK_LOW(id,p,n) + 1)
#define BLOCK_OWNER(index,p,n) \
(((p)*((index)+1)-1)/(n))

// Declarações de funções
void preenche_aleatorio_LR(int nU, int nI, int nF);
void ler_entrada(FILE *file_pointer);
void criar_estruturas_matriz();
void liberar_estruturas_matriz();
void atualizar();
void loop();
void imprimir_recomendacoes();

void multiplicar_nao_zeros(double *L, double *RT, double *B, int *A, int nNaoZero, int nCaracteristicas, int nItens, int block_low);
void copiar_matriz(double *m1, double *m2, int n, int m);
void multiplicar_matriz(double *X, double *Y, double *Z, int n, int m, int p);
double *transpor_matriz(double *M, int n, int m);

int *criar_matriz_compacta(int n);
double *criar_matriz_double(int r, int c);

void liberar_matriz_int(int *matriz);
void liberar_matriz_double(double *matriz);

void imprimir_matriz_int(int *M, int r, int c);
void imprimir_matriz_double(double *M, int r, int c);

// ---------------------
// OPERAÇÕES DE MATRIZES
// ---------------------

void multiplicar_nao_zeros(double *L, double *RT, double *B, int *A, int nElementos, int nCaracteristicas, int nItens, int block_low) {
    for (int n = 0; n < nElementos; n++) {
        int i = A[POS(n,0,3)];
        int j = A[POS(n,1,3)];
        int i2 = i - block_low;

        double sum = 0;
        for (int  k= 0; k < nCaracteristicas; k++) {
            sum += L[POS(i2,k,nCaracteristicas)] * RT[POS(j,k,nCaracteristicas)];
        }
        B[POS(i2,j,nItens)] = sum;
    }
}

void multiplicar_matriz(double *X, double *Y, double *Z, int n, int m, int p) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            double sum = 0;
            for (int k = 0; k < p; k++)
                sum += X[POS(i,k,p)] * Y[POS(j,k,p)];
            Z[POS(i,j,m)] = sum;
        }
    }
}

double *transpor_matriz(double *M, int n, int m) {
    double *MT = criar_matriz_double(m, n);
    
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++) {
            MT[POS(j,i,n)] = M[POS(i,j,m)];
        }

    return MT;
}

// -------------------
// CRIAÇÃO DE MATRIZES
// -------------------

int *criar_matriz_compacta(int n) {
    int *M = (int *)calloc(3 * n, sizeof(int));
    return M;
}

double *criar_matriz_double(int r, int c) {
    double *M = (double *)calloc(c * r, sizeof(double)); 
    return M;
}

void liberar_matriz_int(int *matriz) {
    free(matriz);
}

void liberar_matriz_double(double *matriz) {
    free(matriz);
}

// -------------------
// IMPRESSÃO DE MATRIZES
// -------------------

void imprimir_matriz_int(int *M, int r, int c) {
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++)
            printf("%d ", M[POS(i,j,c)]);
        printf("\n");
    }
    printf("\n");
}

void imprimir_matriz_double(double *M, int r, int c) {
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++)
            printf("%f ", M[POS(i,j,c)]);
        printf("\n");
    }
    printf("\n");
}

// Variáveis Globais
int iteracoes, nCaracteristicas, nUsuarios, nItens, nNaoZero, numThreads, max, block_size, id, nproc, nElementos = 0;
int *A, *recomendacoes;
double *L, *R, *RT, *B, *Lsum, *RTsum;
double alpha;

// Principal
int main(int argc, char **argv) {
    
    double start, end;

    if(argc != 2) {
        printf("Número inválido de argumentos fornecidos");
        return -1;
    }

    MPI_Init(&argc, &argv);
    MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    
    FILE *file_pointer;
    file_pointer = fopen(argv[1], "r");
    fscanf(file_pointer, "%d", &iteracoes);
    fscanf(file_pointer, "%lf", &alpha);
    fscanf(file_pointer, "%d", &nCaracteristicas);
    fscanf(file_pointer, "%d", &nUsuarios);
    fscanf(file_pointer, "%d", &nItens);
    fscanf(file_pointer, "%d", &nNaoZero);
    if (id) fclose(file_pointer);

    block_size = BLOCK_SIZE(id, nproc, nUsuarios);
    
    criar_estruturas_matriz();

    if(!id) {
        
        // Ler toda a entrada e criar as estruturas necessárias
        ler_entrada(file_pointer);

        // Preencher as matrizes aleatoriamente
        preenche_aleatorio_LR(nUsuarios, nItens, nCaracteristicas);

        // Obter melhor desempenho por causa do cache 
        // trabalhando nos mesmos blocos de memória -> acertos de cache 
        RT = transpor_matriz(R, nCaracteristicas, nItens);

    } else {
        MPI_Recv(&nElementos, 1, MPI_INT, 0, id, MPI_COMM_WORLD, &status);

        A = criar_matriz_compacta(nElementos);

        MPI_Recv(A, nElementos * 3, MPI_INT, 0, id, MPI_COMM_WORLD, &status);

        // Receber L
        MPI_Recv(L, BLOCK_SIZE(id, nproc, nUsuarios) * nCaracteristicas, MPI_DOUBLE, 0, id, MPI_COMM_WORLD, &status);
    }

    // Broadcast RT
    MPI_Bcast(RT, nItens * nCaracteristicas, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    loop();

    imprimir_recomendacoes();

    liberar_estruturas_matriz();

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();
    
    return 0;
}

void ler_entrada(FILE *file_pointer) {

    int *buffer = criar_matriz_compacta(((block_size + 2) * nItens) * 3);

    int proc_number = 0;
    int elem_number = 0;
    int row_number = 0;
    int prev_n = 0;

    for (int i = 0; i < nNaoZero; i++) {
        int n, m;
        double v;

        int high = BLOCK_HIGH(proc_number, nproc, nUsuarios);

        fscanf(file_pointer, "%d", &n);
        fscanf(file_pointer, "%d", &m);
        fscanf(file_pointer, "%lf", &v);

        if (n <= high) {
            if (proc_number == 0) {
                A[POS(i,0,3)] = n;
                A[POS(i,1,3)] = m;
                A[POS(i,2,3)] = v;
                nElementos++;
            } else {
                buffer[POS(elem_number,0,3)] = n;
                buffer[POS(elem_number,1,3)] = m;
                buffer[POS(elem_number,2,3)] = v;
                elem_number++;
            }

            if (i == nNaoZero - 1) {
                MPI_Send(&elem_number, 1, MPI_INT, proc_number, proc_number, MPI_COMM_WORLD);
                MPI_Send(buffer, elem_number * 3, MPI_INT, proc_number, proc_number, MPI_COMM_WORLD);
            }
        } else {
            if (proc_number) {
                MPI_Send(&elem_number, 1, MPI_INT, proc_number, proc_number, MPI_COMM_WORLD);
                MPI_Send(buffer, elem_number * 3, MPI_INT, proc_number, proc_number, MPI_COMM_WORLD);
                elem_number = 0;
            }
            buffer[POS(elem_number,0,3)] = n;
            buffer[POS(elem_number,1,3)] = m;
            buffer[POS(elem_number,2,3)] = v;
            elem_number++;
            proc_number++;
        }
    }
    
    liberar_matriz_int(buffer);
    fclose(file_pointer);
}

void criar_estruturas_matriz() {
    if(id) {
        RT = criar_matriz_double(nItens, nCaracteristicas);
    } else {
        A = criar_matriz_compacta(nNaoZero);
        R = criar_matriz_double(nCaracteristicas, nItens);
    }

    L = criar_matriz_double(block_size, nCaracteristicas);
    Lsum = criar_matriz_double(block_size, nCaracteristicas);
    RTsum = criar_matriz_double(nItens, nCaracteristicas);
    B = criar_matriz_double(block_size, nItens); 

}

void liberar_estruturas_matriz() {
    free(A);
    free(B);
    free(L);
    free(Lsum);
    free(RT);
    free(RTsum);
    if (!id) {
        free(R);
    }
    free(recomendacoes);
}

void preenche_aleatorio_LR(int nU, int nI, int nF) {
    srandom(0);

    double *buffer = criar_matriz_double(block_size + 1, nF);

    for (int i = 0; i < block_size; i++)
        for (int j = 0; j < nF; j++)
            L[POS(i,j,nF)] = RAND01 / (double) nF;

    for (int i = block_size, p = 1; i < nU; i++) {
        for (int j = 0; j < nF; j++)
            buffer[POS(i - BLOCK_LOW(p, nproc, nU),j,nF)] = RAND01 / (double) nF; 

        if (BLOCK_OWNER(i + 1,nproc,nU) > p) {
                MPI_Send(buffer, BLOCK_SIZE(p, nproc, nU) * nF, MPI_DOUBLE, p, p, MPI_COMM_WORLD);
                p++;
        }
    }

    for (int i = 0; i < nF; i++)
        for (int j = 0; j < nI; j++)
            R[POS(i,j,nI)] = RAND01 / (double) nF;
}

void atualizar() {
    int i, j, n, k, i2;

    for (n = 0; n < nElementos; n++) {
        i = A[POS(n,0,3)];
        j = A[POS(n,1,3)];
        i2 = i - BLOCK_LOW(id, nproc, nUsuarios);

        for (k = 0; k < nCaracteristicas; k++) {
            Lsum[POS(i2,k,nCaracteristicas)] += alpha * ( 2 * ( A[POS(n,2,3)] - B[POS(i2,j,nItens)] ) * ( -RT[POS(j,k,nCaracteristicas)] ) );
            RTsum[POS(j,k,nCaracteristicas)] += alpha * ( 2 * ( A[POS(n,2,3)] - B[POS(i2,j,nItens)] ) * ( -L[POS(i2,k,nCaracteristicas)] ) );
        }
    }

    MPI_Allreduce(MPI_IN_PLACE, RTsum, nCaracteristicas * nItens, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    for (int i = 0; i < max(block_size, nItens); i++)
        for (int j = 0; j < nCaracteristicas; j++) {
            if (i < block_size) {
                L[POS(i,j,nCaracteristicas)] -= Lsum[POS(i,j,nCaracteristicas)];
                Lsum[POS(i,j,nCaracteristicas)] = 0;
            }
            if (i < nItens) {
                RT[POS(i,j,nCaracteristicas)] -= RTsum[POS(i,j,nCaracteristicas)];
                RTsum[POS(i,j,nCaracteristicas)] = 0;
            }
        }
}

void loop() {
    
    for (int i = 0; i < iteracoes; i++) {
        multiplicar_nao_zeros(L, RT, B, A, nElementos, nCaracteristicas, nItens, BLOCK_LOW(id, nproc, nUsuarios));
        atualizar();
    }

    multiplicar_matriz(L, RT, B, block_size, nItens, nCaracteristicas);
}

void imprimir_recomendacoes() {
    recomendacoes = (int *) malloc((block_size + 1) * sizeof(int));

    for (int user = 0, i = 0, index = 0; user < block_size; user++) {
        double max = -1;
        int recomendacao;

        for (int item = 0; item < nItens; item++) {
            if (index < nElementos && ( A[POS(index,0,3)] - BLOCK_LOW(id, nproc, nUsuarios) ) == user && A[POS(index,1,3)] == item) {
                index++;
                continue;
            }

            if (B[POS(user,item,nItens)] > max) {
                max = B[POS(user,item,nItens)];
                recomendacao = item;
            }
        }
        recomendacoes[i++] = recomendacao;
    }

    if (id) {   
        MPI_Send(recomendacoes, block_size, MPI_INT, 0, id, MPI_COMM_WORLD);
    } else {
        for (int i = 0; i < block_size; i++)
                printf("%d\n", recomendacoes[i]);

        MPI_Status status;
        for (int i = 1; i < nproc; i++) {
            int proc_block_size = BLOCK_SIZE(i, nproc, nUsuarios);
            // Processo 0 recebe recomendações
            MPI_Recv(recomendacoes, BLOCK_SIZE(i, nproc, nUsuarios), MPI_INT, i, i, MPI_COMM_WORLD, &status);
            for (int i = 0; i < proc_block_size; i++)
                printf("%d\n", recomendacoes[i]);
        }
    }
}
