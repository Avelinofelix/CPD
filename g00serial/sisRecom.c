
/*      Membros do grupo
1-Avelino Félix Adelino-20210568
Obs: Estou sozinho Sr Professor



   */

#include <stdio.h>
#include <stdlib.h>



#define ALEATORIO ((double)random() / (double)RAND_MAX)
#define max(x,y) ( (x) > (y) ? (x) : (y) )



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
void multiplicar_nao_zeros(double **L, double **RT, double **B, int **A, int nNaoZero, int nCaracteristicas);
void copiar_matriz(double **m1, double **m2, int n, int m);
void multiplicar_matriz(double **X, double **Y, double **Z, int n, int m, int p);
double **transpor_matriz(double **M, int n, int m);




int iteracoes, nCaracteristicas, nUsuarios, nItens, nNaoZero, numThreads, max;
int **A;
double alpha;
double ***privLsum, ***privRTsum;
double **L, **R, **RT, **B;




int main(int argc, char **argv) {
    if(argc != 2) {
        printf("Número inválido de argumentos fornecidos.\n");
        printf("Uso: %s <nome_do_arquivo_de_entrada>\n", argv[0]);
        return -1;
    }

ler_entrada(argv[1]);

preenche_aleatorio_LR(nUsuarios, nItens, nCaracteristicas);

RT = transpor_matriz(R, nCaracteristicas, nItens);

loop();

imprimir_recomendacoes();

liberar_estruturas_matrizes();

return 0;

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

    for (int i= 0; i < n; i++) 
        M[i] = (int *) calloc(3, sizeof(int));
    return M;
}

double **criar_matriz_double(int r, int c) {
    double **M = (double **)malloc(r * sizeof(double *)); 

    for (int i= 0; i < r; i++) 
        M[i] = (double *) calloc(c, sizeof(double));

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
void multiplicar_nao_zeros(double **L, double **RT, double **B, int **A, int nNaoZero, int nCaracteristicas) {
    for (int n = 0; n < nNaoZero; n++) {
        int i = A[n][0];
        int j = A[n][1];

        double soma = 0;
        for (int  k = 0; k < nCaracteristicas; k++) {
            soma += L[i][k] * RT[j][k];
        }
        B[i][j] = soma;
    }
}

void imprimir_matriz_double(double **M, int r, int c) {
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++)
            printf("%f ", M[i][j]);
        printf("\n");
    }
    printf("\n");
}


void ler_entrada(const char *filename) {
    FILE *ponteiro_arquivo;
    ponteiro_arquivo = fopen(filename, "r");
    if(ponteiro_arquivo == NULL) {
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
}

void liberar_estruturas_matrizes() {
    liberar_matriz_int(A, nNaoZero);
    liberar_matriz_double(B, nUsuarios);
    liberar_matriz_double(L, nUsuarios);
    liberar_matriz_double(R, nCaracteristicas);
    liberar_matriz_double(RT, nItens);
}

void preenche_aleatorio_LR(int nU, int nI, int nF) {
    srandom(0);
    int i, j;
    for (i = 0; i < nU; i++)
        for (j = 0; j < nF; j++)
            L[i][j] = ALEATORIO / (double) nF;

    for (i = 0; i < nF; i++)
        for (j = 0; j < nI; j++)
            R[i][j] = ALEATORIO / (double) nF;
}

double **Lnew;

double **Lnew;
double **RTnew;

void atualizar() {
    int i, j, n, k;

    // Certifique-se de que Lnew e RTnew estão alocados
    Lnew = criar_matriz_double(nUsuarios, nCaracteristicas);
    RTnew = criar_matriz_double(nItens, nCaracteristicas);
    
    // Copie L e RT para Lnew e RTnew, respectivamente
    copiar_matriz(L, Lnew, nUsuarios, nCaracteristicas);
    copiar_matriz(RT, RTnew, nItens, nCaracteristicas);

    // Atualize Lnew e RTnew
    for (n = 0; n < nNaoZero; n++) {
        i = A[n][0];
        j = A[n][1];

        for (k = 0; k < nCaracteristicas; k++) {
            Lnew[i][k] -= alpha * (2 * (A[n][2] - B[i][j]) * (-RT[j][k]));
            RTnew[j][k] -= alpha * (2 * (A[n][2] - B[i][j]) * (-L[i][k]));
        }
    }

    // Troque os ponteiros para atualizar L e RT
    double **mAux;
    mAux = Lnew;    Lnew = L;       L = mAux;
    mAux = RTnew;   RTnew = RT;     RT = mAux;

    // Libere a memória alocada para Lnew e RTnew
    liberar_matriz_double(Lnew, nUsuarios);
    liberar_matriz_double(RTnew, nItens);
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
    for (int i = 0; i < iteracoes; i++) {
        multiplicar_nao_zeros(L, RT, B, A, nNaoZero, nCaracteristicas);
        atualizar();
    }
    multiplicar_matriz(L, RT, B, nUsuarios, nItens, nCaracteristicas);
}
