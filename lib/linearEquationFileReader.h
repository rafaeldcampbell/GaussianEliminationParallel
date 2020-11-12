#ifndef LINEAREQUATIONFILEREADER
#define LINEAREQUATIONFILEREADER
#include <stdlib.h>

int getMatrixAandB(int* n /*receberá a quantidade de equações lidas*/,
               double*** A /*receberá a matriz A*/,
               double** B /*receberá a matriz B*/)
/*
* Recebe os endereços de n (&(int*)), A (&(double**)) e B(&(double*)).
* Define os valores de n, A[n][n] e B[n] com as matrizes do arquivo randomMatrix.txt
* */
{
    FILE *arq = fopen("../src/randomMatrix.txt", "r"); //abre o arquivo
    if(arq == NULL) {
        fprintf(stderr, "Erro ao abrir o randomMatrix.txt.");
        return -1;
    }

    fscanf(arq, "%d\n", n); //le a primeira linha, definindo a quantidade de equações
    int i,j; //variáveis para loops

    /*
    * Alocando espaço para matriz A
    * */
    *A = (double**) malloc((*n) * sizeof(double*));
    for(i = 0; i < (*n); i++)
    {
        (*A)[i] = (double*) malloc((*n) * sizeof(double));
    }

    /*
    * Lendo valores do arquivo para matriz A
    * */
    for(i = 0; i < (*n); i++) //para cada linha
    {
        for(j = 0; j < (*n); j++) //para cada valor na linha
        {
            fscanf(arq, "%lf ", &((*A)[i][j])); //le um valor e guarda na matriz
        }
        fscanf(arq, "\n"); //le a quebra de linha
    }
    fscanf(arq, "\n"); //le a quebra de linha

     /*
    * Alocando espaço para matriz B
    * */
    *B = (double*) malloc((*n) * sizeof(double));
    
    /*
    * Lendo os campos da matriz B
    * */
    for(i = 0; i < (*n); i++) //para cada linha
    {
        fscanf(arq, "%lf\n", &((*B)[i])); //le a quebra de linha
    }

    fclose(arq);
    return 0; //retorna a matriz
};

#endif