#ifndef LINEAREQUATIONFILEREADER
#define LINEAREQUATIONFILEREADER
#include <stdlib.h>

double** getMatrixA(int* n /*receberá a quantidade de equações lidas*/)
/*
* Retorna uma matriz A[n][n] com os valores do arquivo randomMatrix.txt
* */
{
    FILE *arq = fopen("../src/randomMatrix.txt", "r"); //abre o arquivo
    fscanf(arq, "%d\n", n); //le a primeira linha, definindo a quantidade de equações
    int i,j; //variáveis para loops

    /*
    * Alocando espaço para matriz A
    * */
    double** A;
    A = (double**) malloc((*n) * sizeof(double*));
    for(i = 0; i < (*n); i++)
    {
        A[i] = (double*) malloc((*n) * sizeof(double));
    }
    
    /*
    * Lendo valores do arquivo
    * */
    for(i = 0; i < (*n); i++) //para cada linha
    {
        for(j = 0; j < (*n); j++) //para cada valor na linha
        {
            fscanf(arq, "%lf ", &A[i][j]); //le um valor e guarda na matriz
        }
        fscanf(arq, "\n"); //le a quebra de linha
    }
    fclose(arq);
    return A; //retorna a matriz
};


double* getMatrixB(int* n /*receberá a quantidade de equações lidas*/)
/*
* Retorna uma matriz A[n][n] com os valores do arquivo randomMatrix.txt
* */
{
    FILE *arq = fopen("../src/randomMatrix.txt", "r"); //abre o arquivo
    fscanf(arq, "%d\n", n); //le a primeira linha, definindo a quantidade de equações
    int i; //variáveis para loops

    /*
    * Alocando espaço para matriz B
    * */
    double* B;
    B = (double*) malloc((*n) * sizeof(double));

    /*
    * Pulando as linhas da matriz A
    * */
    for(i = 0; i < (*n); i++) //para cada linha
    {
        fscanf(arq, "%*[^\n]\n"); //pula uma linha do arquivo
    }
    fscanf(arq, "\n"); //le a quebra de linha

    /*
    * Lendo os campos da matriz B
    * */
    for(i = 0; i < (*n); i++) //para cada linha
    {
        fscanf(arq, "%lf\n", &B[i]); //le a quebra de linha
    }
    fclose(arq);
    return B; //retorna a matriz
};

#endif