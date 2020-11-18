#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <string.h>
#include "../lib/linearEquationFileReader.h"

#define NTHREADS 4

typedef struct Compare { double val; int index; } compare;  //estrutura de apoio para pivoteamento
//declarando um reduction customizado para pivoteamento
#pragma omp declare reduction(maximum : compare : omp_out = omp_in.val > omp_out.val ? omp_in : omp_out)

/*
* Soma todas as colunas já calculadas -- multiplicando pelas variáveis obtidas -- e subtrai
* do total (B). O valor resultante, dividido pelo coeficiente não calculado, representa o
* valor de X daquela coluna.
* */
void backSubstitution(int m, double **a, double *b, double* k)
{
    int i = 0, j = 0;
    double sum;
    for(i=m-1; i >= 0; i--) //para cada linha, ao contrario
    {
        sum = 0;
        for(j=i+1; j < m; j++) //para cada coluna
        {
            sum += a[i][j] * k[j];
        }
        k[i] = (b[i] - sum) / a[i][i];
    }
}

void gaussElimination(int m, double **a, double *b){
    // A[linhas][colunas]
    // A -- matriz de coeficientes
    // B -- vetor de solução
    int i = 0, j = 0, k = 0;
    for(i=0;i<m-1;i++){ //para cada coluna

        //Pivoteamento parcial
        compare max; //guarda o valor absoluto e o index respectivo
        max.val = fabs(a[i][i]);
        max.index = i;
        #pragma omp parallel for private(k) reduction(maximum:max) //faz uma redução buscando o index do maior valor
        for(k=i+1;k<m;k++){ // para cada linha
            //Faz a troca das equações caso encontre um valor absoluto maior
            if(max.val<fabs(a[k][i])){
                max.val = fabs(a[k][i]);
                max.index = k;
            }
        }

        //troca as linhas, caso o pivot atual não seja o mais expressivo
        if(max.index != i)
        {
            for(j=0;j<m;j++){  //troca os valores de A          
                double temp;
                temp=a[i][j];
                a[i][j]=a[max.index][j];
                a[max.index][j]=temp;
            }
            double temp = b[max.index]; //troca o valor de B
            b[max.index] = b[i];
            b[i] = temp;
        }
        //Begin Gauss Elimination
        double term;
        #pragma omp parallel for shared(a,b) private(term, k,j)
        for(k=i+1;k<m;k++){ //para cada linha
            term = a[k][i]/ a[i][i]; //calcula o termo pelo qual a linha pivo sera multiplicada
            for(j=0;j<m;j++){ //para cada coluna
                a[k][j] = a[k][j] - (term * a[i][j]); //faz a subtituição da linha por uma equivalente
            }
            b[k] = b[k] - (term*b[i]); //substitui o valor de B correspondente
        }
    }         
}

int main(int argc, char **argv)
{
    omp_set_num_threads(NTHREADS); //define numero de threads

    int printMatrix = 0;
    int printResult = 0;
    if(argc > 1 && !strcmp(argv[1], "-p")) //parâmetro "-p" determina a impressão das matrizes
        printMatrix = 1;
    if(argc > 1 && !strcmp(argv[1], "-r")) //parâmetro "-r" determina a impressão apenas do resultado
        printResult = 1;

    int m;
    double **a;
    double *b;
    getMatrixAandB(&m, &a, &b); //realiza a leitura do arquivo gerado
    double *x = malloc(m*sizeof(double)); //aloca espaço pro vetor de variáveis

    double time = -omp_get_wtime(); //inicia a contagem de tempogc

    if(printMatrix){ //imprime a matriz A extendida
        printf("*** Matriz A|B****\n");
        for(int i=0; i < m; i++)
        {
            for(int j=0; j < m; j++)
                printf("%lf ", a[i][j]);
            printf("%lf\n", b[i]);
        }
        printf("\n");
    }

    gaussElimination(m, a, b); //realiza a eliminação gaussiana
    backSubstitution(m, a, b, x); //realiza a substituição regressiva

    if(printMatrix){ //imprime a matriz A extendida
        printf("*** Matriz Triangulada ****\n");
        for(int i=0; i < m; i++)
        {
            for(int j=0; j < m; j++)
                printf("%lf ", a[i][j]);
            printf("%lf\n", b[i]);
        }
        printf("\n");
    }

    if(printMatrix || printResult){ //imprime o vetor X calculado
        printf("*** Vetor X ***\n");
        for(int j=0; j < m; j++)
        {
            printf("X%d = %lf\n", j, x[j]);
        }
    }

    /*
    * Desalocando todos os vetores utilizados
    * */
    for(int j=0; j < m; j++){
        free(a[j]);
    }
    free(a);
    free(b);
    free(x);

    time += omp_get_wtime();
    printf("**** Cálculo realizado em %lf segundos\n", time);


    return 0;

}