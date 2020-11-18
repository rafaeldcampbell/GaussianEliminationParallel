#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "../lib/linearEquationFileReader.h"

/*
* Soma todas as colunas já calculadas -- multiplicando pelas variáveis obtidas -- e subtrai
* do total (B). O valor resultante, dividido pelo coeficiente não calculado, representa o
* valor de X daquela coluna.
* */
void backSubstitution(int m, double **a, double *b, double* x)
{
    int i = 0, j = 0;
    double sum;
    for(i=m-1; i >= 0; i--) //para cada linha, ao contrario
    {
        sum = 0;
        for(j=i+1; j < m; j++) //para cada coluna
        {
            sum += a[i][j] * x[j];
        }
        x[i] = (b[i] - sum) / a[i][i];
    }
}

void gaussElimination(int m, double **a, double *b){
    // A[linhas][colunas]
    // A -- matriz de coeficientes
    // B -- vetor de solução
    int i = 0, j = 0, k = 0;
    for(i=0;i<m-1;i++){ //para cada coluna
        //Pivoteamento parcial

        for(k=i+1;k<m;k++){ // para cada linha
            //Faz a troca das equações caso encontre um valor absoluto maior
            if(fabs(a[i][i])<fabs(a[k][i])){
                //troca as linhas
                for(j=0;j<m;j++){                
                    double temp;
                    temp=a[i][j];
                    a[i][j]=a[k][j];
                    a[k][j]=temp;
                }
                double temp = b[k];
                b[k] = b[i];
                b[i] = temp;
            }
        }
        //Begin Gauss Elimination
        for(k=i+1;k<m;k++){ //para cada linha
            double  term = a[k][i]/ a[i][i];
            for(j=0;j<m;j++){ //para cada coluna
                a[k][j] = a[k][j] - (term * a[i][j]);
            }
            b[k] = b[k] - (term*b[i]);
        }
    }         
}

int main(int argc, void **argv)
{
    int printMatrix = 0;
    int printResult = 0;
    if(argc > 1 && !strcmp(argv[1], "-p")) //parâmetro "-p" determina a impressão de todas as matrizes
        printMatrix = 1;
    if(argc > 1 && !strcmp(argv[1], "-r")) //parâmetro "-r" determina a impressão apenas do resultado
        printResult = 1;

    int m;
    double **a;
    double *b;
    getMatrixAandB(&m, &a, &b); //realiza a leitura do arquivo
    double *x = malloc(m*sizeof(double)); //aloca espaço para o vetor de variáveis


    double time;
    clock_t Ticks[2];
    Ticks[0] = clock(); //marca inicio da contagem de tempo


    if(printMatrix){ //imprime matriz estendida A
        printf("*** Matriz A|B ****\n");
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

    if(printMatrix){ //imprime a matriz triangulada
        printf("*** Matriz Triangulada ****\n");
        for(int i=0; i < m; i++)
        {
            for(int j=0; j < m; j++)
                printf("%lf ", a[i][j]);
            printf("%lf\n", b[i]);
        }
        printf("\n");
    }

    if(printMatrix || printResult){ //imprime a matriz de variáveis
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

    Ticks[1] = clock(); //marca fim da contagem
    time = ((double) (Ticks[1] - Ticks[0])) / CLOCKS_PER_SEC; //calcula total
    printf("**** Cálculo realizado em %lf segundos\n", time);

    return 0;

}