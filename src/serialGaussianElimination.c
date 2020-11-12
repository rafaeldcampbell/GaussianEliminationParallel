#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../lib/linearEquationFileReader.h"


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
    int m;
    double **a;
    double *b;
    getMatrixAandB(&m, &a, &b);
    double *x = malloc(m*sizeof(double));

    printf("*** Matriz ****\n");
    for(int i=0; i < m; i++)
    {
        for(int j=0; j < m; j++)
            printf("%lf ", a[i][j]);
        printf("%lf\n", b[i]);
    }
    printf("\n");

    gaussElimination(m, a, b);
    backSubstitution(m, a, b, x);

    printf("*** Vetor X ***\n");
    for(int j=0; j < m; j++)
    {
        printf("X%d = %lf\n", j, x[j]);
    }

    return 0;

}