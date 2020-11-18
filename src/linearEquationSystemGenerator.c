#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#define VAL_RANGE 100 //intervalo de valores possíveis
#define DECIMAL_FACTOR 100 //duas casas decimais

int main(int argc, void **argv)
{
    int equationCount = 0; //numero de equações do sistema
    char *eptr; //EndPointer utilizado no método strtod
    int i,j; //usados nos loops

    if(argc < 3){ //testa para caso haja erro na entrada dos parâmetros
        printf("Use \"./linearEquationSystemGenerator <x1> <x2> (<x3> ... <xn>)\", definindo os valores para ao menos 2 variáveis.\n\n");
    }
    double* X;
    if(!strcmp(argv[1], "-n"))
    {
        equationCount = atoi(argv[2]);
        X = (double*) malloc(equationCount*sizeof(double)); //aloca espaço para o vetor de variáveis
        for(i = 0; i < equationCount; i++){ //percorre os argumentos de entrada e constroi o vetor de variáveis
            X[i] = (double) i+1; //converte o valor para Double
        }
    }
    else
    {
        equationCount = (argc - 1); //define o numero de equações que o sistema terá (baseado no numero de variáveis)
        X = (double*) malloc(equationCount*sizeof(double)); //aloca espaço para o vetor de variáveis
        for(i = 0; i < equationCount; i++){ //percorre os argumentos de entrada e constroi o vetor de variáveis
            X[i] = strtod(argv[i + 1], &eptr); //converte o valor para Double
        }
    }

    double** matrixA = (double**) malloc(equationCount*sizeof(double*)); //aloca espaço para a matriz de constantes
    for(int j=0; j<equationCount; j++)
        matrixA[j] = (double*) malloc(equationCount*sizeof(double));
    double* matrixB = (double*) malloc(equationCount*sizeof(double)); //aloca espaço para a matriz resultado
    int rangeWithDecimal = VAL_RANGE * DECIMAL_FACTOR; //deve ser dividido por DECIMAL_FACTOR depois
    double sum = 0.0; //calcula o valor de B para cada equação
    double constant; //calcula cada um dos coeficientes

    srand(time(NULL)); //define semente para geração aleatória

    /*
    * Gera os valores das matrizes A e B
    * */
    for(i = 0; i < equationCount; i++) //gera cada equação
    {
        for(j = 0; j < equationCount; j++) //gera cada constante de uma equação
        {
            constant = (double) (rand() % rangeWithDecimal) / DECIMAL_FACTOR;
            matrixA[i][j] = constant;
            sum += constant * X[j];
        }
        matrixB[i] = sum; //define o valor de B para a equação i
        sum = 0.0;
    }

    
    /*
    * Abre o arquivo para guardar os valores calculados
    * */
    FILE *arq;
    if((arq = fopen("randomMatrix.txt", "w")) == NULL) //Abre o arquivo
    {
        //printf("Erro na abertura do arquivo randomMatrix.txt\n");
        exit(1);    
    }

    /*
    * Imprime as matrizes geradas na tela e no arquivo randomMatrix.txt
    * */
    fprintf(arq, "%d\n", equationCount); //primeira linha define o número de equações

    printf("--- Matriz de coeficientes A ----\n");
    for(i = 0; i < equationCount; i++)
    {
        for(j = 0; j < equationCount; j++)
        {
            printf("%lf\t", matrixA[i][j]); //imprime um valor na tela
            fprintf(arq, "%0.5lf ", matrixA[i][j]); //guarda o coeficiente no arquivo
        }
        printf("\n");
        fprintf(arq, "\n"); //adiciona uma quebra de linha no arquivo
    }
    fprintf(arq, "\n"); //adiciona uma linha em branco separando as matrizes

    printf("\n--- Matriz de resultado B ----\n");
    for(i = 0; i < equationCount; i++){
        printf("%lf\n", matrixB[i]);
        fprintf(arq, "%0.5lf ", matrixB[i]); //guarda o valor B no arquivo
    }
    fclose(arq); //fecha o arquivo
}