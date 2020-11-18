#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <string.h>
#include <time.h>
#include "../lib/linearEquationFileReader.h"

#define NUM_THREADS 4

//estrutura para se passar pelo pthreads ao se chamar uma função
struct gauss_elimination_param{
	int col;
	int threadNumber;
	int matrixSize;
	double **matrix;
	double* arrayB;
};

void *partialPivotingPthread(void *param){
	//recupera as informações recebidas como parametro pelo pthreads_create
	struct gauss_elimination_param* threadParam = (struct gauss_elimination_param*) param;
	int i = threadParam->col;
	int thread = threadParam->threadNumber;
	double **a = threadParam->matrix;
	int m = threadParam->matrixSize;
	double *b = threadParam->arrayB;
	
	int k, j;
	//Para separar a coluna de acordo com o número de threads, se usa a mesma ideia do cycling striped;
	for(k=i+1 + thread;k<m;k = k+ NUM_THREADS){ // para cada linha
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
}

//função para ser chamada pelo pthreads
void *gaussPthread( void *param){

	//recupera as informações recebidas como parametro pelo pthreads_create
	struct gauss_elimination_param* threadParam = (struct gauss_elimination_param*) param;
	int i = threadParam->col;
	int thread = threadParam->threadNumber;
	double **a = threadParam->matrix;
	int m = threadParam->matrixSize;
	double *b = threadParam->arrayB;
	
	int k, j;
	
	//ELIMINAÇÃO GAUSSIANA utilizando cycling striped. Começa pulando "+thread" e salta pelo número de threads	
	for(k=i + 1 +thread; k<m;k = k + NUM_THREADS){ //para cada linha
        double  term = a[k][i]/ a[i][i];
        
        for(j=0;j<m;j++){ //para cada coluna
            
            a[k][j] = a[k][j] - (term * a[i][j]);	
        }
        
        b[k] = b[k] - (term*b[i]);
    }
}

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
    pthread_t threadArray[NUM_THREADS];

    int i = 0, j = 0, k = 0;
    for(i=0;i<m-1;i++){ //para cada coluna
        //Cria o vetor de parametros. Um parametro para cada thread que será inicializado
        struct gauss_elimination_param* param = malloc(NUM_THREADS * sizeof(struct gauss_elimination_param));
        
        int threadOpen, threadJoin;
        //Iniciliza os threads contendoa a função para fazer a eliminação gaussiana
        for(threadOpen = 0; threadOpen < NUM_THREADS; threadOpen++){
        	//passa os parametros para a estrutura
        	param[threadOpen].col = i;
        	param[threadOpen].threadNumber = threadOpen;
        	param[threadOpen].matrix = a;
        	param[threadOpen].matrixSize = m;
        	param[threadOpen].arrayB = b;
        	//inicializa um thread
        	pthread_create(&threadArray[threadOpen], NULL, partialPivotingPthread, &param[threadOpen]);
        }
        
        //espera todos os threads terminarem
        for(threadJoin = 0; threadJoin < NUM_THREADS; threadJoin++){
        	pthread_join(threadArray[threadJoin], NULL);
        }
        
        //Iniciliza os threads contendoa a função para fazer a eliminação gaussiana
        for(threadOpen = 0; threadOpen < NUM_THREADS; threadOpen++){
        	//passa os parametros para a estrutura
        	param[threadOpen].col = i;
        	param[threadOpen].threadNumber = threadOpen;
        	//inicializa um thread
        	pthread_create(&threadArray[threadOpen], NULL, gaussPthread, &param[threadOpen]);
        }
        
        //espera todos os threads terminarem
        for(threadJoin = 0; threadJoin < NUM_THREADS; threadJoin++){
        	pthread_join(threadArray[threadJoin], NULL);
        }
        
        free(param);
    }         
}



int main(int argc, void **argv)
{

    double time;
    clock_t Ticks[2];
    Ticks[0] = clock(); //marca inicio da contagem de tempo

    int printMatrix = 0;
    int printResult = 0;
    if(argc > 1 && !strcmp(argv[1], "-p")) //parâmetro "-P" determina a impressão de todas as matrizes
        printMatrix = 1;
    if(argc > 1 && !strcmp(argv[1], "-r")) //parâmetro "-r" determina a impressão apenas do resultado
        printResult = 1;

    int m;
    double **a;
    double *b;
    getMatrixAandB(&m, &a, &b); //realiza a leitura do arquivo
    double *x = malloc(m*sizeof(double)); //aloca espaço para o vetor de variáveis

    if(printMatrix){ //imprime a matriz A estendida
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

    if(printMatrix){ //imprime matriz triangulada
        printf("*** Matriz Triangulada ****\n");
        for(int i=0; i < m; i++)
        {
            for(int j=0; j < m; j++)
                printf("%lf ", a[i][j]);
            printf("%lf\n", b[i]);
        }
        printf("\n");
    }

    if(printMatrix || printResult){ //imprime matriz de variáveis X
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
