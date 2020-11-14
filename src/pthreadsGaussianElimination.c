#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>

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
	
	printf("Chamou partial pivoting thread %i\n ", thread);
	
	int k, j;
	//Pivoteamento. Para separar a coluna de acordo com o número de threads, se usa a mesma ideia do cycling striped;
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
	
	printf("Chamou thread gauss elimination %i\n ", thread);
	
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
        
        double beginnnig = clock();
        
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
        
        double end = clock();
        
        printf("Tempo de duração eliminação gaussiana com pivoteamento pthreads %f\n", ((double) (end - beginnnig)) / CLOCKS_PER_SEC);

        free(param);
        
    }         
}



int main(int argc, void **argv)
{
    int m = 3;
    double **a = malloc(m*sizeof(double));
    for(int i = 0; i < m; i++)
    {
        a[i] = malloc(m*sizeof(double));
    }
    a[0][0] = 13.0;
    a[0][1] = 17.0;
    a[0][2] = 20.0;
    a[1][0] = 10.0;
    a[1][1] = 12.0;
    a[1][2] = 9.0;
    a[2][0] = 7.0;
    a[2][1] = 15.0;
    a[2][2] = 11.0;
    
    double *b = malloc(m*sizeof(double));
    b[0] = 9.0;
    b[1] = 7.0;
    b[2] = 2;
    double *k = malloc(m*sizeof(double));

    printf("*** Matriz ****\n");
    for(int i=0; i < m; i++)
    {
        for(int j=0; j < m; j++)
            printf("%lf ", a[i][j]);
        printf("%lf\n", b[i]);
    }
    printf("\n");

    gaussElimination(m, a, b);
    backSubstitution(m, a, b, k);

    printf("*** Vetor K ***\n");
    for(int j=0; j < m; j++)
    {
        printf("X%d = %lf\n", j, k[j]);
    }

    return 0;

}
