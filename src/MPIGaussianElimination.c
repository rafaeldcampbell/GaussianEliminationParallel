#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../lib/linearEquationFileReader.h"
#include <mpi.h>

//estrutura utilizada no findPivotIndex
typedef struct pivotValueIndex{ double val; int index;} PIVOTVALUEINDEX;

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

/*
* Realiza a troca de duas linhas (from e to) nas matrizes A e B.
* Esse processo, por ocorrer apenas uma vez por iteração, não será paralelizado
* */
void switchLines(double ***a, double **b, int from, int to)
{
    //trocando linhas de A
    double *temp = (*a)[from];
    (*a)[from] = (*a)[to];
    (*a)[to] = temp;

    //trocando linhas de B
    double tempB = (*b)[from];
    (*b)[from] = (*b)[to];
    (*b)[to] = tempB;
}

/*
* Atualiza a linha "where" de A e B com os valores de NewA e NewB.
* Esse processo, sendo uma simples atribuição, não será paralelizado.
* */
void updateLine(int m, double ***a, double **b, double *newA, double newB, int where)
{
    //trocando linhas de A
    memcpy((*a)[where], newA, m*sizeof(double));

    //trocando linhas de B
    (*b)[where] = newB;
}

/*
* Determina a linha do valor mais expressivo da coluna "column", considerando apenas os valores a partir do atual pivot ("column").
* Esse processo será paralelizado nos mesmos moldes de uma busca regular.
* */
int findPivotIndex(int m, int column, double **a, int size, int rank)
{
    int sizeSearchSubvector = m - column; //indica entre quantos elementos será feita a busca
    int chunk = sizeSearchSubvector / size;
    PIVOTVALUEINDEX *subvector;
    PIVOTVALUEINDEX *mysubvector = (PIVOTVALUEINDEX*) malloc(chunk*sizeof(PIVOTVALUEINDEX)); //alocando espaço para o vetor onde será feita a busca    
    
    if(rank == 0) //mestre
    {
        subvector = (PIVOTVALUEINDEX*) malloc(sizeSearchSubvector*sizeof(PIVOTVALUEINDEX)); //alocando espaço para o vetor onde será feita a busca
        for(int i=column; i < m; i++){ //constroi o vetor onde será realizada a busca
            subvector[i-column].val = abs(a[i][column]); //valor
            subvector[i-column].index = i; //index
        }
    }

    MPI_Scatter(subvector, chunk, MPI_DOUBLE_INT, mysubvector, chunk, MPI_DOUBLE_INT, 0, MPI_COMM_WORLD); //envia uma parte (chunk) igual do vetor
    
    PIVOTVALUEINDEX bigger = mysubvector[0];
    PIVOTVALUEINDEX pivot;
    for(int i = 0; i < chunk; i++) //procura, dentro de sua parte, o maior valor
    {
        if(abs(mysubvector[i].val) > abs(bigger.val)) //se encontra um valor maior, atualiza as referências
        {
            bigger = mysubvector[i]; //maior elemento
        }
    }
    if(rank == 0) //o processo mestre procura entre os que sobraram
    {
        for(int i=(size*chunk); i<sizeSearchSubvector; i++)
        {
            if(abs(subvector[i].val) > abs(bigger.val)) //se encontra um valor maior, atualiza as referências
            {
                bigger = subvector[i]; //maior elemento
            }
        }
    }

    MPI_Reduce(&bigger, &pivot, 1, MPI_DOUBLE_INT, MPI_MAXLOC, 0, MPI_COMM_WORLD); //define o maior valor
    MPI_Bcast(&pivot, 1, MPI_DOUBLE_INT, 0, MPI_COMM_WORLD); //garante que todos estejam com o mesmo pivô
    free(mysubvector);
    if(rank == 0)
        free(subvector);
    return pivot.index;
}

/*
* Executa a eliminação gaussiana com pivoteamento parcial.
* Esse processo foi paralelizado usando o método "cyclic striped".
* */
void gaussElimination(int m, double **a, double *b, int size, int rank){
    // A[linhas][colunas]
    // A -- matriz de coeficientes
    // B -- vetor de solução
    int c, i, j, pivotIndex, currentLine; //variáveis de apoio
    int loop; //determina quantos ciclos de envio serão realizados
    double *recvLine = (double*) malloc(m*sizeof(double)); //linha de apoio para receber resultado da normalização
    double *pivotLine; //linha de apoio que receberá o a linha pivot para normalização
    if(rank != 0) //o mestre apenas aponta e envia, logo, não precisa alocar.
        pivotLine = (double*) malloc(m*sizeof(double)); //os escravos, como recebem por mensagem, precisam de espaço alocado
    double bValue, bPivot, term; //valores utilizados no cálculo da linha normalizada
    MPI_Status status;

    for(c=0;c<m-1;c++){ //para cada coluna
        pivotIndex = findPivotIndex(m, c, a, size, rank); //busca do pivoteamento parcial
        loop = ceil((double) (m - c) / (size-1)); //determina quantos ciclos de envio serão realizados

        if(rank == 0){ //mestre
            if(pivotIndex != c) //caso o pivot real seja diferente do atual
                switchLines(&a, &b, c, pivotIndex); //realizando a troca do pivoteamento parcial

            pivotLine = a[c]; //definindo a linha pivot
            bPivot = b[c]; //definindo o B pivot
            /*
            * Nesta etapa, dá-se início ao processo de cyclic Striped. 
            * As linhas serão enviadas aos processos escravos de forma alternada, em loops.
            * Primeiro, será enviado o índex da linha. Caso seja válido, envia-se a linha.
            * Ao final de um loop de envios, executa-se um loop de recebimentos seguido das devidas trocas no vetor original.
            * */
            for(i=0; i<loop; i++) //executa "loop" vezes um ciclo de envios para todos os processos escravos
            {
                for(int rankDest=1; rankDest<size; rankDest++) //executa cada envio do ciclo
                {
                    currentLine = ((i*(size-1)) + rankDest) + c;//determina a linha que será enviada

                    MPI_Send(&currentLine, 1, MPI_INT, rankDest, 0, MPI_COMM_WORLD); //envia o índex da linha que será enviada
                    if(currentLine < m) //caso a linha atual seja válida
                    {
                        term = a[currentLine][c] / a[c][c]; //definindo o coeficiente
                        MPI_Send(&term, 1, MPI_DOUBLE, rankDest, 1, MPI_COMM_WORLD); //envia o coeficiente
                        MPI_Send(a[currentLine], m, MPI_DOUBLE, rankDest, 2, MPI_COMM_WORLD); //envia a linha
                        MPI_Send(pivotLine, m, MPI_DOUBLE, rankDest, 2, MPI_COMM_WORLD); //envia a linha pivo
                        MPI_Send(&(b[currentLine]), 1, MPI_DOUBLE, rankDest, 3, MPI_COMM_WORLD); //envia o valor de B
                        MPI_Send(&bPivot, 1, MPI_DOUBLE, rankDest, 3, MPI_COMM_WORLD); //envia o valor de B pivo
                    }
                }
                for(j=1; j<size; j++) //recebe os valores de volta
                {
                    MPI_Recv(&currentLine, 1, MPI_INT, j, 0, MPI_COMM_WORLD, &status); //recebe o índex da linha
                    if(currentLine < m) //caso seja uma linha válida
                    {
                        MPI_Recv(recvLine, m, MPI_DOUBLE, j, 2, MPI_COMM_WORLD, &status); //recebe a linha atualizada
                        MPI_Recv(&bValue, 1, MPI_DOUBLE, j, 3, MPI_COMM_WORLD, &status); //recebe B atualizado
                        updateLine(m, &a, &b, recvLine, bValue, currentLine); //atualiza a matriz original
                    }
                }
            }
        }
        else //trabalhadores
        {
            for(i=0; i<loop; i++) //executa "loop" vezes um ciclo de envios para todos os processos escravos
            {
                MPI_Recv(&currentLine, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status); //recebe o índex da linha
                if(currentLine < m) //caso seja uma linha válida
                {
                    MPI_Recv(&term, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status); //recebe o coeficiente
                    MPI_Recv(recvLine, m, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, &status); //recebe a linha
                    MPI_Recv(pivotLine, m, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, &status); //recebe a linha pivot
                    MPI_Recv(&bValue, 1, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD, &status); //recebe B
                    MPI_Recv(&bPivot, 1, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD, &status); //recebe B pivot
                    for(j=0;j<m;j++){ //calcula cada coluna
                        recvLine[j] = recvLine[j] - (term * pivotLine[j]); //atualiza as colunas de A
                    }
                    bValue =  bValue - (term*bPivot); //atualiza B
                }
                MPI_Send(&currentLine, 1, MPI_INT, 0, 0, MPI_COMM_WORLD); //envia o índex da linha que será enviada
                
                if(currentLine < m) //caso a linha atual seja válida
                {
                    MPI_Send(recvLine, m, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD); //envia a linha atualizada
                    MPI_Send(&bValue, 1, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD); //envia o valor de B atualizado
                }
            }
        }
    }
    free(recvLine);
    if(rank != 0) //apenas escravos alocam espaço para esse vetor, logo, apenas eles precisam liberar.
        free(pivotLine);
}

int main(int argc, char **argv)
{

    int rank,  size;
    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &size );


    MPI_Barrier(MPI_COMM_WORLD); //sincroniza processos para iniciar contagem do tempo
    double time = -MPI_Wtime();

    int printMatrix = 0;
    int printResult = 0;
    if(argc > 1 && !strcmp(argv[1], "-p"))
        printMatrix = 1;
    if(argc > 1 && !strcmp(argv[1], "-r")) //parâmetro "-r" determina a impressão apenas do resultado
        printResult = 1;

    int m; //tamanho da matriz
    double **a; //matriz de constantes
    double *b; //vetor de resultado
    double *x; //vetor de variáveis

    if(rank == 0)
        getMatrixAandB(&m, &a, &b);
    
    MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD); //garantindo que todos recebam m
    
    if(rank == 0 && printMatrix)
    {
        printf("*** Matriz A|B ****\n"); //imprime a matriz estendida no estado atual
        for(int i=0; i < m; i++)
        {
            for(int j=0; j < m; j++)
                printf("\t%lf ", a[i][j]);
            printf("\t%lf\n", b[i]);
        }
        printf("\n");
    }

    gaussElimination(m, a, b, size, rank); //realiza a eliminação gaussiana
    
    if(rank == 0)
    {
        x = malloc(m*sizeof(double)); //aloca espaço para o vetor de X
        backSubstitution(m, a, b, x);

        if(printMatrix)
        {
            printf("*** Matriz Triangulada ****\n"); //imprime a matriz já triangulada
            for(int i=0; i < m; i++)
            {
                for(int j=0; j < m; j++)
                    printf("\t%lf ", a[i][j]);
                printf("\t%lf\n", b[i]);
            }
            printf("\n");
        }

        if(printMatrix || printResult)
        {
            printf("*** Vetor X ***\n");
            for(int j=0; j < m; j++)
                printf("X%d = %lf\n", j, x[j]);
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
    }

    MPI_Barrier(MPI_COMM_WORLD); //sincroniza processos para fim da contagem do tempo
    time += MPI_Wtime();

    if(rank == 0) //imprime o tempo
        printf("**** Cálculo realizado em %lf segundos\n", time);

    MPI_Finalize();
    return 0;

}