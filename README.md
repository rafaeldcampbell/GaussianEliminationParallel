# Eliminação Gaussiana com Pivoteamento Parcial
### Implementações paralelas para resolução de sistemas de equações lineares

Como parte da avaliação da disciplina de Laboratório de Programação Paralela -- da graduação em Ciência da Computação da Universidade Federal Fluminense (UFF) --, foi proposto a implementação de um algoritmo paralelo para resolução de sistemas de equações lineares por meio da técnica de eliminação gaussiana com pivoteamento parcial. Neste projeto, serão feitas quatro implementações (uma serial e três paralelas).

* Sequencial
* Pthreads
* MPI
* OpenMP

### Resolução de sistemas
Equações lineares são, sobretudo, equações polinomiais que podem ser escritas na forma <img src="https://render.githubusercontent.com/render/math?math=$a_1X_1 + a_2X_2 + ... + a_nX_n = b$"> , onde <img src="https://render.githubusercontent.com/render/math?math=$a_1,\ a_2,\ ...,\ a_n$"> e <img src="https://render.githubusercontent.com/render/math?math=$b$"> são constantes e <img src="https://render.githubusercontent.com/render/math?math=$X_1,\ X_2,\ ...,\ X_n$">  são incógnitas. Um sistema de equações lineares é, essencialmente, um conjunto de equações lineares que diz respeito ao mesmo conjunto de incógnitas, sendo normalmente escrito como:

<img src="https://render.githubusercontent.com/render/math?math=\begin{equation}\begin{cases}a_{11}x_1 + a_{12}x_2 + ... + a_{1n}x_n = b_1 \\a_{21}x_1 + a_{22}x_2 + ... + a_{2n}x_n = b_2 \\a_{31}x_1 + a_{32}x_2 + ... + a_{3n}x_n = b_3\end{cases}\end{equation}">

As técnicas de resolução se concentram em determinar valores para <img src="https://render.githubusercontent.com/render/math?math=$X_1,\ X_2,\ ...,\ X_n$"> de modo a satisfazer todas as equações de um sistema.

### A Eliminação Gaussiana
A resolução de sistemas pelo método da eliminação gaussiana tem, essencialmente, três passos:

1. Obter uma matriz aumentada <img src="https://render.githubusercontent.com/render/math?math=$[A|b]$">, onde o sistema de equação seja <img src="https://render.githubusercontent.com/render/math?math=$AX = b$">.
2. Obter uma matriz equivalente <img src="https://render.githubusercontent.com/render/math?math=$[\dot A|\dot b]$">, onde <img src="https://render.githubusercontent.com/render/math?math=$\dot A$"> seja triangular superior.
3. Resolver o sistema de equações <img src="https://render.githubusercontent.com/render/math?math=$\dot AX = \dot b$"> por substituição regressiva.

Entretanto, esse método tem três condições para garantir uma solução única:
* As matrizes devem ser densas e o pivô deve ser não-nulo.
* As equações devem ser linearmente independentes.
* O número de equações deve ser igual ao número de variáveis.

A segunda e a terceira condição serão tratadas pela geração das matrizes, mas a primeira adiciona um critério de dificuldade. Como forma de garantir que o pivô seja não nulo -- e evitar uma possível variação de escala --, usa-se a técnica de pivoteamento parcial, que consiste em dois passos:
1. A cada passo (ou coluna operada), escolhe-se o elemento de maior valor absoluto a partir da linha pivô.
2. Caso o elemento não seja o atual pivô, faz-se a permutação das linhas.

### Executando os algoritmos
Para executar cada um dos algoritmos, basta compilar e executar.

* Serial
```
gcc serialGaussianElimination.c -o serialGaussianElimination -lm && ./serialGaussianElimination
```

* Pthreads
```
gcc pthreadsGaussianElimination.c -o pthreadsGaussianElimination -lm -lpthread && ./pthreadsGaussianElimination
```

* MPI
```
mpicc MPIGaussianElimination.c -o MPIGaussianElimination -lm && mpirun -n <numberOfProcesses> MPIGaussianElimination
```

* OpenMP
```
gcc openMPGaussianElimination.c -o openMPGaussianElimination -lm -fopenmp && ./openMPGaussianElimination
```

### Parâmetros e resultado
Por padrão, todas as implementações retornaram apenas o tempo de cálculo.
```
./<implementacao>
```
>>>
\**** Cálculo realizado em 0.000513 segundos  
>>>

Alternativamente, pode-se adicionar o parâmetro **"-r"** para exibir o resultado ou **"-p"** para exibir as matrizes e o resultado.

```
./<implementacao> -r
```
>>>
\*** Vetor X ***  
    X0 = 2.200000  
    X1 = 4.100000  
    X2 = 5.000000  
    X3 = 2.200000  
 \**** Cálculo realizado em 0.000513 segundos
>>>

```
./<implementacao> -p
```
>>>    
\*** Matriz A|B ****  
59.220000 5.110000 68.360000 89.050000 688.945000  
55.790000 11.560000 83.860000 22.280000 638.450000  
13.080000 67.990000 98.950000 89.580000 999.361000  
93.270000 44.600000 30.620000 42.130000 633.840000  
  
*** Matriz Triangulada ****  
93.270000 44.600000 30.620000 42.130000 633.840000   
0.000000 61.735384 94.655912 83.671772 910.472534  
0.000000 0.000000 88.723794 17.569225 482.271264  
0.000000 0.000000 0.000000 77.021511 169.447324  
  
*** Vetor X ***  
X0 = 2.200000  
X1 = 4.100000  
X2 = 5.000000  
X3 = 2.200000  
**** Cálculo realizado em 0.001034 segundos  
>>>

### Mais informação
Para explicações mais detalhadas sobre execução, código e motivação, consulte o material disponível em [Relatorio_final.pdf](./Relatorio_final.pdf).
