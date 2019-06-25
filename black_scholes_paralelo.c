#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <pthread.h>

void *func(void *variavel);

double trials[100000];
double S, E, r, sigma, T, t, aux;
int M, i, k, j;

void blackScholes()
{

    int THREADS_MAX = 10;
    pthread_t threads[THREADS_MAX]; //declara as threads
    double thread_args[THREADS_MAX];

    FILE *arquivo = fopen("entrada_blackscholes.txt", "r");
    if (arquivo == NULL)
    {
        printf("Erro ao abrir o arquivo");
        return;
    }

    double variaveis[6];

    for (int i = 0; i < 6; i++)
    {
        fscanf(arquivo, "%lf\n", &variaveis[i]);
    }

    double S = variaveis[0];
    double E = variaveis[1];
    double r = variaveis[2];
    double sigma = variaveis[3];
    double T = variaveis[4];
    int M = variaveis[5];
    double aux1, aux2, med, d_p, amplitude, min, max;
    double trials[M];
    int gseed = time(NULL);

    for (k = 0; k < M / THREADS_MAX; k++)
    {
        j = THREADS_MAX * k;
        for (i = 0; i < THREADS_MAX; i++)
        {
            thread_args[i] = j;
            j = j + 1.0;
            //printf("Valor enviado para a thread: %lf\n", thread_args[i]);
            pthread_create(&threads[i], NULL, func, (void *)&thread_args[i]);
        }
        for (i = 0; i < THREADS_MAX; i++)
        {
            pthread_join(threads[i], NULL);
        }
    }

    double aux_mean = 0;

    for (i = 0; i < M; i++){
        aux_mean += trials[i];
    }
    med = aux_mean / (double)M;

    double aux_stddev = 0;
    for (i = 0; i < M; i++){
        aux_stddev += pow(trials[i] - med, 2);
    }
    d_p = sqrt(aux_stddev / (double)M);

    amplitude = 1.96 * d_p / sqrt(M);
    min = med - amplitude;
    max = med + amplitude;

    FILE *arquivoSaida = fopen("blackScholes_out.txt", "w");
    if (arquivoSaida == NULL)
    {
        printf("Erro ao escrever no arquivo");
        return;
    }

    fprintf(arquivoSaida, "S = %lf \nE = %lf \nr = %lf \nsigma = %lf \nT = %lf \nM = %d \nConfidence interval: (%lf, %lf)", S, E, r, sigma, T, M, min, max);
}

int main(int argc, char const *argv[])
{

    if (argc != 3)
    { //Testo se está passando o número correto de parâmetros
        printf("Erro na passagem de parametros !");
        return -1;
    }
    FILE *arq_entrada = fopen(argv[1], "r"); //Abro arquivo de entrada para leitura
    FILE *arq_saida = fopen(argv[2], "a+");  //Abro o arquivo de saída para escrita

    if (arq_entrada == NULL || arq_saida == NULL)
    { //Testo se os arquivos existem/ foram abertos corretamento
        printf("Erro ao abrir arquivos !");
        return -1;
    }

    fclose(arq_entrada);
    fclose(arq_saida);

    srand(time(NULL));
    blackScholes();

    return 0;
}

void *func(void *variavel)
{
    double *var = (double *)variavel;
    int i = (int)*var;
    double aux = ((double)rand() / (double)RAND_MAX);
    double t = S * exp(((r - (0.5 * pow(sigma, 2))) * T) + (sigma * sqrt(T) * aux));

    if (t - E > 0)
    {
        trials[i] = exp((-r) * T) * (t - E);
    }
    else
    {
        trials[i] = 0;
    }
    return NULL;
}