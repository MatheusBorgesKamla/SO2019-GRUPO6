#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <math.h>
#include <pthread.h>
#define N 100000
int control_t1 = 0;
int control_t2 = 0;
int control_t3 = 0;
int control_t4 = 0;

pthread_mutex_t mutex;

typedef struct Arg_borwein{

    mpf_t aux, aux1, aux2, aux3, aux4, aux5, aux_k, soma, sub;

} Arg_Borwein;

void *termo_t1(void *arg);
void *termo_t2(void *arg);
void *termo_t3(void *arg);
void *termo_t4(void *arg);
void *termo_t5(void *arg);

int main(int argc, char const *argv[])
{

    if (argc != 3)
    {
        printf("Erro na passagem de parametros !");
        return -1;
    }
    FILE *arq_entrada = fopen(argv[1], "r");
    FILE *arq_saida = fopen(argv[2], "w");

    if (arq_entrada == NULL || arq_saida == NULL)
    {
        printf("Erro ao abrir arquivos !");
        return -1;
    }

    Arg_Borwein aux;
    mpf_init(aux.aux);
    mpf_init(aux.aux1);
    mpf_init(aux.aux2);
    mpf_init(aux.aux3);
    mpf_init(aux.aux4);
    mpf_init(aux.aux5);
    mpf_init(aux.aux_k);
    mpf_init(aux.soma);
    mpf_init(aux.sub);

    pthread_mutex_init(&mutex, NULL);
    pthread_t t[5];

    pthread_create(&t[0], NULL, termo_t1, &aux);
    pthread_create(&t[1], NULL, termo_t2, &aux);
    pthread_create(&t[2], NULL, termo_t3, &aux);
    pthread_create(&t[3], NULL, termo_t4, &aux);
    pthread_create(&t[4], NULL, termo_t5, &aux);

    pthread_join(t[0], NULL);
    pthread_join(t[1], NULL);
    pthread_join(t[2], NULL);
    pthread_join(t[3], NULL);
    pthread_join(t[4], NULL);



    gmp_printf("Aproximacao Borwein: %.6Ff \n", aux.soma);

    fclose(arq_entrada);
    fclose(arq_saida);

    return 0;
}

void *termo_t1(void *arg)
{
    Arg_Borwein *aux = (Arg_Borwein *)arg;
    int i = 0;
    while (i < N)
    {
        if (control_t1 == 0)
        {
            pthread_mutex_lock(&mutex);

            //aux = 8*i + 1;
            mpf_set_ui(aux->aux, 8 * i);
            mpf_add_ui(aux->aux, aux->aux, 1);

            //aux2 = 4/aux;
            mpf_set_ui(aux->aux_k, 4);
            mpf_div(aux->aux2, aux->aux_k, aux->aux);

            control_t1 = 1;
            i++;

            pthread_mutex_unlock(&mutex);

        }
    }
    return (void *)aux->aux2;
}

void *termo_t2(void *arg)
{
    Arg_Borwein *aux = (Arg_Borwein *)arg;
    int i = 0;
    while (i < N)
    {
        if (control_t2 == 0)
        {
            pthread_mutex_lock(&mutex);
            //aux = 8*i + 4;
            mpf_set_ui(aux->aux, 8 * i);
            mpf_add_ui(aux->aux, aux->aux, 4);

            //aux3 = 2/aux;
            mpf_set_ui(aux->aux_k, 2);
            mpf_div(aux->aux3, aux->aux_k, aux->aux);

            control_t2 = 1;
            i++;
            pthread_mutex_unlock(&mutex);
        }
    }
    return (void *)aux->aux3;
}

void *termo_t3(void *arg)
{
    Arg_Borwein *aux = (Arg_Borwein *)arg;
    int i = 0;
    while (i < N)
    {
        if (control_t3 == 0)
        {
            pthread_mutex_lock(&mutex);
            //aux = 8*i + 5;
            mpf_set_ui(aux->aux, 8 * i);
            mpf_add_ui(aux->aux, aux->aux, 5);

            //aux4 = 1/aux;
            mpf_set_ui(aux->aux_k, 1);
            mpf_div(aux->aux4, aux->aux_k, aux->aux);

            control_t3 = 1;
            i++;
            pthread_mutex_unlock(&mutex);
        }
    }
    return (void *)aux->aux4;
}

void *termo_t4(void *arg)
{
    Arg_Borwein *aux = (Arg_Borwein *)arg;
    int i = 0;
    while (i < N)
    {
        if (control_t4 == 0)
        {
            pthread_mutex_lock(&mutex);
            //aux = 8*i + 6;
            mpf_set_ui(aux->aux, 8 * i);
            mpf_add_ui(aux->aux, aux->aux, 6);

            //aux5 = 1/aux;
            mpf_set_ui(aux->aux_k, 1);
            mpf_div(aux->aux5, aux->aux_k, aux->aux);

            control_t4 = 1;
            i++;
            pthread_mutex_unlock(&mutex);
        }
    }
    return (void *)aux->aux5;
}

void *termo_t5(void *arg)
{
    Arg_Borwein *aux = (Arg_Borwein *)arg;
    int i = 0;
    while (i < N)
    {
        if (control_t1 == 1 && control_t2 == 1 && control_t3 == 1 && control_t4 == 1)
        {
            pthread_mutex_lock(&mutex);
            //aux1 = (1/pow(16,i))*(aux2 - aux3 - aux4 - aux5);
            mpf_sub(aux->sub, aux->aux2, aux->aux3);
            mpf_sub(aux->sub, aux->sub, aux->aux4);
            mpf_sub(aux->sub, aux->sub, aux->aux5);
            mpf_set_ui(aux->aux1, 16);
            mpf_pow_ui(aux->aux1, aux->aux1, i);
            mpf_set_ui(aux->aux_k,1);
            mpf_div(aux->aux1, aux->aux_k, aux->aux1);
            mpf_mul(aux->aux1, aux->aux1, aux->sub);

            //soma += aux1;
            mpf_add(aux->soma, aux->soma, aux->aux1);


            i++;

            control_t1 = 0;
            control_t2 = 0;
            control_t3 = 0;
            control_t4 = 0;
            pthread_mutex_unlock(&mutex);
        }
    }
    return (void *)aux->soma;
}