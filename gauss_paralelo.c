#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>
#include <time.h>
#include <pthread.h>
#define n_gauss 100000
//Variáveis testes utilizadas para fazer com que as threads sejam executadas na ordem correta p/ compartilhamento d memória
int test_ab = 0, test_ta = 0, test_tp = 0, test_tb = 0;
//Mutex utilizada para garantir que o compartilhamento de memória vai ocorrer corretamente
pthread_mutex_t mutex;
//Struct utilizada para passagem de parâmetros para as threads
typedef struct arg_gauss
{
    mpf_t a0, b0, p0, t0;
    mpf_t a_ant, b_ant, p_ant, t_ant;
    mpf_t a_prox, b_prox, p_prox, t_prox;
} Arg_Gauss;

void *gauss_a(void *arg)
{
    Arg_Gauss *aux = (Arg_Gauss *)arg;
    int i = 1;
    while (i < n_gauss)
    {
        if (test_ab == 0 && test_ta == 0)
        {
            pthread_mutex_lock(&mutex);
            if(i != 1){
                mpf_set(aux->a0, aux->a_ant);
                mpf_set(aux->a_ant, aux->a_prox);
                mpf_add(aux->a_prox, aux->a_ant, aux->b_prox);
                mpf_div_ui(aux->a_prox, aux->a_prox, 2);
            }else{
                mpf_add(aux->a_prox, aux->a_ant, aux->b_ant);
                mpf_div_ui(aux->a_prox, aux->a_prox, 2);
            }
            pthread_mutex_unlock(&mutex);
            test_ab = 1;
            test_ta = 1;
            i++;
        }
    }
}
void *gauss_b(void *arg)
{
    Arg_Gauss *aux = (Arg_Gauss *)arg;
    int i = 1;
    mpf_t temp;
    mpf_init(temp);
    while (i < n_gauss)
    {
        if (test_ab == 1 && test_tb == 0)
        {
            pthread_mutex_lock(&mutex);
            if(i != 1){
                //Calculo do p
                mpf_set(aux->p0, aux->p_ant);
                mpf_set(aux->p_ant, aux->p_prox);
                //Calculo do b
                mpf_set(aux->b0, aux->b_ant);
                mpf_set(aux->b_ant, aux->b_prox);
            }
            //Calculo do p
            mpf_mul_ui(aux->p_prox, aux->p_ant, 2);
            //Calculo do b
            mpf_mul(temp, aux->a_ant, aux->b_ant);
            mpf_sqrt(aux->b_prox,temp);
            pthread_mutex_unlock(&mutex);
            test_ab = 0;
            test_tb = 1;
            i++;
        }
    }
}
/*void *gauss_p(void *arg)
{
    Arg_Gauss *aux = (Arg_Gauss *)arg;
    int i = 1;
    while (i < n_gauss)
    {
        if (test_tp == 0)
        {
            if(i != 1){
                mpf_set(aux->p0, aux->p_ant);
                mpf_set(aux->p_ant, aux->p_prox);
            }
            mpf_mul_ui(aux->p_prox, aux->p_ant, 2);
            test_tp = 1;
            i++;
        }
    }
}*/

/*void *gausst(void *arg)
{
    Arg_Gauss *aux = (Arg_Gauss *)arg;
    mpf_t temp;
    mpf_init(temp);
    int i = 0;
    while (i < n_gauss)
    {
        if (test_ta == 1 && test_tb == 1 && test_tp == 1)
        {
            pthread_mutex_lock(&mutex);
            mpf_sub(temp, aux->a_ant, aux->a_prox);
            mpf_pow_ui(temp, temp, 2);
            mpf_mul(temp, aux->p_ant, temp);
            mpf_sub(aux->t_prox, aux->t_ant, temp);
            mpf_set(aux->t_ant, aux->t_prox);
            pthread_mutex_unlock(&mutex);
            if (i < n_gauss - 2)
            {
                test_ta = 0;
                test_tb = 0;
                test_tp = 0;
            }
            //printf("Interacao T: %d\n", i);
            /*gmp_printf("a0: %.6Ff \n",aux->a0);
            gmp_printf("b0: %.6Ff \n",aux->b0);
            gmp_printf("p0: %.6Ff \n",aux->p0);
            gmp_printf("t0: %.6Ff \n",aux->t0);
            gmp_printf("a_ant: %.6Ff \n",aux->a_ant);
            gmp_printf("b_ant: %.6Ff \n",aux->b_ant);
            gmp_printf("p_ant: %.6Ff \n",aux->p_ant);
            gmp_printf("t_ant: %.6Ff \n",aux->t_ant);
            gmp_printf("a_prox: %.6Ff \n",aux->a_prox);
            gmp_printf("b_prox: %.6Ff \n",aux->b_prox);
            gmp_printf("p_prox: %.6Ff \n",aux->p_prox);
            gmp_printf("t_prox: %.6Ff \n",aux->t_prox);
            printf("---------------\n");
        
            i++;
        }
    }
}*/
void *gausst(void *arg)
{
    Arg_Gauss *aux = (Arg_Gauss *)arg;
    mpf_t temp;
    mpf_init(temp);
    int i = 0;
    while (i < n_gauss)
    {
        if (test_ta == 1 && test_tb == 1)
        {
            pthread_mutex_lock(&mutex);
            mpf_sub(temp, aux->a_ant, aux->a_prox);
            mpf_pow_ui(temp, temp, 2);
            mpf_mul(temp, aux->p_ant, temp);
            mpf_sub(aux->t_prox, aux->t_ant, temp);
            mpf_set(aux->t_ant, aux->t_prox);
            pthread_mutex_unlock(&mutex);
            if (i < n_gauss - 2)
            {
                test_ta = 0;
                test_tb = 0;
            }
            i++;
        }
    }
}

void gauss(mpf_t pi)
{
    Arg_Gauss arg;
    mpf_t temp;

    pthread_mutex_init(&mutex, NULL);
    mpf_set_default_prec(pow(10,5));
    //Valores inicias do método (n=0) e dos termos n=1 para conseguir paralelizar o t
    mpf_init(pi);
    mpf_init_set_ui(arg.a0, 1);
    mpf_init_set_ui(arg.p0, 1);
    mpf_init_set_d(arg.b0, 1 / sqrt(2));
    mpf_init_set_d(arg.t0, 0.25);
    mpf_init_set_d(arg.a_ant, (1 + (1 / sqrt(2))) / 2);
    mpf_init_set_d(arg.b_ant, sqrt(1 * (1 / sqrt(2))));
    mpf_init_set_d(arg.p_ant, 2);
    mpf_init(arg.t_ant);
    mpf_init(temp);

    mpf_sub(temp, arg.a_ant, arg.a0);
    mpf_pow_ui(temp, temp, 2);
    mpf_mul(temp, arg.p0, temp);
    mpf_sub(arg.t_ant, arg.t0, temp);
    mpf_init(arg.a_prox);
    mpf_init(arg.b_prox);
    mpf_init(arg.p_prox);
    mpf_init(arg.t_prox);

    pthread_t *t = malloc(sizeof(pthread_t) * 4);

    if (pthread_create(&t[0], NULL, gauss_a, &arg))
    {
        perror("Erro na thread 1\n");
        exit(EXIT_FAILURE);
    }
    //printf("Iniciou Gauss A\n");
    if (pthread_create(&t[1], NULL, gauss_b, &arg))
    {
        perror("Erro na thread 2\n");
        exit(EXIT_FAILURE);
    }
    /*//printf("Iniciou Gauss B\n");
    if (pthread_create(&t[2], NULL, gauss_p, &arg))
    {
        perror("Erro na thread 3\n");
        exit(EXIT_FAILURE);
    }*/
    //printf("Iniciou Gauss P\n");
    if (pthread_create(&t[3], NULL, gausst, &arg))
    {
        perror("Erro na thread 4\n");
        exit(EXIT_FAILURE);
    }
    //printf("Iniciou Gauss T\n");
    pthread_join(t[0], NULL);
    //printf("Terminou Gauss A\n");
    pthread_join(t[1], NULL);
    //printf("Terminou Gauss B\n");
    pthread_join(t[2], NULL);
    //printf("Terminou Gauss C\n");
    pthread_join(t[3], NULL);
    //printf("Terminou Gauss T\n");

    mpf_t aux1, aux2;
    mpf_init(aux1);
    mpf_init(aux2);

    mpf_add(aux1, arg.a_prox, arg.b_prox);
    mpf_pow_ui(aux1, aux1, 2);
    mpf_set(aux2, arg.t_prox);
    mpf_mul_ui(aux2, aux2, 4);
    mpf_div(pi, aux1, aux2);
}

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

    mpf_t pi_gauss;
    gauss(pi_gauss);
    gmp_printf("Aproximacao Gauss-Legendre: %.6Ff \n", pi_gauss);
    fclose(arq_entrada);
    fclose(arq_saida);
    return 0;
}
