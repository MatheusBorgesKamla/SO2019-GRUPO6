#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>
#include <time.h>
#include <pthread.h>
#define N_threads 10
#define iteracoes 1000000000

int rand_i = 0;

void *mc(void *arg){
    rand_i += 1000;
    unsigned int gseed = time(NULL) - rand_i; // Varíavel gseed sempre vai pegar um valor deferente para cada execução fazendo com que varie a sequência de nros aleatorios gerados
    srand(gseed); // Função que vai inicializar a rand_r utilizando o gseed para determinar a sequência de aleatórios gerados
    int teste;
    int *count = (int *)malloc(sizeof(int)); 
    

    mpf_set_default_prec(pow(10,5));
    mpf_t x,y,z; //Variáveis x e y representam as coordenadas dos pontos gerados aleatoriamente
    mpf_init(x);
    mpf_init(y);
    mpf_init(z);

    for(int i=0; i<iteracoes/N_threads; i++){
        //Para cada interação gero um valor aletório para as coordenadas com a função rand_r (foi checado sua documentação e ela é thread-safe) e divido pelo valor
        //máximo aleatório que posso gerar (RAND_MAX) para sempre gerar os números entre 0 e 1
        mpf_set_d(x,(((double)rand_r(&gseed)) / ((double)RAND_MAX)));
        mpf_set_d(y,(((double)rand_r(&gseed)) / ((double)RAND_MAX)));
        mpf_pow_ui(x,x,2);//x²
        mpf_pow_ui(y,y,2);//y²
        mpf_add(z,x,y);//z = x² + y² - calculo sua distância a partir da origem para se está contido na circunferência de raio 1
        teste = mpf_cmp_d(z,1.0); // Se z<1 - teste recebe um valor negativo, se z == 1 - teste recebe 0, se z>1 - teste recebe um valor positivo
        if(teste <= 0){
            *count += 1; // Contabilizo meu contador se o ponto estiver contido
        }
    }
    //printf("Valor de count = %d \n",*count);
    return count;
}

int main(int argc, char const *argv[])
{
    if(argc != 3){
        printf("Erro na passagem de parametros !");
        return -1;
    }
    FILE *arq_entrada = fopen(argv[1],"r");
    FILE *arq_saida = fopen(argv[2],"a+");

    if(arq_entrada == NULL || arq_saida == NULL){
        printf("Erro ao abrir arquivos !");
        return -1;
    }

    mpf_t pi;
    pthread_t *t = malloc(sizeof(pthread_t)*N_threads);
    int *counters;
    int i, cont_total;

   for(i=0; i<N_threads; i++){
        if(pthread_create(&t[i], NULL, mc, NULL)){
            printf("Erro ao criar thread %d",i);
		    exit(EXIT_FAILURE);
        }
    }

    for(i=0; i<N_threads; i++){
        if(pthread_join(t[i], &counters)){
            printf("Erro ao criar thread %d",i);
		    exit(EXIT_FAILURE);
        }
        //printf("Valor de count = %d \n",*counters);
        cont_total += *counters;
    }
    mpf_set_d(pi,(double)cont_total / iteracoes * 4);
    gmp_printf("Aproximacao Monte Carlo: %.6Ff \n",pi);
    gmp_fprintf(arq_saida,"Aproximacao Monte Carlo: %.6Ff \n",pi);
    fclose(arq_entrada);
    fclose(arq_saida);
    
    return 0;
}
