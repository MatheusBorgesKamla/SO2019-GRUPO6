#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>
#include <time.h>

//Funçao que aplica o método de Monte Carlo
void monte_carlo(int n, mpf_t pi){
    unsigned int gseed = time(NULL); // Varíavel gseed sempre vai pegar um valor deferente para cada execução fazendo com que varie a sequência de nros aleatorios gerados
    srand(gseed); // Função que vai inicializar a rand_r utilizando o gseed para determinar a sequência de aleatórios gerados
    int teste, count = 0; 

    mpf_set_default_prec(pow(10,5));
    mpf_t x,y,z; //Variáveis x e y representam as coordenadas dos pontos gerados aleatoriamente
    mpf_init(pi);
    mpf_init(x);
    mpf_init(y);
    mpf_init(z);

    for(int i=0; i<n; i++){
        //Para cada interação gero um valor aletório para as coordenadas com a função rand_r (foi checado sua documentação e ela é thread-safe) e divido pelo valor
        //máximo aleatório que posso gerar (RAND_MAX) para sempre gerar os números entre 0 e 1
        mpf_set_d(x,(((double)rand_r(&gseed)) / ((double)RAND_MAX)));
        mpf_set_d(y,(((double)rand_r(&gseed)) / ((double)RAND_MAX)));
        mpf_pow_ui(x,x,2);//x²
        mpf_pow_ui(y,y,2);//y²
        mpf_add(z,x,y);//z = x² + y² - calculo sua distância a partir da origem para se está contido na circunferência de raio 1
        teste = mpf_cmp_d(z,1.0); // Se z<1 - teste recebe um valor negativo, se z == 1 - teste recebe 0, se z>1 - teste recebe um valor positivo
        if(teste <= 0)
            count++; // Contabilizo meu contador se o ponto estiver contido
    }
    mpf_set_d(pi,(double)count / n * 4);
}

int main(int argc, char const *argv[]) {
    if(argc != 3){ //Testo se está passando o número correto de parâmetros
        printf("Erro na passagem de parametros !");
        return -1;
    }
    FILE *arq_entrada = fopen(argv[1],"r"); //Abro arquivo de entrada para leitura
    FILE *arq_saida = fopen(argv[2],"a+"); //Abro o arquivo de saída para escrita

    if(arq_entrada == NULL || arq_saida == NULL){ //Testo se os arquivos existem/ foram abertos corretamento
        printf("Erro ao abrir arquivos !");
        return -1;
    }

    mpf_t pi_mc;
    monte_carlo(pow(10,9),pi_mc);
    
    //Printando no terminal
    gmp_printf("Aproximacao Monte Carlo: %.6Ff \n",pi_mc);

    //Printando no arquivo de saída
    gmp_fprintf(arq_saida,"Aproximacao Monte Carlo: %.6Ff \n",pi_mc);

    fclose(arq_entrada);
    fclose(arq_saida);

    return 0;
}