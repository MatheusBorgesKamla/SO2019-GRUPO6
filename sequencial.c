#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>
#include <time.h>
 
//Funçao que aplica o método de Gauss
void gauss(int iteracao, mpf_t pi){
    //Valores inicias do método
    mpf_set_default_prec(pow(10,5));
    mpf_t a, b, t, p;
    mpf_init_set_ui(a,1); // a = 1
    mpf_init_set_d(b,1/sqrt(2.0));// b = 1/raiz(2)
    mpf_init_set_d(t,0.25);// t = 1/4
    mpf_init_set_ui(p,1);// p = 1
    mpf_t a_prox, b_prox, t_prox, p_prox, aux; //constantes que receberão a próxima interação e uma auxiliar
    mpf_init(pi);
    mpf_init(a_prox);
    mpf_init(b_prox);
    mpf_init(t_prox);
    mpf_init(p_prox);
    mpf_init(aux);
    //Para cada interação temos os seguintes passos     
    for(int i = 0; i < iteracao; i++){
        mpf_add(aux,b,a);
        mpf_div_ui(a_prox,aux,2);//a[n+1] = (a[n] + b[n])/2

        mpf_mul(aux,a,b);
        mpf_sqrt(b_prox,aux);//b[n+1] = raiz(a[n]*b[n])
        //mpf_clear(aux);

        mpf_sub(aux,a,a_prox);
        mpf_pow_ui(aux,aux,2);
        mpf_mul(aux,p,aux);
        mpf_sub(t_prox,t,aux);// t[n+1] = t[n] - p[n](a[n] - a[n+1])²
        //mpf_clear(aux);

        mpf_mul_ui(p_prox,p,2);// p[n+1] = 2*p[n]

        //atualiza para a proxima interação
        mpf_set(a,a_prox);
        mpf_set(b,b_prox);
        mpf_set(t,t_prox);
        mpf_set(p,p_prox);
    }
    mpf_add(aux,a,b);
    mpf_pow_ui(aux,aux,2);
    mpf_mul_ui(t,t,4);
    mpf_div(pi,aux,t);//(a+b)²/4*t
}

//Funçao que aplica o método de Borwein
void bbp(int k, mpf_t pi){
    mpf_set_default_prec(pow(10,5));
    mpf_t aux, aux1, aux2, aux3, aux4, aux5, aux_k, soma, sub;
    mpf_init(pi);
    mpf_init(aux);
    mpf_init(aux1);
    mpf_init(aux2);
    mpf_init(aux3);
    mpf_init(aux4);
    mpf_init(aux5);
    mpf_init(aux_k);
    mpf_init(soma);
    mpf_init(sub);
    //Para cada interação temos um Somatório de 0 até N do termo : (1/16^k)((4/(8k + 1)) - 2/(8k + 5) - 1/(8k + 5) - 1/(8k + 6))
    for (int i = 0; i <= k; i++)
    {
        //aux = 8*i + 1;
        mpf_set_ui(aux,8*i);
        mpf_add_ui(aux,aux,1);

        //aux2 = 4/aux;
        mpf_set_ui(aux_k,4);
        mpf_div(aux2,aux_k,aux);

        //aux = 8*i + 4;
        mpf_set_ui(aux,8*i);
        mpf_add_ui(aux,aux,4);

        //aux3 = 2/aux;
        mpf_set_ui(aux_k,2);
        mpf_div(aux3,aux_k,aux);

        //aux = 8*i + 5;
        mpf_set_ui(aux,8*i);
        mpf_add_ui(aux,aux,5);

        //aux4 = 1/aux;
        mpf_set_ui(aux_k,1);
        mpf_div(aux4,aux_k,aux);

        //aux = 8*i + 6;
        mpf_set_ui(aux,8*i);
        mpf_add_ui(aux,aux,6);

        //aux5 = 1/aux;
        mpf_set_ui(aux_k,1);
        mpf_div(aux5,aux_k,aux);

        //aux1 = (1/pow(16,i))*(aux2 - aux3 - aux4 - aux5);
        mpf_sub(sub,aux2,aux3);
        mpf_sub(sub,sub,aux4);
        mpf_sub(sub,sub,aux5);
        mpf_set_ui(aux1,16);
        mpf_pow_ui(aux1,aux1,i);
        mpf_div(aux1,aux_k,aux1);
        mpf_mul(aux1,aux1,sub);
        
        //soma += aux1;
        mpf_add(soma,soma,aux1);
    }
    //pi = soma;
    mpf_set(pi,soma);
}

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
    if(argc != 3){
        printf("Erro na passagem de parametros !");
        return -1;
    }
    FILE *arq_entrada = fopen(argv[1],"r");
    FILE *arq_saida = fopen(argv[2],"w");

    if(arq_entrada == NULL || arq_saida == NULL){
        printf("Erro ao abrir arquivos !");
        return -1;
    }

    mpf_t pi_gauss, pi_bbp, pi_mc;
    gauss(pow(10,5),pi_gauss);
    bbp(pow(10,5),pi_bbp);
    monte_carlo(pow(10,9),pi_mc);
 
    gmp_printf("Aproximacao Gauss-Legendre: %.6Ff \n",pi_gauss);
    gmp_printf("Aproximacao Borwein: %.6Ff \n",pi_bbp);
    gmp_printf("Aproximacao Monte Carlo: %.6Ff \n",pi_mc);

    gmp_fprintf(arq_saida,"Aproximacao Gauss-Legendre: %.6Ff \n",pi_gauss);
    gmp_fprintf(arq_saida,"Aproximacao Borwein: %.6Ff \n",pi_bbp);
    gmp_fprintf(arq_saida,"Aproximacao Monte Carlo: %.6Ff \n",pi_mc);

    fclose(arq_entrada);
    fclose(arq_saida);

    return 0;
}