
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>

//Funçao que aplica o método de Gauss
void gauss(int n_iteracao, mpf_t pi){
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
    for(int i = 0; i < n_iteracao; i++){
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

    mpf_t pi_gauss; //Variável que guardará a aproximação de pi
    gauss(pow(10,5),pi_gauss);
    
    //Printando no terminal
    //gmp_printf("Aproximacao Gauss-Legendre: %.6Ff \n",pi_gauss);

    //Printando no arquivo de saída
    gmp_fprintf(arq_saida,"Aproximacao Gauss-Legendre: %.6Ff \n",pi_gauss);

    fclose(arq_entrada);
    fclose(arq_saida);

    return 0;
}