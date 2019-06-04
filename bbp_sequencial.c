#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>

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

    mpf_t pi_bbp;
    bbp(pow(10,5),pi_bbp);
    
    //Printando no terminal
    gmp_printf("Aproximacao Borwein: %.6Ff \n",pi_bbp);

    //Printando no arquivo de saída
    gmp_fprintf(arq_saida,"Aproximacao Borwein: %.6Ff \n",pi_bbp);

    fclose(arq_entrada);
    fclose(arq_saida);

    return 0;
}