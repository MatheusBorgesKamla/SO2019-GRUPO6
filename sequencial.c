#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>
 
//Funçao que aplica o método de gauss
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
        mpf_set_ui(aux1,pow(16,i));
        mpf_div(aux1,aux_k,aux1);
        mpf_mul(aux1,aux1,sub);
        
        //soma += aux1;
        mpf_add(soma,soma,aux1);
    }
    //pi = soma;
    mpf_set(pi,soma);
}


int main() {
    int n = pow(10,2);
    
    mpf_t pi_gauss;
    gauss(n,pi_gauss);
    //int pre = mpf_get_prec(pi_gauss);
    gmp_printf("Aproximacao Gauss-Legendre: %.6Ff \n",pi_gauss);
    
    mpf_t pi_bbp;
    bbp(n,pi_bbp);
    gmp_printf("Aproximacao BBP: %.6Ff \n",pi_bbp);


    return 0;
}