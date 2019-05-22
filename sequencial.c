#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>
 
//Funçao que aplica o método de gauss
void gauss(int iteracao, mpf_t pi){
    //Valores inicias do método
    mpf_set_default_prec(pow(10,5));
    mpf_t a, b, t, p;
    mpf_init_set_ui(a,1);
    mpf_init_set_d(b,1/sqrt(2.0));
    mpf_init_set_d(t,0.25);
    mpf_init_set_ui(p,1);
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
    mpf_div(pi,aux,t);//(a+b)²/4t*/
}

int main() {
    int i, n = pow(10,5);
    mpf_t pi;
    gauss(n,pi);
    int pre = mpf_get_prec(pi);
    gmp_printf("Aproximacao Gauss-Legendre: %.6Ff \n",pi);
    return 0;
}