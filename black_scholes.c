#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

double media(double trials[], unsigned int m){
    double total = 0;
    double media;
    for (int i = 0; i < m; i++)
    {
        total += trials[i];
    }
    //printf("Média é = %lf", total/m);
    media = total / m;
    return media;
}

double desvio_padrao(double trials[], unsigned int m, double media){
    double p, dp;

    for(int i=0; i<m; i++){
        p = p + pow(trials[i]-media,2);
    }
    dp=sqrt(p/(m-1));
    return dp;

}

void blackscholes(){
    FILE *arquivo = fopen("entrada_blackscholes.txt", "r");
    if(arquivo==NULL){
        printf("Erro ao abrir o arquivo");
        return;
    }

    double variaveis[6];

    for(int i=0; i<6; i++){
        fscanf(arquivo, "%lf\n", &variaveis[i]);
    }

    double s = variaveis[0];
    double e = variaveis[1];
    double r = variaveis[2];
    double sigma = variaveis[3];
    double tempo = variaveis[4];
    int m = variaveis[5];
    double aux1, aux2, med, d_p, amplitude, min, max;
    double trials[m];
    int gseed = time(NULL);

    srand(gseed);
    for (int i = 0; i < m; i++)
    {
        aux1 = ((r - 0.5 * pow(sigma, 2)) * tempo + sigma * sqrt(tempo) * (((double)rand_r(&gseed)) / ((double)RAND_MAX)));
        aux2 = s * exp(aux1);
        if ((aux2 - e) < 0)
        {
            trials[i] = 0;
        }
        else
        {
            trials[i] = exp(-r * tempo) * (aux2 - e);
        }
        //printf("%lf \n", trials[i]);
    }

    med = media(trials, m);
    d_p = desvio_padrao(trials, m, med);
    //printf("%lf", d_p);
    amplitude = 1.96*d_p/(sqrt(m));
    //printf("%lf", amplitude);
    min = med-amplitude;
    max = med+amplitude;
    //printf("%lf, %lf", max, min);

    FILE *arquivoSaida = fopen("blackScholes_out.txt", "w");
    if(arquivoSaida == NULL){
        printf("Erro ao escrever no arquivo");
        return;
    }

    fprintf(arquivoSaida, "S = %lf \nE = %lf \nr = %lf \nsigma = %lf \nT = %lf \nM = %d \nConfidence interval: (%lf, %lf)", s, e, r, sigma, tempo, m, min, max);
}

int main()
{

    blackscholes();

    return 0;
}