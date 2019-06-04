/* Program to compute Pi using Monte Carlo methods */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <pthread.h>

unsigned int gseed = 45162;

int montecarlo()
{
    int niter = 1000000000;
    double x, y, z;
    int i, count = 0; /* # of points in the 1st quadrant of unit circle */
    double pi;

    /* initialize random numbers */
    for (i = 0; i < niter; i++)
    {
        x = (((double)rand_r(&gseed)) / ((double)RAND_MAX));
        y = (((double)rand_r(&gseed)) / ((double)RAND_MAX));
        z = x * x + y * y;
        if (z <= 1)
            count++;
    }
    pi = (double)count / niter * 4;
    printf("# of trials= %d , estimate of pi is %lf \n", niter, pi);
    return 1;
}

int main(int argc, char *argv)
{
    gseed = time(NULL);
    srand(gseed);
    montecarlo();

    //rand_number();
    //printf("%u\n", aux);
    return 0;
}