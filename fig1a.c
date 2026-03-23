#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---------- parameters ---------- */

const double gam = 1.978e5;
const double tau   = 3.5e-9;
const double N0    = 0.175;
const double k0    = 0.17;
const double a     = 0.07;

/* ---------- bifurcation parameter range ---------- */

const double F_MIN = 184500.0;
const double F_MAX = 186000.0;
const int N_FREQ   = 500;

/* ---------- Integration ---------- */

const double DT      = 2e-8;
const double T_TOTAL = 0.2;

/* ---------- Maxima storage ---------- */

const long MAX_MAXIMA = 200000;   // upper bound for detected maxima

/* ---------- RHS ---------- */

void rhs(double t, double y[3], double dydt[3], double f)
{
    double I = y[0];
    double N = y[1];
    double z = y[2];

    dydt[0] = ((N - k0*(1 + a*cos(2*M_PI*f*z))) * I)/tau;
    dydt[1] = (N0 - N)*gam - I*N;
    dydt[2] = 1.0;
}

/* ---------- RK4 ---------- */

void rk4_step(double t, double y[3], double h, double f)
{
    double k1[3],k2[3],k3[3],k4[3],yt[3];

    rhs(t,y,k1,f);

    for(int i=0;i<3;i++) yt[i]=y[i]+0.5*h*k1[i];
    rhs(t+0.5*h,yt,k2,f);

    for(int i=0;i<3;i++) yt[i]=y[i]+0.5*h*k2[i];
    rhs(t+0.5*h,yt,k3,f);

    for(int i=0;i<3;i++) yt[i]=y[i]+h*k3[i];
    rhs(t+h,yt,k4,f);

    for(int i=0;i<3;i++)
        y[i]+= (h/6.0)*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
}

/* ---------- MAIN ---------- */

int main()
{
    FILE *fp = fopen("bif_laser_peaks(f=184.5_186).dat","w");
    if(!fp){ printf("Cannot open file\n"); return 1; }

    long N_steps = (long)(T_TOTAL/DT);

    double *I      = malloc(sizeof(double)*N_steps);
    double *Idot   = malloc(sizeof(double)*N_steps);
    double *maxbuf = malloc(sizeof(double)*MAX_MAXIMA);

    if(!I || !Idot || !maxbuf){
        printf("Memory allocation failed\n");
        return 1;
    }

    double N_S = k0;
    double I_S = gam*(N0/k0 - 1.0);

    for(int k = 0; k < N_FREQ; k++)
    {
        double f = F_MIN + k*(F_MAX - F_MIN)/(N_FREQ - 1);

        double y[3]={ I_S*(1+1e-4), N_S*(1+1e-4), 0.0 };
        double t=0;

        /* integrate */
        for(long i=0;i<N_steps;i++)
        {
            rk4_step(t,y,DT,f);
            I[i]=y[0];
            t+=DT;
        }

        /* numerical derivative */
        for(long i=1;i<N_steps-1;i++)
            Idot[i]=(I[i+1]-I[i-1])/(2*DT);

        /* collect ALL maxima */
        long n_max = 0;

        for(long i=2;i<N_steps-2;i++)
        {
            if(Idot[i-1]>0 && Idot[i+1]<0)
            {
                if(n_max < MAX_MAXIMA)
                {
                    maxbuf[n_max++] = I[i];
                }
            }
        }

        /* write last 10% maxima */
        long start_index = (long)(0.9 * n_max);

        for(long i = start_index; i < n_max; i++)
            fprintf(fp,"%lf %lf\n",f,maxbuf[i]);

        fprintf(fp,"\n");
    }

    fclose(fp);

    free(I);
    free(Idot);
    free(maxbuf);

    printf("Done.\n");
    return 0;
}




