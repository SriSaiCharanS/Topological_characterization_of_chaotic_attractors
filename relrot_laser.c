#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.14159265358979323846

#define MAXPTS 200000
#define MAXP   100

#define ARG(x,y) atan2((double)(y),(double)(x))
//#error COMPILING THIS FILE

int main(int argc, char *argv[])
{
    printf("Program started\n");

    if (argc != 6) {
        printf("Usage: %s orbitA orbitB PA PB T\n", argv[0]);
        return 1;
    }

    const char *fileA = argv[1];
    const char *fileB = argv[2];
    int PA = atoi(argv[3]);
    int PB = atoi(argv[4]);
    double T = atof(argv[5]);

    FILE *fa = fopen(fileA,"r");
    FILE *fb = fopen(fileB,"r");

    if(!fa || !fb){
        printf("File open error\n");
        return 1;
    }

    float A[3][MAXPTS], B[3][MAXPTS];
    float rx[MAXPTS], ry[MAXPTS];

    int I[MAXP], J[MAXP];
    float RR[MAXP][MAXP] = {{0.0f}};

    double t; /*dummy;*/

    int n;
    int countA = 0, countB = 0;

    /* =========================
       READ ORBIT A
       ========================= */

    double prev_phase = -1;

    for(n=0; n<MAXPTS; ++n){

        if(fscanf(fa,"%lf %f %f",
                  &t,&A[2][n],&A[1][n])!=3)
            break;

        double phase = fmod(t,T);
        if(phase < 0) phase += T;

        A[0][n] = phase/T;

        /* robust Poincare detection:
           detect wraparound of phase */
        if(prev_phase >= 0 && phase < prev_phase){
            if(countA < MAXP)
                I[countA++] = n;
        }

        prev_phase = phase;
    }

    int M = n;
    fclose(fa);

    /* =========================
       READ ORBIT B
       ========================= */

    prev_phase = -1;

    for(n=0; n<MAXPTS; ++n){

        if(fscanf(fb,"%lf %f %f",
                  &t,&B[2][n],&B[1][n])!=3)
            break;

        double phase = fmod(t,T);
        if(phase < 0) phase += T;

        B[0][n] = phase/T;

        if(prev_phase >= 0 && phase < prev_phase){
            if(countB < MAXP)
                J[countB++] = n;
        }

        prev_phase = phase;
    }

    int N = n;
    fclose(fb);

    printf("Orbit A points: %d   section hits: %d\n",M,countA);
    printf("Orbit B points: %d   section hits: %d\n",N,countB);

    if(countA < PA || countB < PB){
        printf("Not enough Poincare hits detected!\n");
        return 1;
    }

    /* =========================
       CORE ALGORITHM
       ========================= */

    int Q = (M < N ? M : N);

    for(int i=0;i<PA;i++){
        for(int j=0;j<PB;j++){

            int mi = I[i];
            int nj = J[j];

            for(int q=0;q<Q;q++){

                rx[q] = B[1][nj] - A[1][mi];
                ry[q] = B[2][nj] - A[2][mi];

                mi = (mi+1) % M;
                nj = (nj+1) % N;
            }

            for(int q=0;q<Q-1;q++){

                double a0 = ARG(rx[q],ry[q]);
                double a1 = ARG(rx[q+1],ry[q+1]);

                double da = a1 - a0;

                /* unwrap angle */
                if(da > PI)  da -= 2*PI;
                if(da < -PI) da += 2*PI;

                RR[i][j] += (float)da;
            }
        }
    }

    /* =========================
       OUTPUT
       ========================= */

    printf("\nROTATION RATES\n\n");
    printf("PA: %d PB: %d\n\n",PA,PB);

    for(int i=0;i<PA;i++){
        for(int j=0;j<PB;j++){

            printf("Index: %d,%d  Rot:%g  RelRot:%g\n",
                i+1,j+1,
                RR[i][j]/(2*PI),
                RR[i][j]/(2*PI*PA*PB));
        }
    }

    return 0;
}

