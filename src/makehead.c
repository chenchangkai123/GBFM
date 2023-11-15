#include "main.h"
#include "par.h"

int xargc;
char **xargv;

int main(int argc, char *argv[])
{
    int i, j, k;

    int ntr;
    float dtr;
    float ftr;
    float dftr;
    int totaltr;

    int nshot;
    float dshot;
    float fshot;

    float tracesShot, traceRec;
    FILE *headRecfp;
    FILE *headShotfp;
    FILE *headInfo;

    initargs(argc, argv);

    if (!getparint("ntr", &ntr))
        ntr = 100;
    if (!getparfloat("dtr", &dtr))
        dtr = 12.5;
    if (!getparfloat("ftr", &ftr))
        ftr = 0.0;
    if (!getparfloat("dftr", &dftr))
        dftr = 50;

    if (!getparint("nshot", &nshot))
        nshot = 150;
    if (!getparfloat("dshot", &dshot))
        dshot = 50;
    if (!getparfloat("fshot", &fshot))
        fshot = 50.0;

    headRecfp = fopen("headRec.dat", "wb");
    headShotfp = fopen("headShot.dat", "wb");
    for (i = 0; i < nshot; i++)
    {
        for (j = 0; j < ntr; j++)
        {
            traceRec = ftr + dftr * i + j * dtr;
            tracesShot = fshot + dshot * i;
            fwrite(&tracesShot, sizeof(float), 1, headShotfp);
            fwrite(&traceRec, sizeof(float), 1, headRecfp);
        }
    }
    fclose(headRecfp);
    fclose(headShotfp);
    k=0;
    headInfo = fopen("headInfo", "w");
    for(i=0;i<nshot;i++)
    {
        fprintf(headInfo,"%d,%f,%f,%d,%d\n",i+1,fshot+dshot*i,ftr+i*dftr,ntr,k);
        k=k+ntr;
        // ishot, shotposition, firsttrace, ntraces, totaltraces
    }
    fclose(headInfo);

    printf("\033[0;32mHead files created\033[0m\n");
    printf("\033[0;32mInfo:\033[0m\n");
    printf("ntr=%d, dtr=%f, ftr=%f, dftr=%f\n", ntr, dtr, ftr, dftr);
    printf("nshot=%d, dshot=%f, fshot=%f\n", nshot, dshot, fshot);
    printf("maxtraces=%d\n", ntr*nshot);
    printf("\033[0;32mFile names:\033[0m\n");
    printf("headRec.dat\n");
    printf("headShot.dat\n");

    return 0;
}