#include "main.h"
#include <time.h>
// #include "mpi.h"
#include "par.h"

int xargc;
char **xargv;

int main(int argc, char *argv[])
{
    int i, j, k;
    int ntr, nt;
    int nshot;
    int nx, nz;
    float dx, dz;

    float alpha;
    double migSum;
    double seiSum;
    float migNum;
    float seiNum;
    float resMigNum;

    FILE *seisfp;
    FILE *migfp;

    FILE *resMigfp;
    FILE *newresMIGfp;

    char *seisname;
    char *migname;
    char *resMigname;
    char *newresMIGname;

    initargs(argc, argv);

    if (!getparstring("seisname", &seisname))
        seisname = "seis.dat";
    if (!getparstring("migname", &migname))
        migname = "mig.dat";
    if (!getparstring("resMigname", &resMigname))
        resMigname = "resMig.dat";
    if (!getparstring("newresMIGname", &newresMIGname))
        newresMIGname = "newresMIG.dat";

    if (!getparint("nt", &nt))
        nt = 750;
    if (!getparint("nshot", &nshot))
        nshot = 150;
    if (!getparint("ntr", &ntr))
        ntr = 100;
    if (!getparint("nz", &nz))
        nz = 750;
    if (!getparint("nx", &nx))
        nx = 737;

    seisfp = fopen(seisname, "rb");
    migfp = fopen(migname, "rb");

    seiNum = 0;
    migNum = 0;
    seiSum = 0;
    migSum = 0;

    for (i = 0; i < ntr * nshot; i++)
    {
        fread(&seiNum, sizeof(float), 1, seisfp);
        seiSum =seiSum+ fabs(seiNum);
        // printf("%f\n",(seiNum*1e10));

    }

    for (i = 0; i < nx * nz; i++)
    {
        fread(&migNum, sizeof(float), 1, migfp);
        migSum = migSum+ fabs(migNum);
    }

    alpha = migSum / seiSum;
    printf("migSum=%f,seisSum=%f,alpha=%f\n",migSum,seiSum,alpha);

    fclose(seisfp);
    fclose(migfp);

    migfp = fopen(migname, "rb");
    resMigfp = fopen(resMigname, "rb");
    newresMIGfp = fopen(newresMIGname, "wb");

    for (i = 0; i < nx * nz; i++)
    {
        fread(&migNum, sizeof(float), 1, migfp);
        fread(&resMigNum, sizeof(float), 1, resMigfp);
        seiNum = resMigNum - migNum * alpha;
        fwrite(&seiNum, sizeof(float), 1, newresMIGfp);
    }

    fclose(migfp);
    fclose(resMigfp);
    fclose(newresMIGfp);

    return 0;
}