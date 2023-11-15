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
    FILE *originShotfp;
    FILE *newShotfp;
    FILE *subtractDatafp;
    char *originShotname;
    char *newShotname;
    char *subtractDataname;
    float originShot;
    float newShot;
    float subtractData;

    initargs(argc, argv);

    if (!getparstring("originShotname", &originShotname))
        originShotname = "shot.dat";
    if (!getparstring("newShotname", &newShotname))
        newShotname = "newShot.dat";
    if (!getparstring("subtractDataname", &subtractDataname))
        subtractDataname = "subtractData.dat";
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

    originShotfp = fopen(originShotname, "rb");
    newShotfp = fopen(newShotname, "rb");
    subtractDatafp = fopen(subtractDataname, "wb");

    for (i = 0; i < ntr * nshot; i++)
    {
        fread(&originShot, sizeof(float), 1, originShotfp);
        fread(&newShot, sizeof(float), 1, newShotfp);
        subtractData = newShot-originShot;
        fwrite(&subtractData, sizeof(float), 1, subtractDatafp);
    }

    fclose(originShotfp);
    fclose(newShotfp);
    fclose(subtractDatafp);

    return 0;
}