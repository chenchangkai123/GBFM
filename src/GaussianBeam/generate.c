#include "main.h"
#include <time.h>
// #include "mpi.h"
#include "par.h"

int xargc;
char **xargv;

int main(int argc, char *argv[])
{
    int i,j,k;
    int nx,nz;
    float num;

    FILE *fp;
    char *name;

    initargs(argc, argv);

    if (!getparstring("name", &name))
        name = "test.dat";
    if (!getparint("nx", &nx))
        nx = 737;
    if (!getparint("nz", &nz))
        nz = 750;
    if (!getparfloat("num", &num))
        num = 0;

    fp = fopen(name, "wb");

    for (i = 0; i < nx * nz; i++)
    {
        fwrite(&num, sizeof(float), 1, fp);
    }

    fclose(fp);

    return 0;
}