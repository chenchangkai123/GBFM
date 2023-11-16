#include <stdio.h>
#include <stdlib.h>
#include "par.h"
#include "SegyRead.h"
#include "main.h"

int xargc;
char **xargv;

int main(int argc, char **argv)
{
    int i,j,k;
    float velocity;
    char *output;
    int nx, nz;

    FILE *outfp;

    initargs(argc, argv);
    if (!getparstring("output", &output))
        output = "model.dat";
    if (!getparint("nx", &nx))
        nx = 100;
    if (!getparint("nz", &nz))
        nz = 100;
    if (!getparfloat("velocity", &velocity))
        velocity = 1500.0;

    outfp=fopen(output, "wb");
    for (i=0; i<nz; i++)
    {
        for (j=0; j<nx; j++)
        {
            fwrite(&velocity, sizeof(float), 1, outfp);
        }
    }
    fclose(outfp);
    printf("nx=%d, nz=%d, velocity=%f\n", nx, nz, velocity);
    printf("Model file %s generated.\n", output);

    return 0;
}