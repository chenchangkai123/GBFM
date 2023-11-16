#include <stdio.h>
#include <stdlib.h>
#include "par.h"
#include "SegyRead.h"
#include "main.h"

int xargc;
char **xargv;

int main(int argc, char **argv)
{
    int i, j, k, ir, ii;
    int nx, ny, nz;
    char *input;
    char *output;
    char *elev;

    char *headRecx;
    char *headShotx;
    char *headRecy;
    char *headShoty;

    float datum;
    int isdepth;

    float dz, scalco, scalel;
    float tempv, tempelev;
    float recx, recy, shotx, shoty;
    int ntr, ns, ntrmax;
    int **dataTemp;
    int **dataTempTo;
    int ibmieee;

    FILE *infp;
    FILE *outfp;
    FILE *elevfp;

    FILE *headRecxfp;
    FILE *headShotxfp;
    FILE *headRecyfp;
    FILE *headShotyfp;

    // float **vdata;
    float elevdata;

    segy *tr;
    tr = (segy *)malloc(sizeof(segy));

    initargs(argc, argv);
    if (!getparstring("input", &input))
        input = "shot.segy";
    if (!getparstring("output", &output))
        output = "shot.cut";
    if (!getparstring("elev", &elev))
        elev = "elev.dat";
    if (!getparfloat("datum", &datum))
        datum = 0.0;
    if (!getparint("isdepth", &isdepth))
        isdepth = 0;

    if (!getparstring("headRecx", &headRecx))
        headRecx = "headRecx.dat";
    if (!getparstring("headShotx", &headShotx))
        headShotx = "headShotx.dat";
    if (!getparstring("headRecy", &headRecy))
        headRecy = "headRecy.dat";
    if (!getparstring("headShoty", &headShoty))
        headShoty = "headShoty.dat";
    if (!getparint("ibmieee", &ibmieee))
        ibmieee = 0;

    infp = fopen(input, "rb");
    outfp = fopen(output, "wb");
    elevfp = fopen(elev, "wb");

    headRecxfp = fopen(headRecx, "wb");
    headShotxfp = fopen(headShotx, "wb");
    headRecyfp = fopen(headRecy, "wb");
    headShotyfp = fopen(headShoty, "wb");

    if (infp == NULL || outfp == NULL || elevfp == NULL ||
        headRecxfp == NULL || headShotxfp == NULL ||
        headRecyfp == NULL || headShotyfp == NULL)
    {
        fprintf(stderr, "Error opening file.\n");
        exit(1);
    }

    fseek(infp, 3600L, SEEK_SET);

    fread(tr, 240, 1, infp);
    scalel = (double)swap_int16(tr->scalel);
    scalco = (double)swap_int16(tr->scalco);
    if (scalco != 0)
    {
        scalco = scalco > 0 ? scalco : 1 / fabs(scalco);
    }
    else
    {
        scalco = 1;
    }
    if (scalel != 0)
    {
        scalel = scalel > 0 ? scalel : 1 / fabs(scalel);
    }
    else
    {
        scalel = 1;
    }

    nz = swap_uint16(tr->ns);
    ntr = 1;

    printf("fx_min=%d,fy_min=%d\n", swap_uint32(tr->gx), swap_uint32(tr->gy));

    fseek(infp, (nz) * sizeof(float), SEEK_CUR);
    while (fread(tr, 240, 1, infp) == 1)
    {
        ntr++;
        fseek(infp, (nz) * sizeof(float), SEEK_CUR);
    }

    // vdata = alloc2float(nz, ntr);
    // elevdata = alloc1float(ntr);
    dataTemp = alloc2int(nz, ntr);
    dataTempTo = alloc2int(nz, ntr);
    fclose(infp);

    infp = fopen(input, "rb");
    fseek(infp, 3600L, SEEK_SET);
    for (ir = 0; ir < ntr; ir++)
    {
        fread(tr, 240, 1, infp);
        scalel = (double)swap_int16(tr->scalel);
        scalco = (double)swap_int16(tr->scalco);
        if (scalco != 0)
        {
            scalco = scalco > 0 ? scalco : 1 / fabs(scalco);
        }
        else
        {
            scalco = 1;
        }
        if (scalel != 0)
        {
            scalel = scalel > 0 ? scalel : 1 / fabs(scalel);
        }
        else
        {
            scalel = 1;
        }
        if (isdepth == 0)
        {
            // elevdata[ir] = swap_uint32(tr->gelev) * scalel;
            elevdata = swap_uint32(tr->gelev) * scalel;
        }
        else
        {
            // elevdata[ir] = datum - swap_uint32(tr->gelev) * scalel;
            elevdata = datum - swap_uint32(tr->gelev) * scalel;
        }
        recx = swap_uint32(tr->gx) * scalco;
        recy = swap_uint32(tr->gy) * scalco;
        shotx = swap_uint32(tr->sx) * scalco;
        shoty = swap_uint32(tr->sy) * scalco;

        fwrite(&elevdata, sizeof(float), 1, elevfp);
        fwrite(&recx, sizeof(float), 1, headRecxfp);
        fwrite(&recy, sizeof(float), 1, headRecyfp);
        fwrite(&shotx, sizeof(float), 1, headShotxfp);
        fwrite(&shoty, sizeof(float), 1, headShotyfp);

        // printf("elevdata[%d]=%f\n", ir, elevdata[ir]);
        fread(&dataTemp[ir][0], sizeof(int), nz, infp);
        // 字节反序
        for (ii = 0; ii < nz; ii++)
            dataTemp[ir][ii] = swap_uint32(dataTemp[ir][ii]);
        
        if (ibmieee == 0)
        {
            ibm2ieee(dataTemp[ir], dataTempTo[ir], nz, 1);
        }
        else
        {
            ieee2ibm(dataTemp[ir], dataTempTo[ir], nz, 1);
        }
        
        // memcpy(&vdata[ir][0], dataTempTo[ir], sizeof(float) * nz);
        fwrite(&dataTempTo[ir][0], sizeof(float), nz, outfp);
    }

    // save vdata and elevdata
    // fwrite(&vdata[0][0], sizeof(float), ntr * nz, outfp);
    // fwrite(&elevdata[0], sizeof(float), ntr, elevfp);

    fclose(infp);
    fclose(outfp);
    fclose(elevfp);
    fclose(headRecxfp);
    fclose(headShotxfp);
    fclose(headRecyfp);
    fclose(headShotyfp);

    // printf("elevdata[0]=%f\n", elevdata[0]);

    printf("ntr=%d, nz=%d\n", ntr, nz);
    printf("job done!\n");

    return (0);
}
