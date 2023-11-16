


/* define rgeo using 3 points */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int defgrid_3p(float *point1, float *point2, float *point3, float *rgeo);

/* compute coordinates for a specific line,cmp number using rgeo */
void trans_lc2xy(float *rgeo, int *line, int *cmp, float *xcoor, float *ycoor);

/* compute line,cmp number for given coordinates using rgeo */
void trans_xy2lc(float *rgeo, float *xcoor, float *ycoor, int *line, int *cmp);
