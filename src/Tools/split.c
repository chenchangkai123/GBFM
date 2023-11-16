#include "split.h"
/* define rgeo using 3 points */
int defgrid_3p(float *point1, float *point2, float *point3, float *rgeo) {


	    int maxxline, maxinline, xlinedirection, flag;
	    float dsx, dsy, dxx, dxy, cdp1x, cdp1y;
	    int s1, s2, s3, c1, c2, c3;
	    int ss1, ss2, ss3, cc1, cc2, cc3;
	    float x1, x2, x3, y1, y2, y3;
	    float a, b, c, d;
	    int n; /* debug printout */

	    s1 = (int) point1[0];
	    c1 = (int) point1[1];
	    x1 = point1[2];
	    y1 = point1[3];

	    s2 = (int) point2[0];
	    c2 = (int) point2[1];
	    x2 = point2[2];
	    y2 = point2[3];

	    s3 = (int) point3[0];
	    c3 = (int) point3[1];
	    x3 = point3[2];
	    y3 = point3[3];

	    /********************************************************************
 		   x = cdpx1 + (s-1)*dsx + (c-1)*dxx
  		   y = cdpy1 + (s-1)*dsy + (c-1)*dxy
	    *********************************************************************/

	    a = (s2 - s1)*1.0;
	    b = (c2 - c1)*1.0;
	    c = (s3 - s2)*1.0;
	    d = (c3 - c2)*1.0;

	    if (a == 0.0 && c == 0.0) {
		printf("Kirdmig------message: Err!;  three points in same line \n");
		exit(0);
	    } else if (a == 0.0 && c != 0.0) {

		if (b == 0.0) {
		    printf("Kirdmig------message: Err!;   first 2 pionts are same \n");
		    exit(0);
		} else {
		    dxx = (x2 - x1) / b;
		    dsx = (x3 - x2 - d * dxx) / c;
		}

	    } else if (a != 0.0 && c == 0.0) {

		if (d == 0.0) {
		    printf("Kirdmig------message: Err!;    last 2 pionts are same \n");
		    exit(0);
		} else {
		    dxx = (x3 - x2) / d;
		    dsx = (x2 - x1 - b * dxx) / a;
		}

	    } else {
		if (b == 0.0 && d == 0.0) {
		    printf("Kirdmig------message: Err!;   Error message!  \n");
		    exit(0);
		} else {
		    dxx = ((x2 - x1) / a - (x3 - x2) / c) / (b / a - d / c);
		    dsx = (x3 - x2 - d * dxx) / c;
		}
	    }

	    if (a == 0.0 && c == 0.0) {
		printf("Kirdmig------message: Err!;  three points in same line \n");
		exit(0);
	    } else if (a == 0.0 && c != 0.0) {

		if (b == 0.0) {
		    printf("Kirdmig------message: Err!;   first 2 pionts are same \n");
		    exit(0);
		} else {
		    dxy = (y2 - y1) / b;
		    dsy = (y3 - y2 - d * dxy) / c;
		}

	    } else if (a != 0.0 && c == 0.0) {

		if (d == 0.0) {
		    printf("Kirdmig------message: Err!;    last 2 pionts are same \n");
		    exit(0);
		} else {
		    dxy = (y3 - y2) / d;
		    dsy = (y2 - y1 - b * dxy) / a;
		}

	    } else {
		if (b == 0.0 && d == 0.0) {
		    printf("Kirdmig------message: Err!;  Error message! \n");
		    exit(0);
		} else {
		    dxy = ((y2 - y1) / a - (y3 - y2) / c) / (b / a - d / c);
		    dsy = (y3 - y2 - d * dxy) / c;
		}
	    }

	    cdp1x = x1 - (s1 - 1) * dsx - (c1 - 1) * dxx;
	    cdp1y = y1 - (s1 - 1) * dsy - (c1 - 1) * dxy;

	    rgeo[0] = cdp1x;
	    rgeo[1] = cdp1y;
	    rgeo[2] = dsx;
	    rgeo[3] = dsy;
	    rgeo[4] = dxx;
	    rgeo[5] = dxy;


	    return 0;

	}

/* compute coordinates for a specific line,cmp number using rgeo */
void trans_lc2xy(float *rgeo, int *line, int *cmp, float *xcoor, float *ycoor)
{
	int ss, sx;
	ss = (*line-1);
	sx = (*cmp-1);
	*xcoor = rgeo[0] + ss*rgeo[2] + sx*rgeo[4];
	*ycoor = rgeo[1] + ss*rgeo[3] + sx*rgeo[5];
}

/* compute line,cmp number for given coordinates using rgeo */
void trans_xy2lc(float *rgeo, float *xcoor, float *ycoor, int *line, int *cmp)
{
	float dx, dy, dline, dcmp, lproj, cproj;

	dline = sqrt(rgeo[2]*rgeo[2] + rgeo[3]*rgeo[3]);
	dcmp = sqrt(rgeo[4]*rgeo[4] + rgeo[5]*rgeo[5]);

	dx = *xcoor-rgeo[0];
	dy = *ycoor-rgeo[1];

	lproj = (dx*rgeo[2] + dy*rgeo[3])/dline;
	cproj = (dx*rgeo[4] + dy*rgeo[5])/dcmp;

	*line = 1 + lproj/dline + 0.5;
	*cmp = 1 + cproj/dcmp + 0.5;
}
