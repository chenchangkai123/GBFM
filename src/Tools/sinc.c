/* Copyright (c) Colorado School of Mines, 1999.*/
/* All rights reserved.                       */

/*********************** self documentation **********************/
/*****************************************************************************
SINC - Return SINC(x) for as floats or as doubles

fsinc		return float value of sinc(x) for x input as a float
dsinc		return double precision sinc(x) for double precision x

******************************************************************************
Function Prototype:
float fsinc (float x);
double dsinc (double x);

******************************************************************************
Input:
x		value at which to evaluate sinc(x)

Returned: 	sinc(x)

******************************************************************************
Notes:
    sinc(x) = sin(PI*x)/(PI*x) 

******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/02/89
*****************************************************************************/
/**************** end self doc ********************************/

#include "common.h"



float fsinc (float x)
/*****************************************************************************
Return sinc(x) = sin(PI*x)/(PI*x) (float version)
******************************************************************************
Input:
x		value at which to evaluate sinc(x)

Returned: 	sinc(x)
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/02/89
*****************************************************************************/
{
	float pix;

	if (x==0.0) {
		return 1.0;
	} else {
		pix = PI*x;
		return sin(pix)/pix;
	}
}

double dsinc (double x)
/*****************************************************************************
Return sinc(x) = sin(PI*x)/(PI*x) (double version)
******************************************************************************
Input:
x		value at which to evaluate sinc(x)

Returned:	sinc(x)
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/02/89
*****************************************************************************/
{
	double pix;

	if (x==0.0) {
		return 1.0;
	} else {
		pix = PI*x;
		return sin(pix)/pix;
	}
}
/*********************** self documentation **********************/
/*****************************************************************************
MKSINC - Compute least-squares optimal sinc interpolation coefficients.

mksinc		Compute least-squares optimal sinc interpolation coefficients.

******************************************************************************
Function Prototype:
void mksinc (float d, int lsinc, float sinc[]);

******************************************************************************
Input:
d		fractional distance to interpolation point; 0.0<=d<=1.0
lsinc		length of sinc approximation; lsinc%2==0 and lsinc<=20

Output:
sinc		array[lsinc] containing interpolation coefficients

******************************************************************************
Notes:
The coefficients are a least-squares-best approximation to the ideal
sinc function for frequencies from zero up to a computed maximum
frequency.  For a given interpolator length, lsinc, mksinc computes
the maximum frequency, fmax (expressed as a fraction of the nyquist
frequency), using the following empirically derived relation (from
a Western Geophysical Technical Memorandum by Ken Larner):

	fmax = min(0.066+0.265*log(lsinc),1.0)

Note that fmax increases as lsinc increases, up to a maximum of 1.0.
Use the coefficients to interpolate a uniformly-sampled function y(i) 
as follows:

            lsinc-1
    y(i+d) =  sum  sinc[j]*y(i+j+1-lsinc/2)
              j=0

Interpolation error is greatest for d=0.5, but for frequencies less
than fmax, the error should be less than 1.0 percent.

******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/02/89
*****************************************************************************/
/**************** end self doc ********************************/

void mksinc (float d, int lsinc, float sinc[])
/*****************************************************************************
Compute least-squares optimal sinc interpolation coefficients.
******************************************************************************
Input:
d		fractional distance to interpolation point; 0.0<=d<=1.0
lsinc		length of sinc approximation; lsinc%2==0 and lsinc<=20

Output:
sinc		array[lsinc] containing interpolation coefficients
******************************************************************************
Notes:
The coefficients are a least-squares-best approximation to the ideal
sinc function for frequencies from zero up to a computed maximum
frequency.  For a given interpolator length, lsinc, mksinc computes
the maximum frequency, fmax (expressed as a fraction of the nyquist
frequency), using the following empirically derived relation (from
a Western Geophysical Technical Memorandum by Ken Larner):

	fmax = min(0.066+0.265*log(lsinc),1.0)

Note that fmax increases as lsinc increases, up to a maximum of 1.0.
Use the coefficients to interpolate a uniformly-sampled function y(i) 
as follows:

            lsinc-1
    y(i+d) =  sum  sinc[j]*y(i+j+1-lsinc/2)
              j=0

Interpolation error is greatest for d=0.5, but for frequencies less
than fmax, the error should be less than 1.0 percent.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/02/89
*****************************************************************************/
{
	int j;
	double s[20],a[20],c[20],work[20],fmax;

	/* compute auto-correlation and cross-correlation arrays */
	fmax = 0.066+0.265*log((double)lsinc);
	fmax = (fmax<1.0)?fmax:1.0;
	for (j=0; j<lsinc; j++) {
		a[j] = dsinc(fmax*j);
		c[j] = dsinc(fmax*(lsinc/2-j-1+d));
	}

	/* solve symmetric Toeplitz system for the sinc approximation */
	stoepd(lsinc,a,c,s,work);
	for (j=0; j<lsinc; j++)
		sinc[j] = s[j];
}


/* these are used by both ints8c and ints8r */
#define LTTABLE 8
#define NNTABLE 513

void ints8c (int nxin, float dxin, float fxin, complex yin[], 
	complex yinl, complex yinr, int nxout, float xout[], complex yout[])
/*****************************************************************************
Interpolation of a uniformly-sampled complex function y(x) via a
table of 8-coefficient sinc approximations; maximum error for frequiencies
less than 0.6 nyquist is less than one percent.
******************************************************************************
Input:
nxin		number of x values at which y(x) is input
dxin		x sampling interval for input y(x)
fxin		x value of first sample input
yin		array[nxin] of input y(x) values:  yin[0] = y(fxin), etc.
yinl		value used to extrapolate yin values to left of yin[0]
yinr		value used to extrapolate yin values to right of yin[nxin-1]
nxout		number of x values a which y(x) is output
xout		array[nxout] of x values at which y(x) is output

Output:
yout		array[nxout] of output y(x):  yout[0] = y(xout[0]), etc.
******************************************************************************
Notes:
Because extrapolation of the input function y(x) is defined by the
left and right values yinl and yinr, the xout values are not restricted
to lie within the range of sample locations defined by nxin, dxin, and
fxin.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/02/89
*****************************************************************************/
{
	static float table[NNTABLE][LTTABLE];
	static int tabled=0;
	int jtable;
	float frac;

	/* tabulate sinc interpolation coefficients if not already tabulated */
	if (!tabled) {
		for (jtable=1; jtable<NNTABLE-1; jtable++) {
			frac = (float)jtable/(float)(NNTABLE-1);
			mksinc(frac,LTTABLE,&table[jtable][0]);
		}
		for (jtable=0; jtable<LTTABLE; jtable++) {
			table[0][jtable] = 0.0;
			table[NNTABLE-1][jtable] = 0.0;
		}
		table[0][LTTABLE/2-1] = 1.0;
		table[NNTABLE-1][LTTABLE/2] = 1.0;
		tabled = 1;
	}

	/* interpolate using tabulated coefficients */
	intt8c(NNTABLE,table,nxin,dxin,fxin,yin,yinl,yinr,nxout,xout,yout);
}

void ints8r (int nxin, float dxin, float fxin, float yin[], 
	float yinl, float yinr, int nxout, float xout[], float yout[])
/*****************************************************************************
Interpolation of a uniformly-sampled real function y(x) via a
table of 8-coefficient sinc approximations; maximum error for frequiencies
less than 0.6 nyquist is less than one percent.
******************************************************************************
Input:
nxin		number of x values at which y(x) is input
dxin		x sampling interval for input y(x)
fxin		x value of first sample input
yin		array[nxin] of input y(x) values:  yin[0] = y(fxin), etc.
yinl		value used to extrapolate yin values to left of yin[0]
yinr		value used to extrapolate yin values to right of yin[nxin-1]
nxout		number of x values a which y(x) is output
xout		array[nxout] of x values at which y(x) is output

Output:
yout		array[nxout] of output y(x):  yout[0] = y(xout[0]), etc.
******************************************************************************
Notes:
Because extrapolation of the input function y(x) is defined by the
left and right values yinl and yinr, the xout values are not restricted
to lie within the range of sample locations defined by nxin, dxin, and
fxin.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/02/89
*****************************************************************************/

{
	static float table[NNTABLE][LTTABLE];
	static int tabled=0;
	int jtable;
	float frac;

	/* tabulate sinc interpolation coefficients if not already tabulated */
	if (!tabled) {
		for (jtable=1; jtable<NNTABLE-1; jtable++) {
			frac = (float)jtable/(float)(NNTABLE-1);
			mksinc(frac,LTTABLE,&table[jtable][0]);
		}
		for (jtable=0; jtable<LTTABLE; jtable++) {
			table[0][jtable] = 0.0;
			table[NNTABLE-1][jtable] = 0.0;
		}
		table[0][LTTABLE/2-1] = 1.0;
		table[NNTABLE-1][LTTABLE/2] = 1.0;
		tabled = 1;
	}

	/* interpolate using tabulated coefficients */
	intt8r(NNTABLE,table,nxin,dxin,fxin,yin,yinl,yinr,nxout,xout,yout);
}

/*****************************************************************************
INTTABLE8 -  Interpolation of a uniformly-sampled complex function y(x)
		via a table of 8-coefficient interpolators

intt8c	interpolation of a uniformly-sampled complex function y(x)
		via a table of 8-coefficient interpolators
intt8r	interpolation of a uniformly-sampled real function y(x) via a
		table of 8-coefficient interpolators

******************************************************************************
Function Prototype:
void intt8c (int ntable, float table[][8],
	int nxin, float dxin, float fxin, complex yin[], 
	complex yinl, complex yinr, int nxout, float xout[], complex yout[]);
void intt8r (int ntable, float table[][8],
	int nxin, float dxin, float fxin, float yin[], 
	float yinl, float yinr, int nxout, float xout[], float yout[]);

******************************************************************************
Input:
ntable		number of tabulated interpolation operators; ntable>=2
table		array of tabulated 8-point interpolation operators
nxin		number of x values at which y(x) is input
dxin		x sampling interval for input y(x)
fxin		x value of first sample input
yin		array of input y(x) values:  yin[0] = y(fxin), etc.
yinl		value used to extrapolate yin values to left of yin[0]
yinr		value used to extrapolate yin values to right of yin[nxin-1]
nxout		number of x values a which y(x) is output
xout		array of x values at which y(x) is output

Output:
yout		array of output y(x) values:  yout[0] = y(xout[0]), etc.

******************************************************************************
NOTES:
ntable must not be less than 2.

The table of interpolation operators must be as follows:

Let d be the distance, expressed as a fraction of dxin, from a particular
xout value to the sampled location xin just to the left of xout.  Then,
for d = 0.0,

table[0][0:7] = 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0

are the weights applied to the 8 input samples nearest xout.
Likewise, for d = 1.0,

table[ntable-1][0:7] = 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0

are the weights applied to the 8 input samples nearest xout.  In general,
for d = (float)itable/(float)(ntable-1), table[itable][0:7] are the
weights applied to the 8 input samples nearest xout.  If the actual sample
distance d does not exactly equal one of the values for which interpolators
are tabulated, then the interpolator corresponding to the nearest value of
d is used.

Because extrapolation of the input function y(x) is defined by the left
and right values yinl and yinr, the xout values are not restricted to lie
within the range of sample locations defined by nxin, dxin, and fxin.

******************************************************************************
AUTHOR:  Dave Hale, Colorado School of Mines, 06/02/89
*****************************************************************************/
/**************** end self doc ********************************/


void intt8c (int ntable, float table[][8],
	int nxin, float dxin, float fxin, complex yin[],
	complex yinl, complex yinr, int nxout, float xout[], complex yout[])
{
	int ioutb,ixout,ixoutn,kyin,ktable,itable,nxinm8;
	float xoutb,xoutf,xouts,xoutn,fntablem1,frac,
		yinlr,yinli,yinrr,yinri,yinir,yinii,sumr,sumi,
		*table00,*pyin,*ptable;
	complex *yin0;

	/* compute constants */
	ioutb = -3-8;
	xoutf = fxin;
	xouts = 1.0/dxin;
	xoutb = 8.0-xoutf*xouts;
	fntablem1 = (float)(ntable-1);
	nxinm8 = nxin-8;
	yin0 = &yin[0];
	table00 = &table[0][0];
	yinlr = yinl.r;  yinli = yinl.i;
	yinrr = yinr.r;  yinri = yinr.i;

	/* loop over output samples */
	for (ixout=0; ixout<nxout; ixout++) {

		/* determine pointers into table and yin */
		xoutn = xoutb+xout[ixout]*xouts;
		ixoutn = (int)xoutn;
		kyin = ioutb+ixoutn;
		pyin = (float*)(yin0+kyin);
		frac = xoutn-(float)ixoutn;
		ktable = frac>=0.0?frac*fntablem1+0.5:(frac+1.0)*fntablem1-0.5;
		ptable = table00+ktable*8;
		
		/* if totally within input array, use fast method */
		if (kyin>=0 && kyin<=nxinm8) {
			yout[ixout].r = 
				pyin[0]*ptable[0]+
				pyin[2]*ptable[1]+
				pyin[4]*ptable[2]+
				pyin[6]*ptable[3]+
				pyin[8]*ptable[4]+
				pyin[10]*ptable[5]+
				pyin[12]*ptable[6]+
				pyin[14]*ptable[7];
			yout[ixout].i = 
				pyin[1]*ptable[0]+
				pyin[3]*ptable[1]+
				pyin[5]*ptable[2]+
				pyin[7]*ptable[3]+
				pyin[9]*ptable[4]+
				pyin[11]*ptable[5]+
				pyin[13]*ptable[6]+
				pyin[15]*ptable[7];
		
		/* else handle end effects with care */
		} else {

			sumr = sumi = 0.0;
			for (itable=0; itable<8; itable++,kyin++) {
				if (kyin<0) {
					yinir = yinlr;
					yinii = yinli;
				} else if (kyin>=nxin) {
					yinir = yinrr;
					yinii = yinri;
				} else {
					yinir = yin[kyin].r;
					yinii = yin[kyin].i;
				}
				sumr += yinir*(*ptable);
				sumi += yinii*(*ptable++);
			}
			yout[ixout].r = sumr;
			yout[ixout].i = sumi;
		}
	}
}


void intt8r (int ntable, float table[][8],
	int nxin, float dxin, float fxin, float yin[], float yinl, float yinr,
	int nxout, float xout[], float yout[])
{
	int ioutb,nxinm8,ixout,ixoutn,kyin,ktable,itable;
	float xoutb,xoutf,xouts,xoutn,frac,fntablem1,yini,sum,
		*yin0,*table00,*pyin,*ptable;

	/* compute constants */
	ioutb = -3-8;
	xoutf = fxin;
	xouts = 1.0/dxin;
	xoutb = 8.0-xoutf*xouts;
	fntablem1 = (float)(ntable-1);
	nxinm8 = nxin-8;
	yin0 = &yin[0];
	table00 = &table[0][0];

	/* loop over output samples */
	for (ixout=0; ixout<nxout; ixout++) {

		/* determine pointers into table and yin */
		xoutn = xoutb+xout[ixout]*xouts;
		ixoutn = (int)xoutn;
		kyin = ioutb+ixoutn;
		pyin = yin0+kyin;
		frac = xoutn-(float)ixoutn;
		ktable = frac>=0.0?frac*fntablem1+0.5:(frac+1.0)*fntablem1-0.5;
		ptable = table00+ktable*8;
		
		/* if totally within input array, use fast method */
		if (kyin>=0 && kyin<=nxinm8) {
			yout[ixout] = 
				pyin[0]*ptable[0]+
				pyin[1]*ptable[1]+
				pyin[2]*ptable[2]+
				pyin[3]*ptable[3]+
				pyin[4]*ptable[4]+
				pyin[5]*ptable[5]+
				pyin[6]*ptable[6]+
				pyin[7]*ptable[7];
		
		/* else handle end effects with care */
		} else {
	
			/* sum over 8 tabulated coefficients */
			for (itable=0,sum=0.0; itable<8; itable++,kyin++) {
				if (kyin<0)
					yini = yinl;
				else if (kyin>=nxin)
					yini = yinr;
				else
					yini = yin[kyin];
				sum += yini*(*ptable++);
			}
			yout[ixout] = sum;
		}
	}
}
/*****************************************************************************
STOEP - Functions to solve a symmetric Toeplitz linear system of equations
	 Rf=g for f

stoepd		solve a symmetric Toeplitz system - doubles
stoepf		solve a symmetric Toeplitz system - floats

******************************************************************************
Function Prototypes:
void stoepd (int n, double r[], double g[], double f[], double a[]);
void stoepf (int n, float r[], float g[], float f[], float a[]);

******************************************************************************
Input:
n		dimension of system
r		array[n] of top row of Toeplitz matrix
g		array[n] of right-hand-side column vector

Output:
f		array[n] of solution (left-hand-side) column vector
a		array[n] of solution to Ra=v (Claerbout, FGDP, p. 57)

******************************************************************************
Notes:
These routines do NOT solve the case when the main diagonal is zero, it
just silently returns.

The left column of the Toeplitz matrix is assumed to be equal to the top
row (as specified in r); i.e., the Toeplitz matrix is assumed symmetric.

******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/02/89
*****************************************************************************/
/**************** end self doc ********************************/
void stoepd (int n, double r[], double g[], double f[], double a[])
/*****************************************************************************
Solve a symmetric Toeplitz linear system of equations Rf=g for f
(double version)
******************************************************************************
Input:
n		dimension of system
r		array[n] of top row of Toeplitz matrix
g		array[n] of right-hand-side column vector

Output:
f		array[n] of solution (left-hand-side) column vector
a		array[n] of solution to Ra=v (Claerbout, FGDP, p. 57)
******************************************************************************
Notes:
This routine does NOT solve the case when the main diagonal is zero, it
just silently returns.

The left column of the Toeplitz matrix is assumed to be equal to the top
row (as specified in r); i.e., the Toeplitz matrix is assumed symmetric.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/02/89
*****************************************************************************/
{
	int i,j;
	double v,e,c,w,bot;

	if (r[0] == 0.0) return;

	a[0] = 1.0;
	v = r[0];
	f[0] = g[0]/r[0];

	for (j=1; j<n; j++) {
		
		/* solve Ra=v as in Claerbout, FGDP, p. 57 */
		a[j] = 0.0;
		f[j] = 0.0;
		for (i=0,e=0.0; i<j; i++)
			e += a[i]*r[j-i];
		c = e/v;
		v -= c*e;
		for (i=0; i<=j/2; i++) {
			bot = a[j-i]-c*a[i];
			a[i] -= c*a[j-i];
			a[j-i] = bot;
		}

		/* use a and v above to get f[i], i = 0,1,2,...,j */
		for (i=0,w=0.0; i<j; i++)
			w += f[i]*r[j-i];
		c = (w-g[j])/v;
		for (i=0; i<=j; i++)
			f[i] -= c*a[j-i];
	}
}

void stoepf (int n, float r[], float g[], float f[], float a[])
/*****************************************************************************
Solve a symmetric Toeplitz linear system of equations Rf=g for f
(float version) 
******************************************************************************
Input:
n		dimension of system
r		array[n] of top row of Toeplitz matrix
g		array[n] of right-hand-side column vector

Output:
f		array[n] of solution (left-hand-side) column vector
a		array[n] of solution to Ra=v (Claerbout, FGDP, p. 57)
******************************************************************************
Notes:
This routine does NOT solve the case when the main diagonal is zero, it
just silently returns.

The left column of the Toeplitz matrix is assumed to be equal to the top
row (as specified in r); i.e., the Toeplitz matrix is assumed symmetric.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/02/89
*****************************************************************************/
{
	int i,j;
	float v,e,c,w,bot;

	if (r[0] == 0.0) return;

	a[0] = 1.0;
	v = r[0];
	f[0] = g[0]/r[0];

	for (j=1; j<n; j++) {
		
		/* solve Ra=v as in Claerbout, FGDP, p. 57 */
		a[j] = 0.0;
		f[j] = 0.0;
		for (i=0,e=0.0; i<j; i++)
			e += a[i]*r[j-i];
		c = e/v;
		v -= c*e;
		for (i=0; i<=j/2; i++) {
			bot = a[j-i]-c*a[i];
			a[i] -= c*a[j-i];
			a[j-i] = bot;
		}

		/* use a and v above to get f[i], i = 0,1,2,...,j */
		for (i=0,w=0.0; i<j; i++)
			w += f[i]*r[j-i];
		c = (w-g[j])/v;
		for (i=0; i<=j; i++)
			f[i] -= c*a[j-i];
	}
}


