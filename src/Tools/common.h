#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <errno.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <limits.h>
#include <float.h>
// #include "omp.h"
#include "par.h"
#ifndef ALLOC_H
#define ALLOC_H
#include "alloc.h"
#endif
#ifndef COMPLEX_H
#define COMPLEX_H
#include "complex.h"
#endif
#ifndef FRANUNI_H
#define FRANUNI_H
#include "franuni.h"  
#endif
#ifndef SINC_H
#define SINC_H
#include "sinc.h"
#endif
#ifndef FFT_H
#define FFT_H
#include "fft.h"
#endif



#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))
#define	MAX(x,y) ((x) > (y) ? (x) : (y))
#define	MIN(x,y) ((x) < (y) ? (x) : (y))
#define PI (3.141592653589793)
#ifndef ABS
#define ABS(x) ((x) < 0 ? -(x) : (x))
#endif

