#ifndef COMMON_H
#define COMMON_H
#include "common.h"
#endif
#ifndef COMPLEX_H
#define COMPLEX_H
#include "complex.h"
#endif

/* Prime Factor FFTs */
int npfa (int nmin);
int npfao (int nmin, int nmax);
int npfar (int nmin);
int npfaro (int nmin, int nmax);
void pfacc (int isign, int n, complex z[]);
void pfarc (int isign, int n, float rz[], complex cz[]);
void pfacr (int isign, int n, complex cz[], float rz[]);
void pfa2cc (int isign, int idim, int n1, int n2, complex z[]);
void pfa2rc (int isign, int idim, int n1, int n2, float rz[], complex cz[]);
void pfa2cr (int isign, int idim, int n1, int n2, complex cz[], float rz[]);
void pfamcc (int isign, int n, int nt, int k, int kt, complex z[]);
