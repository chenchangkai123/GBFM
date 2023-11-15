#ifndef COMMON_H
#define COMMON_H
#include "common.h"
#endif

#ifndef COMPLEX_H
#define COMPLEX_H
#include "complex.h"
#endif

#ifndef ALLOC
#define ALLOC
void *aalloc1 (size_t n1, size_t size);
void *alloc1(int n1,int size);
void *realloc1(void *v,int n1,int size);
void **alloc2(int n1,int n2,int size);
void ***alloc3(int n1,int n2,int n3,int size);
void ****alloc4(int n1,int n2,int n3,int n4,int size);
void *****alloc5(int n1,int n2,int n3,int n4,int n5,int size);
void free1 (void *p);
void free2 (void **p);
void free3 (void ***p);
void free4 (void ****p);
void free5 (void *****p);

int *alloc1int(int n1);
void free1int(int *p);
int *realloc1int(int *v, int n1);

int **alloc2int(int n1, int n2);
void free2int(int **p);

int ***alloc3int(int n1, int n2, int n3);
void free3int(int ***p);

float *alloc1float(int n1);
float *realloc1float(float *v, int n1);
void free1float(float *p);

float **alloc2float(int n1, int n2);
void free2float(float **p);

float ***alloc3float(int n1, int n2, int n3);
void free3float(float ***p);

double *alloc1double(int n1);
double *realloc1double(double *v, int n1);
void free1double(double *p);

double **alloc2double(int n1, int n2);
void free2double(double **p);

double ***alloc3double(int n1, int n2, int n3);
void free3double(double ***p);

complex *alloc1complex(int n1);

void free1complex(complex *p);

complex **alloc2complex(int n1, int n2);

void free2complex(complex **p);

complex ***alloc3complex(int n1, int n2, int n3);

void free3complex(complex ***p);

complex ****alloc4complex(int n1, int n2, int n3, int n4);

void free4complex(complex ****p);

#endif

