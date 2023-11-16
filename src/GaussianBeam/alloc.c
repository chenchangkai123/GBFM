/*****************************************************

void *alloc1(int n1,int size)
void *realloc1(void *v,int n1,int size)
void **alloc2(int n1,int n2,int size)
void ***alloc3(int n1,int n2,int n3,int size)
void free1 (void *p)
void free2 (void **p)
void free3 (void ***p)
int *alloc1int(int n1)
void free1int(int *p)
int *realloc1int(int *v, int n1)
int **alloc2int(int n1, int n2)
void free2int(int **p)
int ***alloc3int(int n1, int n2, int n3)
void free3int(int ***p)
float *alloc1float(int n1)
float *realloc1float(float *v, int n1)
void free1float(float *p)
float **alloc2float(int n1, int n2)
void free2float(float **p)
float ***alloc3float(int n1, int n2, int n3)
void free3float(float ***p)
double *alloc1double(int n1)
double *realloc1double(double *v, int n1)
void free1double(double *p)
double **alloc2double(int n1, int n2)
void free2double(double **p)
double ***alloc3double(int n1, int n2, int n3)
void free3double(double ***p)
complex *alloc1complex(int n1)
void free1complex(complex *p)
complex **alloc2complex(int n1, int n2)
void free2complex(complex **p)
complex ***alloc3complex(int n1, int n2, int n3)
void free3complex(complex ***p)

******************************************************/

/* ***************************************************
Author :     LEON LI and ZHANGJT
Data   :     2005.3.24
		    China University of Petroleum(East China)

All the functions can be found in su software by
cwp, so it can be used freely if you want to use 
these function to program free softare or commercial
softare 
*******************************************************/
#include "alloc.h"
/* allocate a 1-d array */

void *aalloc1 (size_t n1, size_t size)
{
	void *p;

	if ((p=malloc(n1*size))==NULL)
		return NULL;
	return p;
}

void *alloc1 (int n1, int size)
{
	void *p;

	if ((p=malloc(n1*size))==NULL)
		return NULL;
	return p;
}

/* re-allocate a 1-d array */
void *realloc1(void *v, int n1, int size)
{
	void *p;

	if ((p=realloc(v,n1*size))==NULL)
		return NULL;
	return p;
}

/* free a 1-d array */
void free1 (void *p)
{
	free(p);
}

/* allocate a 2-d array */
void **alloc2 (int n1, int n2, int size)
{
	int i2;
	void **p;

	if ((p=(void**)malloc(n2*sizeof(void*)))==NULL) 
		return NULL;
	if ((p[0]=(void*)malloc(n2*n1*size))==NULL) {
		free(p);
		return NULL;
	}
	for (i2=0; i2<n2; i2++)
		p[i2] = (char*)p[0]+size*n1*i2;
	return p;
}

/* free a 2-d array */
void free2 (void **p)
{
	free(p[0]);
	free(p);
}

/* allocate a 3-d array */
void ***alloc3 (int n1, int n2, int n3, int size)
{
	int i3,i2;
	void ***p;

	if ((p=(void***)malloc(n3*sizeof(void**)))==NULL)
		return NULL;
	if ((p[0]=(void**)malloc(n3*n2*sizeof(void*)))==NULL) {
		free(p);
		return NULL;
	}
	if ((p[0][0]=(void*)malloc(n3*n2*n1*size))==NULL) {
		free(p[0]);
		free(p);
		return NULL;
	}

	for (i3=0; i3<n3; i3++) {
		p[i3] = p[0]+n2*i3;
		for (i2=0; i2<n2; i2++)
			p[i3][i2] = (char*)p[0][0]+size*n1*(i2+n2*i3);
	}
	return p;
}

/* free a 3-d array */
void free3 (void ***p)
{
	free(p[0][0]);
	free(p[0]);
	free(p);
}
/* allocate a 4-d array */
void ****alloc4 (int n1, int n2, int n3, int n4, int size)
{
	int i4,i3,i2;
	void ****p;

	if ((p=(void****)malloc(n4*sizeof(void***)))==NULL)
		return NULL;
	if ((p[0]=(void***)malloc(n4*n3*sizeof(void**)))==NULL) {
		free(p);
		return NULL;
	}
	if ((p[0][0]=(void**)malloc(n4*n3*n2*sizeof(void*)))==NULL) {
		free(p[0]);
		free(p);
		return NULL;
	}
	if ((p[0][0][0]=(void*)malloc(n4*n3*n2*n1*size))==NULL) {
		free(p[0][0]);
		free(p[0]);
		free(p);
		return NULL;
	}
	for (i4=0; i4<n4; i4++) {
		p[i4] = p[0]+i4*n3;
		for (i3=0; i3<n3; i3++) {
			p[i4][i3] = p[0][0]+n2*(i3+n3*i4);
			for (i2=0; i2<n2; i2++)
				p[i4][i3][i2] = (char*)p[0][0][0]+
						size*n1*(i2+n2*(i3+n3*i4));
		}
	}
	return p;
}

/* free a 4-d array */
void free4 (void ****p)
{
	free(p[0][0][0]);
	free(p[0][0]);
	free(p[0]);
	free(p);
}
void *****alloc5 (int n1, int n2, int n3, int n4, int n5, int size)
{
	int i5,i4,i3,i2;
	long long ntotal,total;
	long long nsamp,samp;
	void *****p;

	if ((p=(void*****)malloc(n5*sizeof(void****)))==NULL)
		return NULL;
	if ((p[0]=(void****)malloc(n5*n4*sizeof(void***)))==NULL) {
		free(p);
		return NULL;
	}
	if ((p[0][0]=(void***)malloc(n5*n4*n3*sizeof(void**)))==NULL) {
		free(p[0]);
		free(p);
		return NULL;
	}
	if ((p[0][0][0]=(void**)malloc(n5*n4*n3*n2*sizeof(void*)))==NULL) {
		free(p[0][0]);
		free(p[0]);
		free(p);
		return NULL;
	}
	total = n5*n4*n3*n2*n1;
	ntotal = total*size;
	//if ((p[0][0][0][0]=(void*)malloc(n5*n4*n3*n2*n1*size))==NULL) {
	if ((p[0][0][0][0]=(void*)malloc(ntotal))==NULL) {
		free(p[0][0][0]);
		free(p[0][0]);
		free(p[0]);
		free(p);
		return NULL;
	}
	for (i5=0; i5<n5; i5++) {
		p[i5] = p[0]+n4*i5;
		for (i4=0; i4<n4; i4++) {
			p[i5][i4] = p[0][0]+n3*(i4+n4*i5);
			for (i3=0; i3<n3; i3++) {
				p[i5][i4][i3] = p[0][0][0]+n2*(i3+n3*(i4+n4*i5));
				for (i2=0; i2<n2; i2++) {
					samp = n1*(i2+n2*(i3+n3*(i4+n4*i5)));
					nsamp = size*samp;
					p[i5][i4][i3][i2] = (char*)p[0][0][0][0]+
						nsamp;
//size*n1*(i2+n2*(i3+n3*(i4+n4*i5)));
				}
			}
		}
	}
	return p;
}

/* free a 5-d array */
void free5 (void *****p)
{
	free(p[0][0][0][0]);
	free(p[0][0][0]);
	free(p[0][0]);
	free(p[0]);
	free(p);
}




/* allocate a 1-d array of ints */
int *alloc1int(int n1)
{
	return (int*)alloc1(n1,sizeof(int));
}

/* re-allocate a 1-d array of ints */
int *realloc1int(int *v, int n1)
{
	return (int*)realloc1(v,n1,sizeof(int));
}

/* free a 1-d array of ints */
void free1int(int *p)
{
	free1(p);
}

/* allocate a 2-d array of ints */
int **alloc2int(int n1, int n2)
{
	return (int**)alloc2(n1,n2,sizeof(int));
}

/* free a 2-d array of ints */
void free2int(int **p)
{
	free2((void**)p);
}

/* allocate a 3-d array of ints */
int ***alloc3int(int n1, int n2, int n3)
{
	return (int***)alloc3(n1,n2,n3,sizeof(int));
}

/* free a 3-d array of ints */
void free3int(int ***p)
{
	free3((void***)p);
}

/* allocate a 1-d array of floats */
float *alloc1float(int n1)
{
	return (float*)alloc1(n1,sizeof(float));
}

/* re-allocate a 1-d array of floats */
float *realloc1float(float *v, int n1)
{
	return (float*)realloc1(v,n1,sizeof(float));
}

/* free a 1-d array of floats */
void free1float(float *p)
{
	free1(p);
}

/* allocate a 2-d array of floats */
float **alloc2float(int n1, int n2)
{
	return (float**)alloc2(n1,n2,sizeof(float));
}

/* free a 2-d array of floats */
void free2float(float **p)
{
	free2((void**)p);
}

/* allocate a 3-d array of floats */
float ***alloc3float(int n1, int n2, int n3)
{
	return (float***)alloc3(n1,n2,n3,sizeof(float));
}

/* free a 3-d array of floats */
void free3float(float ***p)
{
	free3((void***)p);
}

/* allocate a 1-d array of doubles */
double *alloc1double(int n1)
{
	return (double*)alloc1(n1,sizeof(double));
}

/* re-allocate a 1-d array of doubles */
double *realloc1double(double *v, int n1)
{
	return (double*)realloc1(v,n1,sizeof(double));
}


/* free a 1-d array of doubles */
void free1double(double *p)
{
	free1(p);
}

/* allocate a 2-d array of doubles */
double **alloc2double(int n1, int n2)
{
	return (double**)alloc2(n1,n2,sizeof(double));
}

/* free a 2-d array of doubles */
void free2double(double **p)
{
	free2((void**)p);
}

/* allocate a 3-d array of doubles */
double ***alloc3double(int n1, int n2, int n3)
{
	return (double***)alloc3(n1,n2,n3,sizeof(double));
}

/* free a 3-d array of doubles */
void free3double(double ***p)
{
	free3((void***)p);
}

/* allocate a 1-d array of complexs */
complex *alloc1complex(int n1)
{
	return (complex*)alloc1(n1,sizeof(complex));
}

/* re-allocate a 1-d array of complexs */
complex *realloc1complex(complex *v, int n1)
{
	return (complex*)realloc1(v,n1,sizeof(complex));
}

/* free a 1-d array of complexs */
void free1complex(complex *p)
{
	free1(p);
}

/* allocate a 2-d array of complexs */
complex **alloc2complex(int n1, int n2)
{
	return (complex**)alloc2(n1,n2,sizeof(complex));
}

/* free a 2-d array of complexs */
void free2complex(complex **p)
{
	free2((void**)p);
}

/* allocate a 3-d array of complexs */
complex ***alloc3complex(int n1, int n2, int n3)
{
	return (complex***)alloc3(n1,n2,n3,sizeof(complex));
}

/* free a 3-d array of complexs */
void free3complex(complex ***p)
{
	free3((void***)p);
}
/* allocate a 4-d array of complexs */
complex ****alloc4complex(int n1, int n2, int n3, int n4)
{
        return (complex****)alloc4(n1,n2,n3,n4,sizeof(complex));
}

/* free a 4-d array of complexs */
void free4complex(complex ****p)
{
        free4((void****)p);
}
