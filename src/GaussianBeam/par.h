/* Copyright (c) Colorado School of Mines, 2010.*/
/* All rights reserved.                       */

/* par.h - include file for getpar, selfdoc, and error handling functions */

#ifndef PAR_H
#define PAR_H


/* INCLUDES */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <limits.h>
#include <float.h>
// #include "alloc.h"

#include <fcntl.h>      /* non-ANSI */
#include <unistd.h>     /* non-ANSI */
#include <sys/types.h>  /* non-ANSI */

/* GLOBAL DECLARATIONS */
extern int xargc; extern char **xargv;

/* ssize_t does not exist on some systems */
#ifndef ssize_t
#define ssize_t int
#endif



/* DEFINES */

/* getpar macros */
#define MUSTGETPARINT(x,y) if(!getparint(x,y)) err("must specify %s=",x)
#define MUSTGETPARFLOAT(x,y) if(!getparfloat(x,y)) err("must specify %s=",x)
#define MUSTGETPARSTRING(x,y) if(!getparstring(x,y)) err("must specify %s=",x)

#define STDIN (0)
#define STDOUT (1)
#define STDERR (2)


/* FUNCTION PROTOTYPES */


/* getpar parameter parsing */
void initargs (int argc, char **argv);
int getparint (char *name, int *p);
int getparuint (char *name, unsigned int *p);
int getparshort (char *name, short *p);
int getparushort (char *name, unsigned short *p);
int getparlong (char *name, long *p);
int getparulong (char *name, unsigned long *p);
int getparfloat (char *name, float *p);
int getpardouble (char *name, double *p);
int getparstring (char *name, char **p);
int getparstringarray (char *name, char **p);
int getnparint (int n, char *name, int *p);
int getnparuint (int n, char *name, unsigned int *p);
int getnparshort (int n, char *name, short *p);
int getnparushort (int n, char *name, unsigned short *p);
int getnparlong (int n, char *name, long *p);
int getnparulong (int n, char *name, unsigned long *p);
int getnparfloat (int n, char *name, float *p);
int getnpardouble (int n, char *name, double *p);
int getnparstring (int n, char *name, char **p);
int getnparstringarray (int n, char *name, char **p);
int getnpar (int n, char *name, char *type, void *ptr);
int countparname (char *name);
int countparval (char *name);
int countnparval (int n, char *name);

/* For ProMAX */
void getPar(char *name, char *type, void *ptr);

/* errors and warnings */
void err (char *fmt, ...);
void syserr (char *fmt, ...);
void warn (char *fmt, ...);

/* string to numeric conversion with error checking */
short eatoh(char *s);
unsigned short eatou(char *s);
int eatoi(char *s);
unsigned int eatop(char *s);
long eatol(char *s);
unsigned long eatov(char *s);
float eatof(char *s);
double eatod(char *s);

/* system calls with error trapping */
int ecreat(char *path, int perms);
int efork(void);
int eopen(char *path, int flags, int perms);
int eclose(int fd);
int eunlink(char *path);
off_t elseek(int fd, off_t offset, int origin);
int epipe(int fd[2]);

ssize_t eread(int fd, char *buf, size_t nbytes);
ssize_t ewrite(int fd, char *buf, size_t nbytes);

/* system subroutine calls with error trapping */
FILE *efopen(const char *file, const char *mode);
FILE *efreopen(const char *file, const char *mode, FILE *stream1);
FILE *efdopen(int fd, const char *mode);
FILE *epopen(char *command, char *type);
int efclose(FILE *stream);
int epclose(FILE *stream);
int efflush(FILE *stream);
int eremove(const char *file);
int erename(const char *oldfile, const char* newfile);
int efseeko(FILE *stream, off_t offset, int origin);
int efseek(FILE *stream, off_t offset, int origin);
long eftell(FILE *stream);
off_t eftello(FILE *stream);
void erewind(FILE *stream);
FILE *etmpfile(void);
char *emkstemp(char *namebuffer);
void *emalloc(size_t size);
void *erealloc(void *memptr, size_t size);
void *ecalloc(size_t count, size_t size);
size_t efread(void *bufptr, size_t size, size_t count, FILE *stream);
size_t efwrite(void *bufptr, size_t size, size_t count, FILE *stream);

/* some SUN users may need to comment out the next two items, */
/* if your system does not have "fgetpos()" and "fsetpos()" defined. */
/* You will also need to comment out the lines defining the functions */
/* efgetpos() and efsetpos() in CWPROOT/src/par/lib/subcalls.c */

/* Modified: 21 June 1995: */
/* so that setting -DSUN_A in Makefile.config should make this unnecessary */
/* CWP: John Stockwell */

#ifndef SUN_A
int efgetpos(FILE *stream, fpos_t *position);
int efsetpos(FILE *stream, const fpos_t *position);
#endif

/* allocation with error trapping */

#endif /* PAR_H */
