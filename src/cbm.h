#include "main.h"

/* Ray types */
/* one step along ray */
typedef struct RayStepStruct
{
	float t;			  /* time */
	float a;			  /* angle */
	float x, z;			  /* x,z coordinates */
	float q1, p1, q2, p2; /* Cerveny's dynamic ray tracing solution */
	int kmah;			  /* KMAH index */
	float c, s;			  /* cos(angle) and sin(angle) */
	float v, dvdx, dvdz;  /* velocity and its derivatives */
} RayStep;

/* one ray */
typedef struct RayStruct
{
	int nrs;	 /* number of ray steps */
	RayStep *rs; /* array[nrs] of ray steps */
	int nc;		 /* number of circles */
	int ic;		 /* index of circle containing nearest step */
	void *c;	 /* array[nc] of circles */
} Ray;

/* size of cells in which to linearly interpolate complex time and amplitude */
#define CELLSIZE 6

/* the dominant frequency of ricker wavelet */
#define DOMINFY 30

/* define the forward modeling and backward migration */
#define MOD 1
#define MIG -1

/* factor by which to oversample time for linear interpolation of traces */
#define NOVERSAMPLE 8

/* number of exponential decay filters */
#define NFILTER 1

/* exp(EXPMIN) is assumed to be negligible */
#define EXPMIN (-5.0)

/* parameters of shot gathers */
typedef struct ShotParStruct
{
	int index;
	int ntr;
	long ishift;
} ShotPar;

/* filtered complex beam data as a function of real and imaginary time */
typedef struct BeamDataStruct
{
	int ntr;	  /* number of real time samples */
	float dtr;	  /* real time sampling interval */
	float ftr;	  /* first real time sample */
	int nti;	  /* number of imaginary time samples */
	float dti;	  /* imaginary time sampling interval */
	float fti;	  /* first imaginary time sample */
	complex **cf; /* array[npx][nti][ntr] of complex data */
} BeamData;

/* one cell in which to linearly interpolate complex time and amplitude */
typedef struct CellStruct
{
	int live;	 /* random number used to denote a live cell */
	int ip;		 /* parameter used to denote ray parameter */
	float tr;	 /* real part of traveltime */
	float ti;	 /* imaginary part of traveltime */
	float ar;	 /* real part of amplitude */
	float ai;	 /* imaginary part of amplitude */
	float angle; /*angle of beam*/
} Cell;

/* structure containing information used to set and fill cells */
typedef struct CellsStruct
{
	int nt;		 /* number of time samples */
	float dt;	 /* time sampling interval */
	float ft;	 /* first time sample */
	int lx;		 /* number of x samples per cell */
	int mx;		 /* number of x cells */
	int nx;		 /* number of x samples */
	float dx;	 /* x sampling interval */
	float fx;	 /* first x sample */
	int lz;		 /* number of z samples per cell */
	int mz;		 /* number of z cells */
	int nz;		 /* number of z samples */
	float dz;	 /* z sampling interval */
	float fz;	 /* first z sample */
	int live;	 /* random number used to denote a live cell */
	int ip;		 /* parameter used to denote ray parameter */
	float wmin;	 /* minimum (reference) frequency */
	float lmin;	 /* minimum beamwidth for frequency wmin */
	Cell **cell; /* cell array[mx][mz] */
	Ray *ray;	 /* ray */
	//	BeamData *bd;	/* complex beam data as a function of complex time */
	//	float **g;	/* array[nx][nz] containing g(x,z) */
} Cells;

/*conjugate-direction descent*/
void cgstep(int iter, float **x, float **g, float **rr, float **gg, int nx, int nz, int ntr, int nt, float **s, float **ss);

/* Input the shot gather and head info*/
void tripd(float *d, float *e, float *b, int n);
void smooth2d(int n1, int n2, float r1, float r2, float **v);
void partall(int type, int nz, int cdpmin, int apernum, int apermin, float **part, float **all);
void inputrace(int is, int nt, float dx, int maxtr, FILE *fp, float *spx, float *rpxmin, float *rpxmax, int *nistr);

/* tapering & filtering */
void hanning_filter(int ftype, float f1, float f2, float f3, float f4, float fw, float dw, int nw, complex *spec);
void hanning_taper(int nz, int nx, int ntaper, float **g);

/* Ray functions */
Ray *makeRay(float x0, float z0, float a0, int nt, float dt, float ft,
			 int nx, float dx, float fx, int nz, float dz, float fz, float **v);
void freeRay(Ray *ray);
int nearestRayStep(Ray *ray, float x, float z);

/* Velocity functions */
void *vel2Alloc(int nx, float dx, float fx,
				int nz, float dz, float fz, float **v);
void vel2Free(void *vel2);
void vel2Interp(void *vel2, float x, float z,
				float *v, float *vx, float *vz, float *vxx, float *vxz, float *vzz);
void convol(float **f, float **imb, int npx, int ntau, float ftau, float dtau);

/* Beam functions */
void beam2trc(int type, float bwh, float fxb, float dxb, int nxb, float fmin,
			  int nt, float dt, float ft, int nx, float dx, float fx, float **f,
			  int ntau, float dtau, float ftau,
			  int npx, float dpx, float fpx, float ***beam, float *head);
void cbm2trc(int type, float bwh, float fxb, float dxb, int nxb, float fmin,
			 int nt, float dt, float ft, int nx, float dx, float fx, float **f,
			 int ntau, float dtau, float ftau,
			 int npx, float dpx, float fpx, float ***beam, float *head);
void accray(Ray *ray, Cell **c, float fmin, float lmin, int lx, int lz, int mx,
			int mz, int live, int ipx, int nt, float dt, float ft,
			int nx, float dx, float fx, int nz, float dz, float fz, float **g, float **v);
void scanimg(int type, float spx, float bcx, float dx, float dz, Cell ***c1, Cell ***c2, int live, int dead, int nt, float dt, float ft,
			 float fmin, int lx, int lz, int nx, int nz, int mx, int mz, int npx, float **f, float **g,
			 float fpx, float dpx, float **imb, float **ampgrid,float **tang);

/* functions defined and used internally */
void csmiggb(int type, float bwh, float fmin, float fmax, float amin, float amax, int live,
					int dead, int nt, float dt, float spx, float rpxmin, float rpxmax, int nx, float dx, int ntr, float dtr, int nz, float dz,
					float **f, float **v, float **g, float *head,float **tang);

/* circle for efficiently finding nearest ray step */
typedef struct CircleStruct
{
	int irsf; /* index of first raystep in circle */
	int irsl; /* index of last raystep in circle */
	float x;  /* x coordinate of center of circle */
	float z;  /* z coordinate of center of circle */
	float r;  /* radius of circle */
} Circle;

/* functions defined and used internally */
Circle *makeCircles(int nc, int nrs, RayStep *rs);

/* functions defined and used internally */
static void xtop(float w,
				 int nx, float dx, float fx, complex *g,
				 int np, float dp, float fp, complex *h);
static BeamData *beamData(int type, int npx, float wmin, int nt, float dt, float ft, float **f);
static void setCell(Cells *cells, int jx, int jz);
static void accCell(Cells *cells, int jx, int jz);
static int cellTimeAmp(Cells *cells, int jx, int jz);
static void cellBeam(int type, float spx, float bcx, float dx, float dz, Cell **cell, float **f, float **g, int jx, int jz, BeamData *bd,
					 int lx, int lz, int nx, int nz, int npx, int ips, int ipr, float fpx, float dpx, float ft, float dt, int nt, float **imb,float **tang);

void csmiggb_irr(int type, float bwh, float fmin, float fmax, float amin, float amax, int live,
						int dead, int nt, float dt, float spx, float rpxmin, float rpxmax, int nx, float dx, int ntr, float dtr, int nz, float dz,
						float **f, float **v, float **g, float *coordx, float *elevtemp,float *elev, float nxline, float dxline,float **tang);
void cbm2trc_irr(int type, float bwh, float fxb, float dxb, int nxb, float fmin,
				 int nt, float dt, float ft, int nx, float dx, float fx, float dxx, float **f,
				 int ntau, float dtau, float ftau,
				 int npx, float dpx, float fpx, float ***beam, float *head, float *elev, float dz, float **v);
void beam2trc_irr(int type, float bwh, float fxb, float dxb, int nxb, float fmin,
				  int nt, float dt, float ft, int nx, float dx, float fx, float dxx, float **f,
				  int ntau, float dtau, float ftau,
				  int npx, float dpx, float fpx, float ***beam, float *head, float *elev, float dz, float **v);

/* hermite polynomials */
static float h00(float x) { return 2.0 * x * x * x - 3.0 * x * x + 1.0; }
static float h01(float x) { return 6.0 * x * x - 6.0 * x; }
static float h02(float x) { return 12.0 * x - 6.0; }
static float h10(float x) { return -2.0 * x * x * x + 3.0 * x * x; }
static float h11(float x) { return -6.0 * x * x + 6.0 * x; }
static float h12(float x) { return -12.0 * x + 6.0; }
static float k00(float x) { return x * x * x - 2.0 * x * x + x; }
static float k01(float x) { return 3.0 * x * x - 4.0 * x + 1.0; }
static float k02(float x) { return 6.0 * x - 4.0; }
static float k10(float x) { return x * x * x - x * x; }
static float k11(float x) { return 3.0 * x * x - 2.0 * x; }
static float k12(float x) { return 6.0 * x - 2.0; }