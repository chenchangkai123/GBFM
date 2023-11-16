#include "main.h"
#include <time.h>
#include "mpi.h"
#include "par.h"

int xargc;
char **xargv;

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

typedef struct headInformation
{
	int is;
	float isp;
	float firec;
	int ntr;
	int tataltr;
} headI;

/* size of cells in which to linearly interpolate complex time and amplitude */
#define CELLSIZE 6

/* the dominant frequency of ricker wavelet */
#define DOMINFY 17

/* define the forward modeling and backward migration */
#define MOD 1
#define MIG -1

/* factor by which to oversample time for linear interpolation of traces */
#define NOVERSAMPLE 4

/* number of exponential decay filters */
#define NFILTER 10

/* exp(EXPMIN) is assumed to be negligible */
#define EXPMIN (-5.0)

/* filtered complex beam data as a function of real and imaginary time */
typedef struct BeamDataStruct
{
	int ntr;	   /* number of real time samples */
	float dtr;	   /* real time sampling interval */
	float ftr;	   /* first real time sample */
	int nti;	   /* number of imaginary time samples */
	float dti;	   /* imaginary time sampling interval */
	float fti;	   /* first imaginary time sample */
	complex ***cf; /* array[npx][nti][ntr] of complex data */
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

/* Input the shot gather and head info*/
void tripd(float *d, float *e, float *b, int n);

void smooth2d(int n1, int n2, float r1, float r2, float **v);

void partall(int type, int nz, int apernum, int apermin, float **part, float **all);

void inputrace(int is, int nt, float dx, int maxtr, FILE *fp, int *sisp, int *firisp, int *nistr);

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
void synData(float bwh, float dxb, float fmin,
			 int nt, float dt, float ft,
			 int nx, float dx, float fx, float **f,
			 int ntau, float dtau, float ftau,
			 int npx, float dpx, float fpx, float **g);

/* Beam functions */
void formBeams(float bwh, float dxb, float fmin,
			   int nt, float dt, float ft,
			   int nx, float dx, float fx, float **f,
			   int ntau, float dtau, float ftau,
			   int npx, float dpx, float fpx, float **g);
void accray(Ray *ray, Cell **c, float fmin, float lmin, int lx, int lz, int mx,
			int mz, int live, int ipx, int nt, float dt, float ft,
			int nx, float dx, float fx, int nz, float dz, float fz, float **g, float **v);
void scanimg(int type, Cell ***c1, Cell ***c2, int live, int dead, int nt, float dt, float ft,
			 float fmin, int lx, int lz, int nx, int nz, int mx, int mz, int npx, float **f, float **g,
			 float fpx, float dpx, float **imb);

/* functions defined and used internally */
static void csmiggb(int type, float bwh, float fmin, float fmax, float amin, float amax, int live,
					int dead, int nt, float dt, int sp, int firp, int nx, float dx, int ntr, float dtr, int nz, float dz,
					float **f, float **v, float **g);

int main(int argc, char *argv[])
{
	int nx, nz, nt, ix, iz, nshot, ishot, ntr, itr, sp, iangle, i, j, tracepad;
	int firp, firstp, lastp; /* first and last receive cdp of one shot gather*/
	float dx, dz, dt, fmin, fmax, amin, amax, vavg, bwh, dtr;
	float **f, **orig, **lsmg, **refl;
	float **vsm;
	float **tempv, **tempg, **tempf, **templsmg;
	float **s, **ss;
	float **origall, **lsmgall;
	char *vfile, *seifile, *imgfile, *reftfile, *orgfile;
	char *orimigfile;
	char *initreflec;
	int live, dead, maxtrace;
	int naper, apermin, apermax, cdpmin, apernum;
	int nangle;
	int iter, nter;
	float ft;
	clock_t start, finish;
	double duration, res = 0.0;
	float **orif;
	double dsum = 0.0, esum = 0.0;
	float beta;
	float **oriv;
	float **tt;
	double *residual;
	int kt, totaltr;
	int np, myid;

	FILE *vfp;
	FILE *sfp, *dfp;
	// FILE *rfp, *ufp, *mfp; /* temp file to hold traces	*/
	FILE *headRecfp, *headShotfp, *headInfofp;
	char *headRec, *headShot, *headInfo;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	initargs(argc, argv);

	if (!getparint("maxtrace", &maxtrace))
		maxtrace = 20000;
	if (!getparint("cdpmin", &cdpmin))
		cdpmin = 0;
	if (!getparfloat("ft", &ft))
		ft = 1.0;
	if (!getparint("naper", &naper))
		naper = 150;
	if (!getparint("nt", &nt))
		nt = 750;
	if (!getparint("nshot", &nshot))
		nshot = 150;
	if (!getparint("ntr", &ntr))
		ntr = 100;
	if (!getparint("nz", &nz))
		nz = 750;
	if (!getparint("nx", &nx))
		nx = 737;
	if (!getparint("nangle", &nangle))
		nangle = 40;
	if (!getparfloat("dt", &dt))
		dt = 0.004;
	if (!getparfloat("dx", &dx))
		dx = 12.5 * ft;
	if (!getparfloat("dz", &dz))
		dz = 4.0 * ft;
	if (!getparfloat("dtr", &dtr))
		dtr = 12.5 * ft;
	if (!getparfloat("fmin", &fmin))
		fmin = 0.04 / dt;
	if (!getparfloat("fmax", &fmax))
		fmax = 5.0 * fmin;
	if (!getparfloat("amax", &amax))
		amax = 50.0;
	if (!getparfloat("amin", &amin))
		amin = -amax;

	if (!getparstring("vfile", &vfile))
		vfile = "mar_vel_737_750_sun.sgy";
	if (!getparstring("seifile", &seifile))
		seifile = "modeled.su";
	if (!getparstring("imgfile", &imgfile))
		imgfile = "mig.dat";
	// if (!getparstring("reftfile", &reftfile))
	// 	reftfile = "reflectivity.dat";
	// if (!getparstring("orgfile", &orgfile))
	// 	orgfile = "modeled.su";
	if (!getparstring("headRec", &headRec))
		headRec = "headRec.dat";
	if (!getparstring("headShot", &headShot))
		headShot = "headShot.dat";
	if (!getparstring("headInfo", &headInfo))
		headInfo = "headInfo";

	/* open files */
	vfp = fopen(vfile, "rb+");
	sfp = fopen(seifile, "rb+");
	dfp = fopen(imgfile, "wb+");
	// rfp = fopen(reftfile, "wb+");
	// ufp = fopen(orgfile, "wb+");
	headInfofp = fopen(headInfo, "r");
	// headRecfp = fopen(headRec, "rb+");
	// headShotfp = fopen(headShot, "rb+");

	headI *headInformation;
	headInformation = (headI *)malloc(sizeof(headI) * nshot);
	for (i = 0; i < nshot; i++)
	{
		fscanf(headInfofp, "%d,%f,%f,%d,%d\n", &headInformation[i].is, &headInformation[i].isp, &headInformation[i].firec, &headInformation[i].ntr, &headInformation[i].tataltr);
		// printf("%d %f %f %d\n",headInformation[i].is,headInformation[i].isp,headInformation[i].firec,headInformation[i].ntr);
	}
	fclose(headInfofp);

	/* random numbers used to denote live and dead cells */
	live = 1 + (int)(1.0e7 * franuni());
	dead = 1 + (int)(1.0e7 * franuni());

	/* allocate workspace */
	// refl = alloc2float(nz, nx);
	orig = alloc2float(nz, nx);
	vsm = alloc2float(nz, nx);
	oriv = alloc2float(nz, nx);
	origall = alloc2float(nz, nx);

	/* load original velocity and close tmpfile */
	for (ix = 0; ix < nx; ix++)
	{
		fread(oriv[ix], sizeof(float), nz, vfp);
	}
	fclose(vfp);

	for (ix = 0, vavg = 0.0; ix < nx; ++ix)
	{
		for (iz = 0; iz < nz; ++iz)
		{
			vavg += oriv[ix][iz];
		}
	}

	vavg /= nx * nz;

	/* get beam half-width */
	bwh = vavg / fmin / 1.2;
	printf("Beam width is %f\n", bwh);

	/* smooth and get the background velocity */
	for (ix = 0; ix < nx; ++ix)
	{
		for (iz = 0; iz < nz; ++iz)
		{
			vsm[ix][iz] = oriv[ix][iz];
		}
	}
	smooth2d(nz, nx, 8.0, 8.0, vsm);

#if 0
	/* calculate the reflectivity and write out */
	for (ix = 0; ix < nx; ix++)
	{
		for (iz = 1; iz < nz; iz++)
		{
			refl[ix][iz] = (oriv[ix][iz] - oriv[ix][iz - 1]) / (oriv[ix][iz] + oriv[ix][iz - 1]);
		}
		for (iz = 0; iz < 1; iz++)
		{
			refl[ix][iz] = refl[ix][iz + 1];
		}
	}

	for (ix = 0; ix < nx; ix++)
	{
		fwrite(refl[ix], sizeof(float), nz, rfp);
	}
	fclose(rfp);
#endif

	// for (ix=0; ix<nx; ix++) {
	// 		for (iz=1; iz<nz; iz++) {
	// 			fread(&refl[ix][iz],sizeof(float),1,initre);
	// 		}
	// 		for (iz=0; iz<1; iz++) {
	// 			refl[ix][iz]=refl[ix][iz+1];
	// 		}
	// 	}

	/* start counting time */
	if (myid == 0)
		start = clock();

	/* alloc memory for iteration */
	memset(&orig[0][0], 0, sizeof(float) * nx * nz);
	memset(&origall[0][0], 0, sizeof(float) * nx * nz);

	/* First step, calculate modeled data and migrate */
	for (ishot = 0 + myid; ishot < nshot; ishot = ishot + np)
	{
		// inputrace(ishot,nt,dx,maxtrace,sfp,&sp,&firp,&ntr);
		sp = headInformation[ishot].isp / dx;
		firp = headInformation[ishot].firec / dx;
		ntr = headInformation[ishot].ntr;
		totaltr = headInformation[ishot].tataltr;

		// printf("sp=%d firp=%d ntr=%d\n", sp, firp, ntr);

		sp = sp - cdpmin;
		firp = firp - cdpmin;
		lastp = firp + NINT((ntr - 1) * dtr / dx);
		printf("firstp=%d lastp=%d\n", firp, lastp);

		if (sp > nx - 1)
		{
			printf("\033[1;35mshot point is out of range\n");
			printf("skip this shot\n\033[0m");
			continue;
		}
		firstp = MIN(sp, firp);

		apermin = MAX(firstp - naper, 0);
		apermax = MIN(lastp + naper, nx - 1);
		apernum = apermax - apermin + 1;

		firp -= apermin;
		sp -= apermin;

		orif = alloc2float(nt, ntr);
		tempv = alloc2float(nz, apernum);
		tempg = alloc2float(nz, apernum);
		memset(&tempg[0][0], 0, sizeof(float) * nz * apernum);

		printf("%d %d %d %d\n", apermin, apernum, firp, sp);

		/* just need part of the velocity */
		partall(1, nz, apernum, apermin, tempv, vsm);

		printf("First Step: shot number migrated is %d; shot point is %d\n", ishot, sp);
		memset(&orif[0][0], 0, sizeof(float) * ntr * nt);

#if 1
		/* input seismic data, then migrate them */
		fseek(sfp, sizeof(float) * totaltr * nt, SEEK_SET);
		for (i = 0; i < ntr; i++)
		{
			fread(orif[i], sizeof(float), nt, sfp);
		}
		fseek(sfp, 0, SEEK_SET);

		/* migrate the data*/
		memset(&tempg[0][0], 0, sizeof(float) * nz * apernum);
		csmiggb(MIG, bwh, fmin, fmax, amin, amax, live, dead, nt, dt, sp, firp,
				apernum, dx, ntr, dtr, nz, dz, orif, tempv, tempg);

		/* accmulate the subimage to image */
		partall(-1, nz, apernum, apermin, tempg, orig);
#endif

#if 0
		/* forward modeling, just need part of the reflectivity model */
		partall(1, nz, apernum, apermin, tempg, refl);

		/* forward model the reflectivity */
		// memset(&orif[0][0],0,sizeof(float)*ntr*nt);
		csmiggb(MOD, bwh, fmin, fmax, amin, amax, live, dead, nt, dt, sp, firp,
				apernum, dx, ntr, dtr, nz, dz, orif, tempv, tempg);

		for (i = 0; i < 0; i++)
		{
			kt = NINT(450 - i * 4.2);
			for (j = 0; j < kt; j++)
				orif[i][j] = 0.0;
		}

		for (ix = 0; ix < ntr; ix++)
		{
			fwrite(orif[ix], sizeof(float), nt, ufp);
		}
#endif

		/* free temp arrays */
		free2float(orif);
		free2float(tempv);
		free2float(tempg);

	} // end ishot

	MPI_Reduce(&orig[0][0], &origall[0][0], nx * nz, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	// write out LSM-image//
	if (myid == 0)
	{
		for (ix = 0; ix < nx; ix++)
		{
			fwrite(origall[ix], sizeof(float), nz, dfp);
		}
	}

	fclose(dfp);
	fclose(sfp);
	// fclose(ufp);

	if (myid == 0)
	{
		finish = clock();
		duration = (double)(finish - start) / CLOCKS_PER_SEC;
		printf("migration time is %f seconds\n", duration);
	}

	/* free workspace */
	free2float(oriv);
	// free2float(refl);
	free2float(orig);
	free2float(vsm);
	free2float(origall);
	free(headInformation);
	MPI_Finalize();

	return EXIT_SUCCESS;
}

static void csmiggb(int type, float bwh, float fmin, float fmax, float amin, float amax, int live,
					int dead, int nt, float dt, int sp, int firp, int nx, float dx, int ntr, float dtr, int nz, float dz,
					float **f, float **v, float **g)
/*****************************************************************************
Migrate data via accumulation of Gaussian beams.
******************************************************************************
	  input:
type            type of computation
bwh		horizontal beam half-width at surface z=0
fmin		minimum frequency (cycles per unit time)
fmax		maximum frequency (cycles per unit time)
amin		minimum emergence angle at surface z=0 (degrees)
amax		maximum emergence angle at surface z=0 (degrees)
nt		number of time samples
dt		time sampling interval (first time assumed to be zero)
nx		number of x samples
dx		x sampling interval
ntr     trace samples
dtr     trace sampling interval
nz		number of z samples
dz		z sampling interval
f		array[nx][nt] containing zero-offset data f(t,x)
v		array[nx][nz] containing half-velocities v(x,z)

Output:
g		array[nx][nz] containing migrated image g(x,z)
*****************************************************************************/
{
	int nxb, npx, ntau, ipx, ix, ixb, ixlo, ixhi, nxw, iz, mx, mz, lx, lz, rayxp;
	int i, ixx, izz, irs, itau;
	int imx, imz;
	float ft, fx, fz, xwh, dxb, fxb, xb, vmin, dpx, fpx, px,
		taupad, dtau, ftau, fxw, pxmin, pxmax,
		a0, x0, z0, bwhc, **b, **imb, x;
	Ray *ray1, *ray2;
	Cell ***cell1, ***cell2;
	float **ampgrid, **ampr, **ampi;
	float ar, ai, tr, ti, pz, ampmax;
	float wmin = 2 * PI * fmin;

	// FILE *cellfp;
	// cellfp = fopen("cell.bin", "wb+");

	// FILE *rayxfp;
	// rayxfp = fopen("x.bin", "wb+");

	// FILE *rayzfp;
	// rayzfp = fopen("z.bin", "wb+");

	/* first t, x, and z assumed to be zero */
	ft = fx = fz = 0.0;

	/* size of coarse gird */
	lx = CELLSIZE;
	lz = CELLSIZE;

	/* number of cells */
	mx = 2 + (nx - 1) / lx;
	mz = 2 + (nz - 1) / lz;

	/* convert minimum and maximum angles to radians */
	amin *= PI / 180.0;
	amax *= PI / 180.0;
	if (amin > amax)
	{
		float atemp = amin;
		amin = amax;
		amax = atemp;
	}

	/* window half-width */
	xwh = 3.0 * bwh;

	/* beam center sampling */
	dxb = NINT((2.0 * bwh * sqrt(fmin / fmax)) / dtr) * dtr;
	nxb = 1 + (ntr - 1) * dtr / dxb;
	fxb = fx + 0.5 * ((ntr - 1) * dtr - (nxb - 1) * dxb);

	/* minimum velocity at surface z=0 */
	for (ix = 1, vmin = v[0][0]; ix < nx; ++ix)
		if (v[ix][0] < vmin)
			vmin = v[ix][0];

	printf("vmin=%f %d %d\n", vmin, mx, mz);

	/* beam sampling */
	pxmin = sin(amin) / vmin;
	pxmax = sin(amax) / vmin;
	//	dpx = 2.0/(2.0*xwh*sqrt(fmin*fmax));
	dpx = 1.0 / (6.0 * bwh * sqrt(fmin * fmax));
	npx = 1 + (pxmax - pxmin) / dpx;
	fpx = pxmin + 0.5 * (pxmax - pxmin - (npx - 1) * dpx);
	taupad = MAX(ABS(pxmin), ABS(pxmax)) * xwh;
	taupad = NINT(taupad / dt) * dt;
	taupad = 0.0;
	ftau = ft - taupad;
	dtau = dt;
	ntau = nt + 2.0 * taupad / dtau;

	fprintf(stderr, "nxb=%d dxb=%g npx=%d dpx=%g\n", nxb, dxb, npx, dpx);

	cell1 = (Cell ***)alloc3(mz, mx, npx, sizeof(Cell));
	memset(&cell1[0][0][0], 0, sizeof(Cell) * mx * mz * npx);

	/*rays form shot position*/
	for (ipx = 0, px = fpx + 0 * dpx; ipx < npx; ++ipx, px += dpx)
	{

		/* emergence angle and location */
		a0 = -asin(px * v[sp][0]);
		x0 = fx + (sp)*dx;
		z0 = fz;

		if (px * v[sp][0] > sin(amax) + 0.01)
			continue;
		if (px * v[sp][0] < sin(amin) - 0.01)
			continue;

		/* beam half-width adjusted for cosine of angle */
		bwhc = bwh * cos(a0);

		/* trace ray */
		ray1 = makeRay(x0, z0, a0, 2 * nt / 3, dt, ft, nx, dx, fx, nz, dz, fz, v);
		accray(ray1, cell1[ipx], fmin, bwhc, lx, lz, mx, mz, live, ipx, ntau, dtau, ftau,
			   nx, dx, fx, nz, dz, fz, g, v);
		freeRay(ray1);

#if 0
		for(imx=0;imx<mx;imx++)
		{	
			for(imz=0; imz<mz; imz++) {
				if(ipx==22) cell1[ipx][imx][imz].ti=-100000.0;

			}

		}

//		for(imx=0;imx<mx;imx++)
//		{	
//			for(imz=0; imz<mz; imz++) {
 //	   		   	fwrite(&cell1[ipx][imx][imz].tr,sizeof(float),1,cellfp);
//			}
//		}
#endif
	}

	/* loop over beam centers */
	for (ix = 0, x = fx + 0 * dtr; ix < ntr; ++ix, x += dtr)
	// for (ixb = 0, xb = fxb + 0 * dxb; ixb < nxb; ++ixb, xb += dxb)
	{

		/* horizontal window */
		// ix = NINT((xb - fx) / dtr);
		// ixlo = MAX(ix - NINT(xwh / dtr), 0);
		// ixhi = MIN(ix + NINT(xwh / dtr), ntr - 1);
		// nxw = 1 + ixhi - ixlo;
		// fxw = fx + (ixlo - ix) * dtr;

		rayxp = firp + NINT((x - fx) / dx); // ray initial point of each beam center//

		/*  Ensure initial point of ray-cdp at the range of velocity model*/
		if ((rayxp >= 0) && (rayxp < nx))
		{

			//			fprintf(stderr,"ixb/nxb = %d/%d  ix = %d ipx=%d ntau=%d\n",ixb,nxb,ixlo,npx,ntau);

			/* allocate space for beams */
			b = alloc2float(ntau, npx);
			imb = alloc2float(ntau, npx);

			/* zero b[npx][ntau] */
			memset(&b[0][0], 0, sizeof(float) * ntau * npx);

			/*if migration, form beams at surface */
			if (type == -1)
			{
				// formBeams(bwh, dxb, fmin, nt, dt, ft, nxw, dtr, fxw, &f[ixlo],
				// 		  ntau, dtau, ftau, npx, dpx, fpx, b);

				for (ipx = 0; ipx < npx; ipx++)
				{
					/* copy to output array */
					for (itau = 0; itau < ntau; itau++)
					{
						b[ipx][itau] = f[ix][itau];
					}
				}
			}

			cell2 = (Cell ***)alloc3(mz, mx, npx, sizeof(Cell));

			/* zero cell and imb */
			memset(&cell2[0][0][0], 0, sizeof(Cell) * mx * mz * npx);
			memset(&imb[0][0], 0, sizeof(float) * ntau * npx);

			/* loop over beams */
			for (ipx = 0, px = fpx + 0 * dpx; ipx < npx; ++ipx, px += dpx)
			{

				/* sine of emergence angle; skip if out of bounds */
				if (px * v[rayxp][0] > sin(amax) + 0.01)
					continue;
				if (px * v[rayxp][0] < sin(amin) - 0.01)
					continue;

				/* emergence angle and location */
				a0 = -asin(px * v[rayxp][0]);
				x0 = (rayxp)*dx;
				z0 = fz;

				/* beam half-width adjusted for cosine of angle */
				bwhc = bwh * cos(a0);

				/* trace ray */

				ray2 = makeRay(x0, z0, a0, 2 * nt / 3, dt, ft, nx, dx, fx, nz, dz, fz, v);

				/* accumulate ray to the coarse grid--cell */
				accray(ray2, cell2[ipx], fmin, bwhc, lx, lz, mx, mz, live, ipx, ntau, dtau, ftau,
					   nx, dx, fx, nz, dz, fz, g, v);

#if 0			
				/* write out ray steps */	
				for(irs=0; irs<ray2->nrs; irs++) {
					//printf("%d %d\n",irs,ray2->nrs);
					fwrite(&ray2->rs[irs].x,sizeof(float),1,rayxfp);
					fwrite(&ray2->rs[irs].z,sizeof(float),1,rayzfp);
				}
#endif

				/* free ray */
				freeRay(ray2);

#if 0
				for(imx=0;imx<mx;imx++)
				{	
					for(imz=0; imz<mz; imz++) {
						if(ipx==1110) cell2[ipx][imx][imz].ti=-100000.0;
					}
				}

				for(imx=0;imx<mx;imx++)
				{	
					for(imz=0; imz<mz; imz++) {
 	   		   			fwrite(&cell2[ipx][imx][imz].ti,sizeof(float),1,cellfp);
					}
				}
#endif

			} // end loop of ipx

			/* main function of model or migration using beams */
			scanimg(type, cell1, cell2, live, dead, ntau, dtau, ftau, fmin, lx, lz, nx, nz, mx, mz, npx, b, g, fpx, dpx, imb);

#if 0
			for(ixx=0;ixx<npx;ixx++)
			{	
 	   		   		fwrite(imb[ixx],sizeof(float),ntau,cellfp);
			}
#endif
			/* if MOD, calculate the filtered beams at beam centers*/
			if (type == 1)
			{
				convol(b, imb, npx, ntau, ftau, dtau);
			}

#if 0
			for(ixx=0;ixx<npx;ixx++)
			{	
 	   		   		fwrite(b[ixx],sizeof(float),ntau,cellfp);
			}
#endif
			/* if MOD, calculate the synthetic data */
			if (type == 1)
			{
				synData(bwh, dxb, fmin, nt, dt, ft, nxw, dtr, fxw, &f[ixlo], ntau, dtau, ftau, npx, dpx, fpx, b);
			}

#if 0
			if(type==-1) {
				for(ixx=0;ixx<npx;ixx++)
				{	
 	 	  	   		fwrite(b[ixx],sizeof(float),ntau,cellfp);
				}
			}
#endif

			free2float(b);
			free2float(imb);
			free3((void ***)cell2);

		} // ipx
	}	  // end loop of ibx

	free3((void ***)cell1);

#if 0
	if(type==1) {
		for(ixx=0;ixx<nx;ixx++)
		{	
 	   	   		fwrite(f[ixx],sizeof(float),ntau,cellfp);
		}
	}
#endif
}

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

Ray *makeRay(float x0, float z0, float a0, int nt, float dt, float ft,
			 int nx, float dx, float fx, int nz, float dz, float fz, float **vxz)
/*****************************************************************************
Trace a ray for uniformly sampled v(x,z).
******************************************************************************
	  nput:
x0		x coordinate of takeoff point
z0		z coordinate of takeoff point
a0		takeoff angle (radians)
nt		number of time samples
dt		time sampling interval
ft		first time sample
nx		number of x samples
dx		x sampling interval
fx		first x sample
nz		number of z samples
dz		z sampling interval
fz		first z sample
vxz		array[nx][nz] of uniformly sampled velocities v(x,z)

Returned:	pointer to ray parameters sampled at discrete ray steps
******************************************************************************
Notes:
The ray ends when it runs out of time (after nt steps) or with the first
step that is out of the (x,z) bounds of the velocity function v(x,z).
*****************************************************************************/
{
	int it, kmah;
	float t, x, z, a, c, s, p1, q1, p2, q2,
		lx, lz, cc, ss,
		v, dvdx, dvdz, ddvdxdx, ddvdxdz, ddvdzdz,
		vv, ddvdndn;
	Ray *ray;
	RayStep *rs;

	void *vel2;

	/* allocate and initialize velocity interpolator */
	vel2 = vel2Alloc(nx, dx, fx, nz, dz, fz, vxz);

	/* last x and z in velocity model */
	lx = fx + (nx - 1) * dx;
	lz = fz + (nz - 1) * dz;

	/* ensure takeoff point is within model */
	if (x0 < fx || x0 > lx || z0 < fz || z0 > lz)
		return NULL;

	/* allocate space for ray and raysteps */
	ray = (Ray *)alloc1(1, sizeof(Ray));
	rs = (RayStep *)alloc1(nt, sizeof(RayStep));

	/* cosine and sine of takeoff angle */
	c = cos(a0);
	s = sin(a0);
	cc = c * c;
	ss = s * s;

	/* velocity and derivatives at takeoff point */
	vel2Interp(vel2, x0, z0, &v, &dvdx, &dvdz, &ddvdxdx, &ddvdxdz, &ddvdzdz);
	ddvdndn = 0;
	vv = v * v;

	/* first ray step */
	rs[0].t = t = ft;
	rs[0].a = a = a0;
	rs[0].x = x = x0;
	rs[0].z = z = z0;
	rs[0].q1 = q1 = 1.0;
	rs[0].p1 = p1 = 0.0;
	rs[0].q2 = q2 = 0.0;
	rs[0].p2 = p2 = 1.0;
	rs[0].kmah = kmah = 0;
	rs[0].c = c;
	rs[0].s = s;
	rs[0].v = v;
	rs[0].dvdx = dvdx;
	rs[0].dvdz = dvdz;

	/* loop over time steps */
	for (it = 1; it < nt; ++it)
	{

		/* variables used for Runge-Kutta integration */
		float h = dt, hhalf = dt / 2.0, hsixth = dt / 6.0,
			  q2old, xt, zt, at, p1t, q1t, p2t, q2t,
			  dx, dz, da, dp1, dq1, dp2, dq2,
			  dxt, dzt, dat, dp1t, dq1t, dp2t, dq2t,
			  dxm, dzm, dam, dp1m, dq1m, dp2m, dq2m;

		/* if ray is out of bounds, break */
		if (x < fx || x > lx || z < fz || z > lz)
			break;

		/* remember old q2 */
		q2old = q2;

		/* step 1 of 4th-order Runge-Kutta */
		dx = v * s;
		dz = v * c;
		da = dvdz * s - dvdx * c;
		dp1 = -ddvdndn * q1 / v;
		dq1 = vv * p1;
		dp2 = -ddvdndn * q2 / v;
		dq2 = vv * p2;
		xt = x + hhalf * dx;
		zt = z + hhalf * dz;
		at = a + hhalf * da;
		p1t = p1 + hhalf * dp1;
		q1t = q1 + hhalf * dq1;
		p2t = p2 + hhalf * dp2;
		q2t = q2 + hhalf * dq2;
		c = cos(at);
		s = sin(at);
		cc = c * c;
		ss = s * s;
		vel2Interp(vel2, xt, zt,
				   &v, &dvdx, &dvdz, &ddvdxdx, &ddvdxdz, &ddvdzdz);
		ddvdndn = 0;
		vv = v * v;

		/* step 2 of 4th-order Runge-Kutta */
		dxt = v * s;
		dzt = v * c;
		dat = dvdz * s - dvdx * c;
		dp1t = -ddvdndn * q1t / v;
		dq1t = vv * p1t;
		dp2t = -ddvdndn * q2t / v;
		dq2t = vv * p2t;
		xt = x + hhalf * dxt;
		zt = z + hhalf * dzt;
		at = a + hhalf * dat;
		p1t = p1 + hhalf * dp1t;
		q1t = q1 + hhalf * dq1t;
		p2t = p2 + hhalf * dp2t;
		q2t = q2 + hhalf * dq2t;
		c = cos(at);
		s = sin(at);
		cc = c * c;
		ss = s * s;
		vel2Interp(vel2, xt, zt,
				   &v, &dvdx, &dvdz, &ddvdxdx, &ddvdxdz, &ddvdzdz);
		ddvdndn = 0;
		vv = v * v;

		/* step 3 of 4th-order Runge-Kutta */
		dxm = v * s;
		dzm = v * c;
		dam = dvdz * s - dvdx * c;
		dp1m = -ddvdndn * q1t / v;
		dq1m = vv * p1t;
		dp2m = -ddvdndn * q2t / v;
		dq2m = vv * p2t;
		xt = x + h * dxm;
		zt = z + h * dzm;
		at = a + h * dam;
		p1t = p1 + h * dp1m;
		q1t = q1 + h * dq1m;
		p2t = p2 + h * dp2m;
		q2t = q2 + h * dq2m;
		dxm += dxt;
		dzm += dzt;
		dam += dat;
		dp1m += dp1t;
		dq1m += dq1t;
		dp2m += dp2t;
		dq2m += dq2t;
		c = cos(at);
		s = sin(at);
		cc = c * c;
		ss = s * s;
		vel2Interp(vel2, xt, zt,
				   &v, &dvdx, &dvdz, &ddvdxdx, &ddvdxdz, &ddvdzdz);
		ddvdndn = 0;
		vv = v * v;

		/* step 4 of 4th-order Runge-Kutta */
		dxt = v * s;
		dzt = v * c;
		dat = dvdz * s - dvdx * c;
		dp1t = -ddvdndn * q1t / v;
		dq1t = vv * p1t;
		dp2t = -ddvdndn * q2t / v;
		dq2t = vv * p2t;
		x += hsixth * (dx + dxt + 2.0 * dxm);
		z += hsixth * (dz + dzt + 2.0 * dzm);
		a += hsixth * (da + dat + 2.0 * dam);
		p1 += hsixth * (dp1 + dp1t + 2.0 * dp1m);
		q1 += hsixth * (dq1 + dq1t + 2.0 * dq1m);
		p2 += hsixth * (dp2 + dp2t + 2.0 * dp2m);
		q2 += hsixth * (dq2 + dq2t + 2.0 * dq2m);
		c = cos(a);
		s = sin(a);
		cc = c * c;
		ss = s * s;
		vel2Interp(vel2, x, z,
				   &v, &dvdx, &dvdz, &ddvdxdx, &ddvdxdz, &ddvdzdz);
		ddvdndn = 0;
		vv = v * v;

		/* update kmah index */
		if ((q2 <= 0.0 && q2old > 0.0) || (q2 >= 0.0 && q2old < 0.0))
			kmah++;

		/* control the rays */
		/* (a-a0), the difference between initial ray angle and current ray angle */
		if (fabs(a) > 90 * PI / 180 || fabs(a - a0) > PI / 3)
			break;
		// if(fabs(a)>90*PI/180) break;

		/* update time */
		t += dt;

		/* save ray parameters */
		rs[it].t = t;
		rs[it].a = a;
		rs[it].x = x;
		rs[it].z = z;
		rs[it].c = c;
		rs[it].s = s;
		rs[it].q1 = q1;
		rs[it].p1 = p1;
		rs[it].q2 = q2;
		rs[it].p2 = p2;
		rs[it].kmah = kmah;
		rs[it].v = v;
		rs[it].dvdx = dvdx;
		rs[it].dvdz = dvdz;
	}

	/* free velocity interpolator */
	vel2Free(vel2);

	/* return ray */
	ray->nrs = it;
	ray->rs = rs;
	ray->nc = 0;
	ray->c = NULL;
	return ray;
}

/*conjugate-direction descent*/
void cgstep(int iter, float **x, float **g, float **rr, float **gg, int nx, int nz, int ntr, int nt, float **s, float **ss)
{
	double sds = 0.0, gdg = 0.0, gds = 0.0, determ = 0.0, gdr = 0.0, sdr = 0.0, alfa = 0.0, beta = 0.0;
	double dotgg = 0.0, dsum, esum;

	int i, j;

	if (iter == 1)
	{

		memset(&s[0][0], 0, sizeof(float) * nz * nx);
		memset(&ss[0][0], 0, sizeof(float) * nt * ntr);
		beta = 0.0;

		for (i = 0; i < ntr; i++)
			for (j = 0; j < nt; j++)
				dotgg += gg[i][j] * gg[i][j];

		if (dotgg == 0.0)
			printf("Error happens beacause grad vanishes identically");

		dsum = esum = 0.0;
		for (i = 0; i < ntr; i++)
			for (j = 0; j < nt; j++)
				dsum += rr[i][j] * gg[i][j];

		for (i = 0; i < ntr; i++)
			for (j = 0; j < nt; j++)
				esum += gg[i][j] * gg[i][j];

		alfa = -dsum / esum;
	}
	else
	{
		for (i = 0; i < ntr; i++)
			for (j = 0; j < nt; j++)
				gdg += gg[i][j] * gg[i][j],
					sds += ss[i][j] * ss[i][j],
					gds += gg[i][j] * ss[i][j];

		if ((gdg == 0.0) || (sds == 0.0))
			printf("gdg = 0.0, or sds = 0.0\n");

		determ = gdg * sds * MAX(1.0 - (gds / gdg) * (gds / sds), 1.0e-12);

		//		printf("%lf %d\n",determ,1230);

		for (i = 0; i < ntr; i++)
			for (j = 0; j < nt; j++)
				gdr -= gg[i][j] * rr[i][j],
					sdr -= ss[i][j] * rr[i][j];
		alfa = (sds * gdr - gds * sdr) / determ,
		beta = (-gds * gdr + gdg * sdr) / determ;
	}

	for (i = 0; i < nx; i++)
		for (j = 0; j < nz; j++)
			s[i][j] = alfa * g[i][j] + beta * s[i][j],
			x[i][j] += s[i][j];

	for (i = 0; i < ntr; i++)
		for (j = 0; j < nt; j++)
			ss[i][j] = alfa * gg[i][j] + beta * ss[i][j],
			rr[i][j] += ss[i][j];
}

/******************************************************************************/
/* Input the data of shot gather and determine source and first receiver point*/
/******************************************************************************/
void inputrace(int is, int nt, float dx, int maxtr, FILE *fp, int *sisp, int *firisp, int *nistr)
{

	int itr, i, ntrace; // ntrΪ\B5\DAis\C5ڵĵ\C0\CA\FD
	int nshot;
	int *head;
	int offsetmin; //\D7\EEСƫ\D2ƾ\E0
	int cdpmin;	   //\D7\EEСƫ\D2ƾ\E0\B6\D4Ӧ\B5\C4CDP\A3\AC\D3\C3\D2\D4\C7\F3\B3\F6\C5ڵ\E3CDP
	offsetmin = 10000;
	head = alloc1int(60);

	for (i = 0; i < 60; i++)
		head[i] = 0;
	rewind(fp);
	for (itr = 0; itr < maxtr; itr++)
	{

		fread(head, sizeof(int), 60, fp);
		fseek(fp, 4 * nt, 1);

		if (head[2] == is)
		{
			fseek(fp, -4 * (nt + 60), 1);
			break;
		}
	}

	if (head[2] != is)
	{
		*nistr = 0;
		return;
	}
	ntrace = 0;
	for (i = 0; i < 60; i++)
		head[i] = 0;
	for (i = itr; i < maxtr; i++)
	{
		fread(head, sizeof(int), 60, fp);
		fseek(fp, 4 * nt, 1);
		if (abs(head[9]) < abs(offsetmin))
		{
			offsetmin = head[9];
			cdpmin = head[5];
		}
		if (head[2] == is)
			ntrace++;
		else
			break;
	}
	*nistr = ntrace;
	// printf("%d\n",*nistr);

	if (offsetmin < 0) //\C8\F4offsetmin<0,\D4\F2\B6϶\A8Ϊ\D3Ҷ˷\C5\C5ڻ\F2\D6м\E4\B7\C5\C5\DA
		*sisp = cdpmin + NINT(abs(offsetmin) / 2 / dx);
	else //\B7\F1\D4\F2Ϊ\D7\F3\B6˷\C5\C5ڻ\F2\D6м\E4\B7\C5\C5\DA
		*sisp = cdpmin - NINT(abs(offsetmin) / 2 / dx);
	//    printf("%d\n",*sisp);
	if (i == maxtr)
		fseek(fp, -4 * (nt + 60) * (*nistr), 1); //\BBص\BD\B8\C3\C5ڵ\C4\CEļ\FEͷ\A3\AC\BD\D3\CF\C2\C0\B4\B6\C1\C8\EB\C5ڼ\C7¼
	else
		fseek(fp, -4 * (nt + 60) * (*nistr + 1), 1);
	fread(head, sizeof(int), 60, fp);
	*firisp = 2 * head[5] - *sisp;
	fseek(fp, -4 * 60, 1);
	free1int(head);
}

void partall(int type, int nz, int apernum, int apermin, float **part, float **all)
{
	int ix, iz;
	if (type == 1)
	{
		for (ix = 0; ix < apernum; ix++)
			for (iz = 0; iz < nz; iz++)
				part[ix][iz] = all[apermin + ix][iz];
	}
	else if (type == -1)
	{
		for (ix = 0; ix < apernum; ix++)
			for (iz = 0; iz < nz; iz++)
				all[apermin + ix][iz] += part[ix][iz];
	}
}

void freeRay(Ray *ray)
/*****************************************************************************
Free a ray.
******************************************************************************
	  nput:
ray		ray to be freed
*****************************************************************************/
{
	if (ray->c != NULL)
		free1((void *)ray->c);
	free1((void *)ray->rs);
	free1((void *)ray);
}

void smooth2d(int n1, int n2, float r1, float r2, float **v)
{

	int nmax, ix, iz;

	float **w;
	float *d, *e, *f;
	int *win;
	float rw = 0.0;

	r1 = r1 * r1 * 0.25;
	r2 = r2 * r2 * 0.25;

	/* allocate space */
	nmax = (n1 < n2) ? n2 : n1;
	win = alloc1int(4);
	//	v = alloc2float(n1,n2);
	//	v0 = alloc2float(n1,n2);
	w = alloc2float(n1, n2);
	d = alloc1float(nmax);
	e = alloc1float(nmax);
	f = alloc1float(nmax);
	rw = rw * rw * 0.25;

	win[0] = 0;
	win[1] = n1;
	win[2] = 0;
	win[3] = n2;

	/* define the window function */
	for (ix = 0; ix < n2; ++ix)
		for (iz = 0; iz < n1; ++iz)
			w[ix][iz] = 0;
	for (ix = win[2]; ix < win[3]; ++ix)
		for (iz = win[0]; iz < win[1]; ++iz)
			w[ix][iz] = 1;

	if (win[0] > 0 || win[1] < n1 || win[2] > 0 || win[3] < n2)
	{
		/*	smooth the window function */
		for (iz = 0; iz < n1; ++iz)
		{
			for (ix = 0; ix < n2; ++ix)
			{
				d[ix] = 1.0 + 2.0 * rw;
				e[ix] = -rw;
				f[ix] = w[ix][iz];
			}
			d[0] -= rw;
			d[n2 - 1] -= rw;
			tripd(d, e, f, n2);
			for (ix = 0; ix < n2; ++ix)
				w[ix][iz] = f[ix];
		}
		for (ix = 0; ix < n2; ++ix)
		{
			for (iz = 0; iz < n1; ++iz)
			{
				d[iz] = 1.0 + 2.0 * rw;
				e[iz] = -rw;
				f[iz] = w[ix][iz];
			}
			d[0] -= rw;
			d[n1 - 1] -= rw;
			tripd(d, e, f, n1);
			for (iz = 0; iz < n1; ++iz)
				w[ix][iz] = f[iz];
		}
	}

	/* solving for the smoothing velocity */
	for (iz = 0; iz < n1; ++iz)
	{
		for (ix = 0; ix < n2 - 1; ++ix)
		{
			d[ix] = 1.0 + r2 * (w[ix][iz] + w[ix + 1][iz]);
			e[ix] = -r2 * w[ix + 1][iz];
			f[ix] = v[ix][iz];
		}
		d[0] -= r2 * w[0][iz];
		d[n2 - 1] = 1.0 + r2 * w[n2 - 1][iz];
		f[n2 - 1] = v[n2 - 1][iz];
		tripd(d, e, f, n2);
		for (ix = 0; ix < n2; ++ix)
			v[ix][iz] = f[ix];
	}
	for (ix = 0; ix < n2; ++ix)
	{
		for (iz = 0; iz < n1 - 2; ++iz)
		{
			d[iz] = 1.0 + r1 * (w[ix][iz + 1] + w[ix][iz + 2]);
			e[iz] = -r1 * w[ix][iz + 2];
			f[iz] = v[ix][iz + 1];
		}
		f[0] += r1 * w[ix][1] * v[ix][0];
		d[n1 - 2] = 1.0 + r1 * w[ix][n1 - 1];
		f[n1 - 2] = v[ix][n1 - 1];
		tripd(d, e, f, n1 - 1);
		for (iz = 0; iz < n1 - 1; ++iz)
			v[ix][iz + 1] = f[iz];
	}

	free2float(w);
	free1float(d);
	free1float(e);
	free1float(f);
	free1int(win);
}

void tripd(float *d, float *e, float *b, int n)
/*****************************************************************************
Given an n-by-n symmetric, tridiagonal, positive definite matrix A and
 n-vector b, following algorithm overwrites b with the solution to Ax = b*/
{
	int k;
	float temp;

	/* decomposition */
	for (k = 1; k < n; ++k)
	{
		temp = e[k - 1];
		e[k - 1] = temp / d[k - 1];
		d[k] -= temp * e[k - 1];
	}

	/* substitution	*/
	for (k = 1; k < n; ++k)
		b[k] -= e[k - 1] * b[k - 1];

	b[n - 1] /= d[n - 1];
	for (k = n - 1; k > 0; --k)
		b[k - 1] = b[k - 1] / d[k - 1] - e[k - 1] * b[k];
}

int nearestRayStep(Ray *ray, float x, float z)
/*****************************************************************************
Determine index of ray step nearest to point (x,z).
******************************************************************************
	  nput:
ray		ray
x		x coordinate
z		z coordinate

Returned:	index of nearest ray step
*****************************************************************************/
{
	int nrs = ray->nrs, ic = ray->ic, nc = ray->nc;
	RayStep *rs = ray->rs;
	Circle *c = (Circle *)ray->c;
	int irs, irsf, irsl, irsmin = 0, update, jc, js, kc;
	float dsmin, ds, dx, dz, dmin, rdmin, xrs, zrs;

	/* if necessary, make circles localizing ray steps */
	if (c == NULL)
	{
		ray->ic = ic = 0;
		ray->nc = nc = sqrt((float)nrs);
		ray->c = c = makeCircles(nc, nrs, rs);
	}

	/* initialize minimum distance and minimum distance-squared */
	dx = x - c[ic].x;
	dz = z - c[ic].z;
	dmin = 2.0 * (sqrt(dx * dx + dz * dz) + c[ic].r);
	dsmin = dmin * dmin;

	/* loop over all circles */
	for (kc = 0, jc = ic, js = 0; kc < nc; ++kc)
	{

		/* distance-squared to center of circle */
		dx = x - c[jc].x;
		dz = z - c[jc].z;
		ds = dx * dx + dz * dz;

		/* radius of circle plus minimum distance (so far) */
		rdmin = c[jc].r + dmin;

		/* if circle could possible contain a nearer ray step */
		if (ds <= rdmin * rdmin)
		{

			/* search circle for nearest ray step */
			irsf = c[jc].irsf;
			irsl = c[jc].irsl;
			update = 0;
			for (irs = irsf; irs <= irsl; ++irs)
			{
				xrs = rs[irs].x;
				zrs = rs[irs].z;
				dx = x - xrs;
				dz = z - zrs;
				ds = dx * dx + dz * dz;
				if (ds < dsmin)
				{
					dsmin = ds;
					irsmin = irs;
					update = 1;
				}
			}

			/* if a nearer ray step was found inside circle */
			if (update)
			{

				/* update minimum distance */
				dmin = sqrt(dsmin);

				/* remember the circle */
				ic = jc;
			}
		}

		/* search circles in alternating directions */
		js = (js > 0) ? -js - 1 : -js + 1;
		jc += js;
		if (jc < 0 || jc >= nc)
		{
			js = (js > 0) ? -js - 1 : -js + 1;
			jc += js;
		}
	}

	/* remember the circle containing the nearest ray step */
	ray->ic = ic;

	if (irsmin < 0 || irsmin >= nrs)
		fprintf(stderr, "irsmin=%d\n", irsmin);

	/* return index of nearest ray step */
	return irsmin;
}

int xxx_nearestRayStep(Ray *ray, float x, float z)
/*****************************************************************************
Determine index of ray step nearest to point (x,z).  Simple (slow) version.
******************************************************************************
	  nput:
ray		ray
x		x coordinate
z		z coordinate

Returned:	index of nearest ray step
*****************************************************************************/
{
	int nrs = ray->nrs;
	RayStep *rs = ray->rs;
	int irs, irsmin = 0;
	float dsmin, ds, dx, dz, xrs, zrs;

	for (irs = 0, dsmin = FLT_MAX; irs < nrs; ++irs)
	{
		xrs = rs[irs].x;
		zrs = rs[irs].z;
		dx = x - xrs;
		dz = z - zrs;
		ds = dx * dx + dz * dz;
		if (ds < dsmin)
		{
			dsmin = ds;
			irsmin = irs;
		}
	}
	return irsmin;
}

Circle *makeCircles(int nc, int nrs, RayStep *rs)
/*****************************************************************************
Make circles used to speed up determination of nearest ray step.
******************************************************************************
	  nput:
nc		number of circles to make
nrs		number of ray steps
rs		array[nrs] of ray steps

Returned:	array[nc] of circles
*****************************************************************************/
{
	int nrsc, ic, irsf, irsl, irs;
	float xmin, xmax, zmin, zmax, x, z, r;
	Circle *c;

	/* allocate space for circles */
	c = (Circle *)alloc1(nc, sizeof(Circle));

	/* determine typical number of ray steps per circle */
	nrsc = 1 + (nrs - 1) / nc;

	/* loop over circles */
	for (ic = 0; ic < nc; ++ic)
	{

		/* index of first and last raystep */
		irsf = ic * nrsc;
		irsl = irsf + nrsc - 1;
		if (irsf >= nrs)
			irsf = nrs - 1;
		if (irsl >= nrs)
			irsl = nrs - 1;

		/* coordinate bounds of ray steps */
		xmin = xmax = rs[irsf].x;
		zmin = zmax = rs[irsf].z;
		for (irs = irsf + 1; irs <= irsl; ++irs)
		{
			if (rs[irs].x < xmin)
				xmin = rs[irs].x;
			if (rs[irs].x > xmax)
				xmax = rs[irs].x;
			if (rs[irs].z < zmin)
				zmin = rs[irs].z;
			if (rs[irs].z > zmax)
				zmax = rs[irs].z;
		}

		/* center and radius of circle */
		x = 0.5 * (xmin + xmax);
		z = 0.5 * (zmin + zmax);
		r = sqrt((x - xmin) * (x - xmin) + (z - zmin) * (z - zmin));

		/* set circle */
		c[ic].irsf = irsf;
		c[ic].irsl = irsl;
		c[ic].x = x;
		c[ic].z = z;
		c[ic].r = r;
	}

	return c;
}

/*****************************************************************************
Functions to support interpolation of velocity and its derivatives.
******************************************************************************
Functions:
vel2Alloc	allocate and initialize an interpolator for v(x,z)
vel2Interp	interpolate v(x,z) and its derivatives
******************************************************************************
Notes:
	  nterpolation is performed by piecewise cubic Hermite polynomials,
so that velocity and first derivatives are continuous.  Therefore,
velocity v, first derivatives dv/dx and dv/dz, and the mixed
derivative ddv/dxdz are continuous.  However, second derivatives
ddv/dxdx and ddv/dzdz are discontinuous.
*****************************************************************************/

/* number of pre-computed, tabulated interpolators */
#define NTABLE 101

/* length of each interpolator in table (4 for piecewise cubic) */
#define LTABLE 4

/* table of pre-computed interpolators, for 0th, 1st, and 2nd derivatives */
static float tbl[3][NTABLE][LTABLE];

/* constants */
static int ix = 1 - LTABLE / 2 - LTABLE, iz = 1 - LTABLE / 2 - LTABLE;
static float ltable = LTABLE, ntblm1 = NTABLE - 1;

/* indices for 0th, 1st, and 2nd derivatives */
static int kx[6] = {0, 1, 0, 2, 1, 0};
static int kz[6] = {0, 0, 1, 0, 1, 2};

/* function to build interpolator tables; sets tabled=1 when built */
static void buildTables(void);
static int tabled = 0;

/* interpolator for velocity function v(x,z) of two variables */
typedef struct Vel2Struct
{
	int nx;	 /* number of x samples */
	int nz;	 /* number of z samples */
	int nxm; /* number of x samples minus LTABLE */
	int nzm; /* number of x samples minus LTABLE */
	float xs, xb, zs, zb, sx[3], sz[3], **vu;
} Vel2;

void *vel2Alloc(int nx, float dx, float fx,
				int nz, float dz, float fz, float **v)
/*****************************************************************************
Allocate and initialize an interpolator for v(x,z) and its derivatives.
******************************************************************************
	  nput:
nx		number of x samples
dx		x sampling interval
fx		first x sample
nz		number of z samples
dz		z sampling interval
fz		first z sample
v		array[nx][nz] of uniformly sampled v(x,z)

Returned:	pointer to interpolator
*****************************************************************************/
{
	Vel2 *vel2;

	/* allocate space */
	vel2 = (Vel2 *)alloc1(1, sizeof(Vel2));

	/* set state variables used for interpolation */
	vel2->nx = nx;
	vel2->nxm = nx - LTABLE;
	vel2->xs = 1.0 / dx;
	vel2->xb = ltable - fx * vel2->xs;
	vel2->sx[0] = 1.0;
	vel2->sx[1] = vel2->xs;
	vel2->sx[2] = vel2->xs * vel2->xs;
	vel2->nz = nz;
	vel2->nzm = nz - LTABLE;
	vel2->zs = 1.0 / dz;
	vel2->zb = ltable - fz * vel2->zs;
	vel2->sz[0] = 1.0;
	vel2->sz[1] = vel2->zs;
	vel2->sz[2] = vel2->zs * vel2->zs;
	vel2->vu = v;

	/* if necessary, build interpolator coefficient tables */
	if (!tabled)
		buildTables();

	return vel2;
}

void vel2Free(void *vel2)
/*****************************************************************************
Free an interpolator for v(x,z) and its derivatives.
******************************************************************************
	  nput:
vel2		pointer to interpolator as returned by vel2Alloc()
*****************************************************************************/
{
	free1(vel2);
}

void vel2Interp(void *vel2, float x, float z,
				float *v, float *vx, float *vz, float *vxx, float *vxz, float *vzz)
/*****************************************************************************
	  nterpolation of a velocity function v(x,z) and its derivatives.
******************************************************************************
	  nput:
vel2		pointer to interpolator as returned by vel2Alloc()
x		x coordinate at which to interpolate v(x,z) and derivatives
z		z coordinate at which to interpolate v(x,z) and derivatives

Output:
v		v(x,z)
vx		dv/dx
vz		dv/dz
vxx		ddv/dxdx
vxz		ddv/dxdz
vzz		ddv/dzdz
*****************************************************************************/
{
	Vel2 *v2 = (Vel2 *)vel2;
	int nx = v2->nx, nz = v2->nz, nxm = v2->nxm, nzm = v2->nzm;
	float xs = v2->xs, xb = v2->xb, zs = v2->zs, zb = v2->zb,
		  *sx = v2->sx, *sz = v2->sz, **vu = v2->vu;
	int i, jx, lx, mx, jz, lz, mz, jmx, jmz, mmx, mmz;
	float ax, bx, *px, az, bz, *pz, sum, vd[6];

	/* determine offsets into vu and interpolation coefficients */
	ax = xb + x * xs;
	jx = (int)ax;
	bx = ax - jx;
	lx = (bx >= 0.0) ? bx * ntblm1 + 0.5 : (bx + 1.0) * ntblm1 - 0.5;
	lx *= LTABLE;
	mx = ix + jx;
	az = zb + z * zs;
	jz = (int)az;
	bz = az - jz;
	lz = (bz >= 0.0) ? bz * ntblm1 + 0.5 : (bz + 1.0) * ntblm1 - 0.5;
	lz *= LTABLE;
	mz = iz + jz;

	/* if totally within input array, use fast method */
	if (mx >= 0 && mx <= nxm && mz >= 0 && mz <= nzm)
	{
		for (i = 0; i < 6; ++i)
		{
			px = &(tbl[kx[i]][0][0]) + lx;
			pz = &(tbl[kz[i]][0][0]) + lz;
			vd[i] = sx[kx[i]] * sz[kz[i]] * (vu[mx][mz] * px[0] * pz[0] + vu[mx][mz + 1] * px[0] * pz[1] + vu[mx][mz + 2] * px[0] * pz[2] + vu[mx][mz + 3] * px[0] * pz[3] + vu[mx + 1][mz] * px[1] * pz[0] + vu[mx + 1][mz + 1] * px[1] * pz[1] + vu[mx + 1][mz + 2] * px[1] * pz[2] + vu[mx + 1][mz + 3] * px[1] * pz[3] + vu[mx + 2][mz] * px[2] * pz[0] + vu[mx + 2][mz + 1] * px[2] * pz[1] + vu[mx + 2][mz + 2] * px[2] * pz[2] + vu[mx + 2][mz + 3] * px[2] * pz[3] + vu[mx + 3][mz] * px[3] * pz[0] + vu[mx + 3][mz + 1] * px[3] * pz[1] + vu[mx + 3][mz + 2] * px[3] * pz[2] + vu[mx + 3][mz + 3] * px[3] * pz[3]);
		}

		/* else handle end effects with constant extrapolation */
	}
	else
	{
		for (i = 0; i < 6; ++i)
		{
			px = &(tbl[kx[i]][0][0]) + lx;
			pz = &(tbl[kz[i]][0][0]) + lz;
			for (jx = 0, jmx = mx, sum = 0.0; jx < 4; ++jx, ++jmx)
			{
				mmx = jmx;
				if (mmx < 0)
					mmx = 0;
				else if (mmx >= nx)
					mmx = nx - 1;
				for (jz = 0, jmz = mz; jz < 4; ++jz, ++jmz)
				{
					mmz = jmz;
					if (mmz < 0)
						mmz = 0;
					else if (mmz >= nz)
						mmz = nz - 1;
					sum += vu[mmx][mmz] * px[jx] * pz[jz];
				}
			}
			vd[i] = sx[kx[i]] * sz[kz[i]] * sum;
		}
	}

	/* set output variables */
	*v = vd[0];
	*vx = vd[1];
	*vz = vd[2];
	*vxx = vd[3];
	*vxz = vd[4];
	*vzz = vd[5];
}

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

/* function to build interpolation tables */
static void buildTables(void)
{
	int i;
	float x;

	/* tabulate interpolator for 0th derivative */
	for (i = 0; i < NTABLE; ++i)
	{
		x = (float)i / (NTABLE - 1.0);
		tbl[0][i][0] = -0.5 * k00(x);
		tbl[0][i][1] = h00(x) - 0.5 * k10(x);
		tbl[0][i][2] = h10(x) + 0.5 * k00(x);
		tbl[0][i][3] = 0.5 * k10(x);
		tbl[1][i][0] = -0.5 * k01(x);
		tbl[1][i][1] = h01(x) - 0.5 * k11(x);
		tbl[1][i][2] = h11(x) + 0.5 * k01(x);
		tbl[1][i][3] = 0.5 * k11(x);
		tbl[2][i][0] = -0.5 * k02(x);
		tbl[2][i][1] = h02(x) - 0.5 * k12(x);
		tbl[2][i][2] = h12(x) + 0.5 * k02(x);
		tbl[2][i][3] = 0.5 * k12(x);
	}

	/* remember that tables have been built */
	tabled = 1;
}

/* Beam subroutines */

/* functions defined and used internally */
static void xtop(float w,
				 int nx, float dx, float fx, complex *g,
				 int np, float dp, float fp, complex *h);
static BeamData *beamData(int type, int npx, float wmin, int nt, float dt, float ft, float **f);
static void setCell(Cells *cells, int jx, int jz);
static void accCell(Cells *cells, int jx, int jz);
static int cellTimeAmp(Cells *cells, int jx, int jz);
static void cellBeam(int type, Cell **cell, float **f, float **g, int jx, int jz, BeamData *bd, int lx, int lz,
					 int nx, int nz, int npx, int ips, int ipr, float fpx, float dpx,
					 float ft, float dt, int nt, float **imb);

/* functions for external use */

void convol(float **f, float **imb, int npx, int ntau, float ftau, float dtau)
{
	int ntfft, nw, ix, iw, ipx, it, itau;
	float wmin, x, dw, fw, w, fftscl,
		amp, phase, scale, a, b, as, bs, es, cfr, cfi,
		*fpad, *gpad, *imfpad, *imgpad;
	complex **cf, **imcf, *cb;
	int nf1, nf2, nf3, nf4;
	float ffw, tmpp;
	float spr;
	float wd = 2.0 * PI * DOMINFY;
	int iwnyq;

	/* fft sampling */
	ntfft = npfar(ntau + 100);
	nw = ntfft / 2 + 1;
	dw = 2.0 * PI / (ntfft * dtau);
	fw = 0.0;
	fftscl = 1.0 / ntfft;
	iwnyq = nw - 1;

	/* allocate workspace */
	fpad = alloc1float(ntfft);
	imfpad = alloc1float(ntfft);
	gpad = alloc1float(ntfft);
	imgpad = alloc1float(ntfft);
	cf = alloc2complex(nw, npx);
	imcf = alloc2complex(nw, npx);
	cb = alloc1complex(ntfft);

	//	smooth2d(ntau, npx, 5.0, 0.0, f);

	/* loop over x */
	for (ipx = 0; ipx < npx; ++ipx)
	{

		/* pad time with zeros */
		for (it = 0; it < ntau; ++it)
			fpad[it] = f[ipx][it],
			imfpad[it] = imb[ipx][it];
		for (it = ntau; it < ntfft; ++it)
			fpad[it] = 0.0,
			imfpad[it] = 0.0;

		/* Fourier transform time to frequency */
		pfarc(1, ntfft, fpad, cf[ipx]);
		pfarc(1, ntfft, imfpad, imcf[ipx]);

		for (iw = 1; iw < iwnyq; iw++)
			imcf[ipx][iw].r *= 2.0,
				imcf[ipx][iw].i *= 2.0;
	}

	for (ipx = 0; ipx < npx; ++ipx)
	{
		for (iw = 0, w = fw; iw < nw; iw++, w += dw)
		{
			spr = 4 * w * w * sqrt(PI) * exp(-w * w / wd / wd) / wd / wd / wd;
			// spr = 1.0;
			cf[ipx][iw].r = cf[ipx][iw].r * w * w * spr * fftscl;
			cf[ipx][iw].i = cf[ipx][iw].i * w * w * spr * fftscl;
			imcf[ipx][iw].r = imcf[ipx][iw].r * w * w * spr * fftscl;
			imcf[ipx][iw].i = imcf[ipx][iw].i * w * w * spr * fftscl;
		}
	}

	/* loop over px */
	for (ipx = 0; ipx < npx; ++ipx)
	{

		for (iw = 0, w = fw; iw < nw; iw++, w += dw)
		{
			cb[iw].r = imcf[ipx][iw].r;
			cb[iw].i = imcf[ipx][iw].i;
		}
		for (iw = nw; iw < ntfft; iw++)
		{
			cb[iw].r = 0.0;
			cb[iw].i = 0.0;
		}

		/* Fourier transform frequency to time */
		pfacr(-1, ntfft, cf[ipx], gpad);
		pfacc(-1, ntfft, cb);

		/* copy to output array */
		for (itau = 0; itau < ntau; ++itau)
			f[ipx][itau] = (gpad[itau] + cb[itau].i);
		//			printf("%f %f\n",gpad[itau]*10000000,imgpad[itau]*10000000);
	}

	free1float(imfpad);
	free1float(fpad);
	free1float(gpad);
	free1float(imgpad);
	free2complex(cf);
	free2complex(imcf);
	free1complex(cb);
}

void synData(float bwh, float dxb, float fmin,
			 int nt, float dt, float ft,
			 int nx, float dx, float fx, float **f,
			 int ntau, float dtau, float ftau,
			 int npx, float dpx, float fpx, float **g)
/*****************************************************************************
Form beams (filtered slant stacks) for later superposition of Gaussian beams.
*****************************************************************************/
{
	int ntpad, ntfft, nw, ix, iw, ipx, it, itau;
	float wmin, pxmax, xmax, x, dw, fw, w, fftscl,
		amp, phase, scale, a, b, as, bs, es, cfr, cfi,
		*fpad, *gpad;
	complex **cf, **cg, *cfx, *cgpx;
	float ffw, tmpp, px;

	/* minimum frequency in radians */
	wmin = 2.0 * PI * fmin;

	/* pad time axis to avoid wraparound */
	pxmax = (dpx < 0.0) ? fpx : fpx + (npx - 1) * dpx;
	xmax = (dx < 0.0) ? fx : fx + (nx - 1) * dx;
	ntpad = ABS(pxmax * xmax) / dt;
	ntpad = 0;

	/* fft sampling */
	ntfft = npfar(MAX(nt + ntpad, ntau));
	nw = ntfft / 2 + 1;
	dw = 2.0 * PI / (ntfft * dt);
	fw = 0.0;
	fftscl = 1.0 / ntfft;

	/* allocate workspace */
	fpad = alloc1float(ntfft);
	gpad = alloc1float(ntfft);
	cf = alloc2complex(nw, npx);
	cg = alloc2complex(nw, nx);
	cfx = alloc1complex(npx);
	cgpx = alloc1complex(nx);

	/* loop over x */
	for (ipx = 0; ipx < npx; ++ipx)
	{

		/* pad time with zeros */
		for (it = 0; it < nt; ++it)
			fpad[it] = g[ipx][it];
		for (it = nt; it < ntfft; ++it)
			fpad[it] = 0.0;

		/* Fourier transform time to frequency */
		pfarc(1, ntfft, fpad, cf[ipx]);
	}

	/* loop over w */
	for (iw = 0, w = fw; iw < nw; ++iw, w += dw)
	{

		for (ipx = 0; ipx < npx; ipx++)
		{
			cfx[ipx].r = cf[ipx][iw].r;
			cfx[ipx].i = cf[ipx][iw].i;
		}

		/*transform p to x*/
		xtop(-w, npx, dpx, fpx, cfx, nx, dx, fx, cgpx);

		/* frequency-dependent amplitude scale factors */
		scale = -0.5 * w / (wmin * bwh * bwh);
		amp = -dx * fftscl * dxb * sqrt(w) / bwh / sqrt(2 * PI * wmin) / (4 * PI * PI); // minus sign is entered here
		a = b = amp;

		for (ix = 0, x = fx; ix < nx; ++ix, x += dx)
		{
			es = exp(scale * x * x);
			as = a * es;
			bs = b * es;
			cgpx[ix].r = cgpx[ix].r * as;
			cgpx[ix].i = cgpx[ix].i * bs;
			//			cfx[ix].r = sqrt(2*w)*((as*cfr-bs*cfi)-(bs*cfr+as*cfi))/2;
			//			cfx[ix].i = sqrt(2*w)*((bs*cfr+as*cfi)+(as*cfr-bs*cfi))/2;
		}

		for (ix = 0; ix < nx; ++ix)
		{
			cg[ix][iw].r = cgpx[ix].r;
			cg[ix][iw].i = cgpx[ix].i;
		}
	}
	/* loop over px */
	for (ix = 0; ix < nx; ++ix)
	{

		/* Fourier transform frequency to time */
		pfacr(-1, ntfft, cg[ix], gpad);

		/* copy to output array */
		for (itau = 0; itau < ntau; ++itau)
			f[ix][itau] += gpad[itau];
	}

	free1float(fpad);
	free1float(gpad);
	free2complex(cf);
	free2complex(cg);
	free1complex(cfx);
	free1complex(cgpx);
}

/* functions for external use */

void formBeams(float bwh, float dxb, float fmin,
			   int nt, float dt, float ft,
			   int nx, float dx, float fx, float **f,
			   int ntau, float dtau, float ftau,
			   int npx, float dpx, float fpx, float **g)
/*****************************************************************************
Form beams (filtered slant stacks) for later superposition of Gaussian beams.
******************************************************************************
	  nput:
bwh		horizontal beam half-width
dxb		horizontal distance between beam centers
fmin		minimum frequency (cycles per unit time)
nt		number of input time samples
dt		input time sampling interval
ft		first input time sample
nx		number of horizontal samples
dx		horizontal sampling interval
fx		first horizontal sample
f		array[nx][nt] of data to be slant stacked into beams
ntau		number of output time samples
dtau		output time sampling interval (currently must equal dt)
ftau		first output time sample
npx		number of horizontal slownesses
dpx		horizontal slowness sampling interval
fpx		first horizontal slowness

Output:
g		array[npx][ntau] containing beams
*****************************************************************************/
{
	int ntpad, ntfft, nw, ix, iw, ipx, it, itau;
	float wmin, pxmax, xmax, x, dw, fw, w, fftscl,
		amp, phase, scale, a, b, as, bs, es, cfr, cfi,
		*fpad, *gpad;
	complex **cf, **cg, *cfx, *cgpx;
	int nf1, nf2, nf3, nf4;
	float ffw, tmpp;
	float spr;
	float wd = 2.0 * PI * DOMINFY;

	/* minimum frequency in radians */
	wmin = 2.0 * PI * fmin;

	/* pad time axis to avoid wraparound */
	pxmax = (dpx < 0.0) ? fpx : fpx + (npx - 1) * dpx;
	xmax = (dx < 0.0) ? fx : fx + (nx - 1) * dx;
	ntpad = ABS(pxmax * xmax) / dt;
	ntpad = 0;

	/* fft sampling */
	ntfft = npfar(MAX(nt + ntpad, ntau));
	nw = ntfft / 2 + 1;
	dw = 2.0 * PI / (ntfft * dt);
	fw = 0.0;
	fftscl = 1.0 / ntfft;

	/* allocate workspace */
	fpad = alloc1float(ntfft);
	gpad = alloc1float(ntfft);
	cf = alloc2complex(nw, nx);
	cg = alloc2complex(nw, npx);
	cfx = alloc1complex(nx);
	cgpx = alloc1complex(npx);

	/* loop over x */
	for (ix = 0; ix < nx; ++ix)
	{

		/* pad time with zeros */
		for (it = 0; it < nt; ++it)
			fpad[it] = f[ix][it]; //, printf("%d %f\n",DOMINFY,fpad[it]*10000000);
		for (it = nt; it < ntfft; ++it)
			fpad[it] = 0.0;

		/* Fourier transform time to frequency */
		pfarc(1, ntfft, fpad, cf[ix]);
	}

	/* loop over w */
	for (iw = 0, w = fw; iw < nw; ++iw, w += dw)
	{

		/* frequency-dependent amplitude scale factors */
		spr = 4 * w * w * sqrt(PI) * exp(-w * w / wd / wd) / wd / wd / wd;
		// spr = 1.0;
		scale = -0.5 * w / (wmin * bwh * bwh);
		amp = -dpx * w * w * spr * fftscl * dxb * sqrt(w) / bwh / sqrt(2 * PI * wmin) / (4 * PI * PI); // minus sign is entered here and dpx

		/* phase shift to account for ft not equal to ftau */
		//		phase = w*(ft-ftau);

		/* apply complex filter */
		a = b = amp;
		for (ix = 0, x = fx; ix < nx; ++ix, x += dx)
		{
			es = exp(scale * x * x);
			as = a * es;
			bs = b * es;
			cfr = cf[ix][iw].r;
			cfi = cf[ix][iw].i;
			//			cfx[ix].r = sqrt(2*w)*((as*cfr-bs*cfi)-(bs*cfr+as*cfi))/2;
			//			cfx[ix].i = sqrt(2*w)*((bs*cfr+as*cfi)+(as*cfr-bs*cfi))/2;
			cfx[ix].r = cf[ix][iw].r * as;
			cfx[ix].i = cf[ix][iw].i * bs;
		}

		/* transform x to p */
		xtop(w, nx, dx, fx, cfx, npx, dpx, fpx, cgpx);
		for (ipx = 0; ipx < npx; ++ipx)
		{
			cg[ipx][iw].r = cgpx[ipx].r;
			cg[ipx][iw].i = cgpx[ipx].i;
		}
	}
	/* loop over px */
	for (ipx = 0; ipx < npx; ++ipx)
	{

		/* Fourier transform frequency to time */
		pfacr(-1, ntfft, cg[ipx], gpad);

		/* copy to output array */
		for (itau = 0; itau < ntau; ++itau)
			g[ipx][itau] = gpad[itau];
	}

	free1float(fpad);
	free1float(gpad);
	free2complex(cf);
	free2complex(cg);
	free1complex(cfx);
	free1complex(cgpx);
}

void accray(Ray *ray, Cell **c, float fmin, float lmin, int lx, int lz, int mx,
			int mz, int live, int ipx, int nt, float dt, float ft,
			int nx, float dx, float fx, int nz, float dz, float fz, float **g, float **v)
/*****************************************************************************
Accumulate contribution of one Gaussian beam.
******************************************************************************
	  nput:
ray		ray parameters sampled at discrete ray steps
fmin		minimum frequency (cycles per unit time)
lmin		initial beam width for frequency wmin
nt		number of time samples
dt		time sampling interval
ft		first time sample
f		array[nt] containing data for one ray f(t)
nx		number of x samples
dx		x sampling interval
fx		first x sample
nz		number of z samples
dz		z sampling interval
fz		first z sample
g		array[nx][nz] in which to accumulate beam

Output:
g		array[nx][nz] after accumulating beam
*****************************************************************************/
{
	int jx, jz;
	int i;
	float wmin;
	RayStep *rs = ray->rs;
	Cells *cells;
	cells = (Cells *)alloc1(1, sizeof(Cells));
	/* frequency in radians per unit time */
	wmin = 2.0 * PI * fmin;

	/* set information needed to set and fill cells */
	cells->nt = nt;
	cells->dt = dt;
	cells->ft = ft;
	cells->lx = lx;
	cells->mx = mx;
	cells->nx = nx;
	cells->dx = dx;
	cells->fx = fx;
	cells->lz = lz;
	cells->mz = mz;
	cells->nz = nz;
	cells->dz = dz;
	cells->fz = fz;
	cells->live = live;
	cells->ip = ipx;
	cells->wmin = wmin;
	cells->lmin = lmin;
	cells->cell = c;
	cells->ray = ray;
	//	cells->bd = bd;
	//	cells->g = g;

	/* cell closest to initial point on ray will be first live cell */
	jx = NINT((rs[0].x - fx) / dx / lx);
	jz = NINT((rs[0].z - fz) / dz / lz);

	/* set first live cell and its neighbors recursively */
	setCell(cells, jx, jz);
	/*	for(jx=0;jx<mx-1;jx++)
			for(jz=0;jz<mz-1;jz++)
				if(c[jx][jz].live==live)
					printf("%f\n",c[jx][jz].angle);*/

	/* free complex beam data */
	//	free2complex(bd->cf);
	//	free1((void*)bd);

	/* free cells */
	//	free2((void**)cells->cell);
	free1((void *)cells);
}

/* functions for internal use only */

static void xtop(float w,
				 int nx, float dx, float fx, complex *g,
				 int np, float dp, float fp, complex *h)
/*****************************************************************************
Slant stack for one frequency w, where slant stack is defined by

		   fx+(nx-1)*dx
	h(p) =   integral   exp(-sqrt(-1)*w*p*x) * g(x) * dx
				fx
******************************************************************************
	  nput:
w		frequency (radians per unit time)
nx		number of x samples
dx		x sampling interval
fx		first x sample
g		array[nx] containing g(x)
np		number of p samples
dp		p sampling interval
fp		first p sample

Output:
h		array[np] containing h(p)
******************************************************************************
Notes:
The units of w, x, and p must be consistent.

Slant stack is performed via FFT and 8-point (tapered-sinc) interpolation.

The Fourier transform over time (t) is assumed to have been performed with
positive sqrt(-1)*w*t in the exponent;  if negative sqrt(-1)*w*t was used
instead, call this function with negative w.
*****************************************************************************/
{
	int nxfft, nk, nka, ix, ik, ip, lwrap;
	float dk, fk, ek, fka, k, p, phase, c, s, x, xshift, temp, *kp;
	complex czero, *gx, *gk, *gka, *hp;
	czero = cmplx(0.0, 0.0);

	/* number of samples required to make wavenumber k periodic */
	lwrap = 8;

	/* wavenumber k sampling */
	nxfft = npfa((nx + lwrap) * 2);
	nk = nxfft;
	dk = 2.0 * PI / (nxfft * dx);
	fk = -PI / dx;
	ek = PI / dx;
	fka = fk - lwrap * dk;
	nka = lwrap + nk + lwrap;

	/* allocate workspace */
	gka = alloc1complex(nka);
	gx = gk = gka + lwrap;
	hp = alloc1complex(np);
	kp = alloc1float(np);

	/* scale g(x) by x sampling interval dx */
	for (ix = 0; ix < nx; ++ix, x += dx)
	{
		gx[ix].r = dx * g[ix].r;
		gx[ix].i = dx * g[ix].i;
	}
	/* pad g(x) with zeros */
	for (ix = nx; ix < nxfft; ++ix)
		gx[ix].r = gx[ix].i = 0.0;

	/* negate every other sample so k-axis will be centered */
	for (ix = 1; ix < nx; ix += 2)
	{
		gx[ix].r = -gx[ix].r;
		gx[ix].i = -gx[ix].i;
	}

	/* Fourier transform g(x) to g(k) */
	pfacc(-1, nxfft, gx);

	/* wrap-around g(k) to avoid interpolation end effects */
	for (ik = 0; ik < lwrap; ++ik)
		gka[ik] = gk[ik + nk - lwrap];
	for (ik = lwrap + nk; ik < lwrap + nk + lwrap; ++ik)
		gka[ik] = gk[ik - lwrap - nk];

	/* phase shift to account for non-centered x-axis */
	xshift = 0.5 * (nx - 1) * dx;
	for (ik = 0, k = fka; ik < nka; ++ik, k += dk)
	{
		phase = k * xshift;
		c = cos(phase);
		s = sin(phase);
		temp = gka[ik].r * c - gka[ik].i * s;
		gka[ik].i = gka[ik].r * s + gka[ik].i * c;
		gka[ik].r = temp;
	}

	/* compute k values at which to interpolate g(k) */
	for (ip = 0, p = fp; ip < np; ++ip, p += dp)
	{
		kp[ip] = w * p;

		/* if outside Nyquist bounds, do not interpolate */
		if (kp[ip] < fk && kp[ip] < ek)
			kp[ip] = fk - 1000.0 * dk;
		else if (kp[ip] > fk && kp[ip] > ek)
			kp[ip] = ek + 1000.0 * dk;
	}

	/* interpolate g(k) to obtain h(p) */
	ints8c(nka, dk, fka, gka, czero, czero, np, kp, hp);

	/* phase shift to account for non-centered x-axis and non-zero fx */
	xshift = -fx - 0.5 * (nx - 1) * dx;
	for (ip = 0; ip < np; ++ip)
	{
		phase = kp[ip] * xshift;
		c = cos(phase);
		s = sin(phase);
		h[ip].r = hp[ip].r * c - hp[ip].i * s;
		h[ip].i = hp[ip].r * s + hp[ip].i * c;
	}

	/* free workspace */
	free1complex(gka);
	free1complex(hp);
	free1float(kp);
}

static BeamData *beamData(int type, int npx, float wmin, int nt, float dt, float ft, float **f)
/*****************************************************************************
Compute filtered complex beam data as a function of real and imaginary time.
******************************************************************************
	  nput:
wmin		minimum frequency (in radians per unit time)
nt		number of time samples
dt		time sampling interval
ft		first time sample
f		array[nt] containing data to be filtered

Returned:	pointer to beam data
*****************************************************************************/
{
	int ntpad, ntfft, nw, iwnyq, ntrfft, ntr, nti, nwr, it, itr, iti, iw, ipx;
	float dw, fw, dtr, ftr, dti, fti, w, ti, scale, *fa;
	complex *ca, *cb, *cfi, ***cf;
	BeamData *bd;
	int nf1, nf2, nf3, nf4;
	float ffw, tmpp;

	/* pad to avoid wraparound in Hilbert transform */
	ntpad = 0;

	/* fft sampling */
	ntfft = npfaro(nt + ntpad, 2 * (nt + ntpad));
	nw = ntfft / 2 + 1;
	dw = 2.0 * PI / (ntfft * dt);
	fw = 0.0;
	iwnyq = nw - 1;

	/* real time sampling (oversample for future linear interpolation) */
	ntrfft = nwr = npfao(NOVERSAMPLE * ntfft, NOVERSAMPLE * ntfft + ntfft);
	dtr = dt * ntfft / ntrfft;
	ftr = ft;
	ntr = (1 + (nt + ntpad - 1) * dt / dtr);
	//	printf("ntr=%d,dtr=%d,ftr=%f\n",ntr,dtr,ftr+(ntr-1)*dtr);

	/* imaginary time sampling (exponential decay filters) */
	nti = NFILTER;
	dti = EXPMIN / (wmin * (nti - 1));
	fti = 0.0;

	if (type == -1)
	{

		/* allocate space for filtered data */
		cf = alloc3complex(ntr, nti, npx);
		for (ipx = 0; ipx < npx; ipx++)
		{
			/* allocate workspace */
			fa = alloc1float(ntfft);
			ca = alloc1complex(nw);
			cb = alloc1complex(ntrfft);

			/* pad data with zeros */
			for (it = 0; it < nt; ++it)
				fa[it] = f[ipx][it];
			for (it = nt; it < ntfft; ++it)
				fa[it] = 0.0;

			/* Fourier transform and scale to make complex analytic signal */
			pfarc(1, ntfft, fa, ca);
			for (iw = 1; iw < iwnyq; ++iw)
			{
				ca[iw].r *= 2.0;
				ca[iw].i *= 2.0;
			}

			/* loop over imaginary time */
			for (iti = 0, ti = fti; iti < nti; ++iti, ti += dti)
			{

				/* apply exponential decay filter */
				for (iw = 0, w = fw; iw < nw; ++iw, w += dw)
				{
					// scale = exp(w*ti);
					scale = 1.0;
					cb[iw].r = (ca[iw].r) * scale;
					cb[iw].i = (ca[iw].i) * scale;
				}
				/* pad with zeros */
				for (iw = nw; iw < nwr; ++iw)
					cb[iw].r = cb[iw].i = 0.0;

				/* inverse Fourier transform and scale */
				pfacc(-1, ntrfft, cb);

				cfi = cf[ipx][iti];
				scale = 1.0 / ntfft;

				for (itr = 0; itr < ntr; ++itr)
				{
					cfi[itr].r = scale * cb[itr].r;
					cfi[itr].i = scale * cb[itr].i;
				}
			}

			/* free workspace */
			free1float(fa);
			free1complex(ca);
			free1complex(cb);
		}
	}

	/* return beam data */
	bd = (BeamData *)alloc1(1, sizeof(BeamData));

	if (type == -1)
	{
		bd->cf = cf;
	}

	bd->ntr = ntr;
	bd->dtr = dtr;
	bd->ftr = ftr;
	bd->nti = nti;
	bd->dti = dti;
	bd->fti = fti;
	//	bd->cf = cf;
	return bd;
}

static void setCell(Cells *cells, int jx, int jz)
/*****************************************************************************
Set a cell by computing its Gaussian beam complex time and amplitude.
	  f the amplitude is non-zero, set neighboring cells recursively.
******************************************************************************
	  nput:
cells		pointer to cells
jx		x index of the cell to set
jz		z index of the cell to set
******************************************************************************
Notes:
To reduce the amount of memory required for recursion, the actual
computation of complex time and amplitude is performed by the cellTimeAmp()
function, so that no local variables are required in this function, except
for the input arguments themselves.
*****************************************************************************/
{
	int ix, iz;
#if 1
	/* if cell is out of bounds, return */
	if (jx < 0 || jx >= cells->mx || jz < 0 || jz >= cells->mz)
		return;

	/* if cell is live, return */
	if (cells->cell[jx][jz].live == cells->live)
		return;

	/* make cell live */
	cells->cell[jx][jz].live = cells->live;
	cells->cell[jx][jz].ip = cells->ip;
	/* compute complex time and amplitude.  If amplitude is
	 * big enough, recursively set neighboring cells. */
	if (cellTimeAmp(cells, jx, jz))
	{
		setCell(cells, jx + 1, jz);
		setCell(cells, jx - 1, jz);
		setCell(cells, jx, jz + 1);
		setCell(cells, jx, jz - 1);
	}
#endif

#if 0
	for(ix=0;ix<mx;ix++) {
		for(iz=0;iz<mz;iz++) {
			if(cellTimeAmp(cells,ix,iz)) {
			cells->cell[ix][iz].live = cells->live;
			cells->cell[ix][iz].ip=cells->ip;
			}
		}
	}
#endif
}

static int cellTimeAmp(Cells *cells, int jx, int jz)
/*****************************************************************************
Compute complex and time and amplitude for a cell.
******************************************************************************
	  nput:
cells		pointer to cells
jx		x index of the cell to set
jz		z index of the cell to set

Returned:	1 if Gaussian amplitude is significant, 0 otherwise
*****************************************************************************/
{
	int lx = cells->lx, lz = cells->lz;
	float dx = cells->dx, fx = cells->fx, dz = cells->dz, fz = cells->fz,
		  wmin = cells->wmin, lmin = cells->lmin;
	Ray *ray = cells->ray;
	Cell **cell = cells->cell;
	int irs, kmah;
	float tmax, xc, zc, t0, x0, z0, x, z, tr, ti, ar, ai,
		v, q1, p1, q2, p2, e, es, scale, phase, vvmr, vvmi, c, s, ov,
		px, pz, pzpz, pxpx, pxpz, dvdx, dvdz, dvds, dvdn,
		wxxr, wxzr, wzzr, wxxi, wxzi, wzzi;
	float ex, ez, ang, halfbw;
	float lastx, lastz, lastc, lasts, projx, projz;
	RayStep *rs;

	/* maximum time */
	tmax = ray->rs[ray->nrs - 1].t;
	lastx = ray->rs[ray->nrs - 1].x;
	lastz = ray->rs[ray->nrs - 1].z;
	lasts = ray->rs[ray->nrs - 1].s;
	lastc = ray->rs[ray->nrs - 1].c;

	/* cell coordinates */
	xc = fx + jx * lx * dx;
	zc = fz + jz * lz * dz;

	/* vectors from last ray step to cell */
	projx = xc - lastx;
	projz = zc - lastz;

	/* ray step nearest to cell */
	irs = nearestRayStep(ray, xc, zc);
	rs = &(ray->rs[irs]);

	/* ray time and coordinates */
	t0 = rs->t;
	x0 = rs->x;
	z0 = rs->z;

	/* real and imaginary parts of v*v*m = v*v*p/q */
	v = rs->v;
	q1 = rs->q1;
	p1 = rs->p1;
	q2 = rs->q2;
	p2 = rs->p2;
	e = wmin * lmin * lmin;
	es = e * e;
	scale = v * v / (q2 * q2 + es * q1 * q1);
	vvmr = (p2 * q2 + es * p1 * q1) * scale;
	vvmi = -e * scale;
	halfbw = sqrt(2 * (es * q1 * q1 + q2 * q2) / fabs(wmin * e * (p2 * q1 - p1 * q2)));

	/* components of slowness vector, px and pz */
	c = rs->c;
	s = rs->s;
	ov = 1.0 / v;
	px = s * ov;
	pz = c * ov;
	pzpz = pz * pz;
	pxpx = px * px;
	pxpz = px * pz;

	/* velocity derivatives along tangent and normal */
	dvdx = rs->dvdx;
	dvdz = rs->dvdz;
	dvds = s * dvdx + c * dvdz;
	dvdn = c * dvdx - s * dvdz;

	/* real part of W matrix */
	wxxr = pzpz * vvmr - 2.0 * pxpz * dvdn - pxpx * dvds;
	wxzr = (pxpx - pzpz) * dvdn - pxpz * (vvmr + dvds);
	wzzr = pxpx * vvmr + 2.0 * pxpz * dvdn - pzpz * dvds;

	/* imaginary part of W matrix */
	wxxi = pzpz * vvmi;
	wxzi = -pxpz * vvmi;
	wzzi = pxpx * vvmi;

	/* vector from ray to cell */
	x = xc - x0;
	z = zc - z0;

	/* real and imaginary parts of complex time */
	tr = t0 + (px + 0.5 * wxxr * x) * x + (pz + wxzr * x + 0.5 * wzzr * z) * z;
	ti = (0.5 * wxxi * x * x + (wxzi * x + 0.5 * wzzi * z) * z) * 4.0;

	/* real and imaginary parts of complex amplitude */
	kmah = rs->kmah;
	scale = pow(v * v / (q2 * q2 + es * q1 * q1), 0.25);
	phase = 0.5 * (atan2(q2, q1 * e) + 2.0 * PI * ((kmah + 1) / 2));
	ar = scale * cos(phase);
	ai = scale * sin(phase);

	/* determine ray directions of the gird */
	ex = s + (x * pz - z * px) * vvmr * c;
	ez = c - (x * pz - z * px) * vvmr * s;

	if ((ex > 0) && (ez > 0))
		ang = atan(ex / ez);
	if ((ex > 0) && (ez < 0))
		ang = atan(ex / ez) + PI;
	if ((ex < 0) && (ez > 0))
		ang = atan(ex / ez);
	if ((ex < 0) && (ez < 0))
		ang = -PI + atan(px / pz);

	/* set cell parameters */
	cell[jx][jz].tr = tr;
	cell[jx][jz].ti = ti;
	cell[jx][jz].ar = ar;
	cell[jx][jz].ai = ai;
	cell[jx][jz].angle = ang;

	/* return 1 if Gaussian amplitude is significant, 0 otherwise */
	return (wmin * ti > EXPMIN && tr <= tmax && sqrt(x * x + z * z) < halfbw * 0.75 && (projx * lasts + projz * lastc) < 0.0) ? 1 : 0;
	// return (wmin*ti>EXPMIN && tr<=tmax && (projx*lasts+projz*lastc)<0.0)?1:0;
	// return (wmin*ti>EXPMIN && tr<=tmax)?1:0;
}

void scanimg(int type, Cell ***c1, Cell ***c2, int live, int dead, int nt, float dt, float ft,
			 float fmin, int lx, int lz, int nx, int nz, int mx, int mz, int npx, float **f,
			 float **g, float fpx, float dpx, float **imb)
/************************************************************************************
Scan over offset ray parameters then image with beams.
*************************************************************************************
	  nput parameters:
************************************************************************************/
{
	int pm, ips, ph, ipr, imx, imz, tp;
	float wmin, tmin, tempt;
	Cell **fncell;
	BeamData *bd;
	wmin = 2.0 * PI * fmin;
	pm = 2 * npx - 1;

	bd = beamData(type, npx, wmin, nt, dt, ft, f);

	/* for each midpoint ray parameter*/
	for (ips = 0; ips < npx; ips++)
	{

		for (ipr = 0; ipr < npx; ipr++)
		{

			fncell = (Cell **)alloc2(mz, mx, sizeof(Cell));
			memset(&fncell[0][0], 0, sizeof(Cell) * mz * mx);

			for (imx = 0; imx < mx; imx++)
			{
				for (imz = 0; imz < mz; imz++)
				{
					// fncell[imx][imz].ip=0;
					// fncell[imx][imz].live=0;
					tmin = 1000.0;
					tp = 1000;

					if ((c1[ips][imx][imz].live == live) && (c2[ipr][imx][imz].live == live))
					{
						tempt = c1[ips][imx][imz].ti + c2[ipr][imx][imz].ti;
						//	  printf("%d %d %f\n",ipr, tp, tempt );

						if ((fabs(tempt) < 5.0))
							tmin = fabs(tempt),
							tp = ips;
					}

					if ((tp < 1000))
					{

						fncell[imx][imz].ai = c1[ips][imx][imz].ar * c2[ipr][imx][imz].ai + c1[ips][imx][imz].ai * c2[ipr][imx][imz].ar;
						fncell[imx][imz].ar = c1[ips][imx][imz].ar * c2[ipr][imx][imz].ar - c1[ips][imx][imz].ai * c2[ipr][imx][imz].ai;
						fncell[imx][imz].ti = -tmin;
						fncell[imx][imz].ai *= exp(-2 * wmin * tmin);
						fncell[imx][imz].ar *= exp(-2 * wmin * tmin);
						fncell[imx][imz].tr = c1[ips][imx][imz].tr + c2[ipr][imx][imz].tr;
						fncell[imx][imz].live = live;
						fncell[imx][imz].angle = c1[ips][imx][imz].angle - c2[ipr][imx][imz].angle;
						fncell[imx][imz].ip = c2[ipr][imx][imz].ip;

					} // end if tp

				} // end loop imz
			}	  // end loop imx

			for (imx = 0; imx < mx - 1; imx++)
			{
				for (imz = 0; imz < mz - 1; imz++)
				{
					if ((fncell[imx][imz].live == live) && (fncell[imx + 1][imz].live == live) &&
						(fncell[imx][imz + 1].live == live) && (fncell[imx + 1][imz + 1].live == live) && (abs(fncell[imx + 1][imz + 1].ip - fncell[imx][imz + 1].ip) < 2) && (abs(fncell[imx][imz + 1].ip - fncell[imx][imz].ip) < 2) && (abs(fncell[imx][imz].ip - fncell[imx + 1][imz].ip) < 2) && (abs(fncell[imx + 1][imz].ip - fncell[imx + 1][imz + 1].ip) < 2))
					{

						// printf("%d %d %d %f\n",ipr,fncell[imx][imz].ip,c2[ipr][imx][imz].ip,fncell[imx][imz].tr);
						cellBeam(type, fncell, f, g, imx, imz, bd, lx, lz, nx, nz, npx, ips, ipr, fpx, dpx, ft, dt, nt, imb);
					} // end if

				} // end loop imz
			}	  // end loop imx

			free2((void **)fncell);

		} // endif ipr
	}	  // endif ips

	/* if MIG, free complex filtered data */
	if (type == -1)
	{
		free3complex(bd->cf);
	}
}

static void cellBeam(int type, Cell **cell, float **f, float **g, int jx, int jz, BeamData *bd, int lx, int lz,
					 int nx, int nz, int npx, int ips, int ipr, float fpx, float dpx,
					 float ft, float dt, int nt, float **imb)
/*****************************************************************************
Accumulate Gaussian beam for one cell.
******************************************************************************
	  nput:
cells		pointer to cells
jx		x index of the cell in which to accumulate beam
jz		z index of the cell in which to accumulate beam
*****************************************************************************/
{
	int ntr = bd->ntr, nti = bd->nti;
	float dtr = bd->dtr, ftr = bd->ftr, dti = bd->dti, fti = bd->fti;
	complex ***cf = bd->cf;
	int kxlo, kxhi, kzlo, kzhi, kx, kz, itr, iti, np;
	float ta00r, ta01r, ta10r, ta11r, ta00i, ta01i, ta10i, ta11i,
		aa00r, aa01r, aa10r, aa11r, aa00i, aa01i, aa10i, aa11i,
		tax0r, tax1r, tax0i, tax1i, aax0r, aax1r, aax0i, aax1i,
		taxzr, taxzi, aaxzr, aaxzi,
		dtax0r, dtax0i, dtax1r, dtax1i, daax0r, daax0i, daax1r, daax1i,
		dtaxzr, dtaxzi, daaxzr, daaxzi;

	float odtr, odti, xdelta, zdelta, trn, tin, trfrac, mtrfrac, tifrac, mtifrac,
		cf0r, cf0i, cf1r, cf1i, cfr, cfi;
	float a00g, a01g, a10g, a11g;

	float da0g, da1g, a0g, a1g, azg, dazg;

	complex *cf0, *cf1;
	int iangle;
	float ktr, dktr;

	/* inverse of time sampling intervals */
	//	ftr=ftr-0.9;
	odtr = 1.0 / dtr;
	odti = 1.0 / dti;

	/* complex time and amplitude for each corner */
	ta00r = cell[jx][jz].tr;
	ta01r = cell[jx][jz + 1].tr;
	ta10r = cell[jx + 1][jz].tr;
	ta11r = cell[jx + 1][jz + 1].tr;
	ta00i = cell[jx][jz].ti;
	ta01i = cell[jx][jz + 1].ti;
	ta10i = cell[jx + 1][jz].ti;
	ta11i = cell[jx + 1][jz + 1].ti;
	aa00r = cell[jx][jz].ar;
	aa01r = cell[jx][jz + 1].ar;
	aa10r = cell[jx + 1][jz].ar;
	aa11r = cell[jx + 1][jz + 1].ar;
	aa00i = cell[jx][jz].ai;
	aa01i = cell[jx][jz + 1].ai;
	aa10i = cell[jx + 1][jz].ai;
	aa11i = cell[jx + 1][jz + 1].ai;

	/* opening angle for each corner */
	a00g = cell[jx][jz].angle;
	a01g = cell[jx][jz + 1].angle;
	a10g = cell[jx + 1][jz].angle;
	a11g = cell[jx + 1][jz + 1].angle;

	/* x and z samples for cell */
	kxlo = jx * lx;
	kxhi = kxlo + lx;
	if (kxhi > nx)
		kxhi = nx;
	kzlo = jz * lz;
	kzhi = kzlo + lz;
	if (kzhi > nz)
		kzhi = nz;

	/* fractional increments for linear interpolation */
	xdelta = 1.0 / lx;
	zdelta = 1.0 / lz;

	/* increments for times and amplitudes at top and bottom of cell */
	dtax0r = (ta10r - ta00r) * xdelta;
	dtax1r = (ta11r - ta01r) * xdelta;
	dtax0i = (ta10i - ta00i) * xdelta;
	dtax1i = (ta11i - ta01i) * xdelta;
	daax0r = (aa10r - aa00r) * xdelta;
	daax1r = (aa11r - aa01r) * xdelta;
	daax0i = (aa10i - aa00i) * xdelta;
	daax1i = (aa11i - aa01i) * xdelta;
	da0g = (a10g - a00g) * xdelta;
	da1g = (a11g - a01g) * xdelta;

	/* times and amplitudes at top-left and bottom-left of cell */
	tax0r = ta00r;
	tax1r = ta01r;
	tax0i = ta00i;
	tax1i = ta01i;
	aax0r = aa00r;
	aax1r = aa01r;
	aax0i = aa00i;
	aax1i = aa01i;
	a0g = a00g;
	a1g = a01g;

	/* loop over x samples */
	for (kx = kxlo; kx < kxhi; ++kx)
	{

		/* increments for time and amplitude */
		dtaxzr = (tax1r - tax0r) * zdelta;
		dtaxzi = (tax1i - tax0i) * zdelta;
		daaxzr = (aax1r - aax0r) * zdelta;
		daaxzi = (aax1i - aax0i) * zdelta;
		dazg = (a1g - a0g) * zdelta;

		/* time and amplitude at top of cell */
		taxzr = tax0r;
		taxzi = tax0i;
		aaxzr = aax0r;
		aaxzi = aax0i;
		azg = a0g;

		/* loop over z samples */
		for (kz = kzlo; kz < kzhi; ++kz)
		{

			/* index of imaginary time */
			iti = tin = (taxzi - fti) * odti;
			//			printf("%f %f %f %d %d %f\n",taxzr,taxzi,dtaxzr,cell[jx][jz].ip,ipr,tin);

			if (fabs(azg) > 120 * PI / 180)
				continue;
			if (iti < 0 || iti >= nti - 1)
				continue;

			np = ipr;

			if (type == 1)
			{
				itr = trn = (taxzr - ft) / dt;
				dktr = (taxzr - itr * dt) / dt;
				if (itr < 0 || itr > nt - 1)
					continue;

#if 0
				/* real part contribution */			
				f[ipr][itr] += (1-dktr)*aaxzr*g[kx][kz];
				f[ipr][itr+1] += dktr*aaxzr*g[kx][kz];

				/* imaginary part contribution */
				imb[ipr][itr] += (1-dktr)*aaxzi*g[kx][kz];
				imb[ipr][itr+1] += dktr*aaxzi*g[kx][kz];
#endif

#if 1
				f[ipr][itr] += aaxzr * g[kx][kz];
				imb[ipr][itr] += aaxzi * g[kx][kz];
#endif
			}

			if (type == -1)
			{
				/* pointers to left and right imaginary time samples */
				cf0 = cf[np][iti];
				cf1 = cf[np][iti + 1];

				/* imaginary time linear interpolation coefficients */
				tifrac = tin - iti;
				mtifrac = 1.0 - tifrac;

				/* index of real time */
				itr = trn = (taxzr - ftr) * odtr;
				if (itr < 0 || itr >= ntr - 1)
					continue;

				/* real time linear interpolation coefficients */
				trfrac = trn - itr;
				mtrfrac = 1.0 - trfrac;

				/* real and imaginary parts of complex beam data */
				cfr = mtrfrac * cf0[itr].r + trfrac * cf0[itr + 1].r;
				cfi = mtrfrac * cf0[itr].i + trfrac * cf0[itr + 1].i;

				//				cfr = cf0[itr].r;
				//				cfi = cf0[itr].i;

				/* accumulate beam */
				g[kx][kz] += (aaxzr * cfr - aaxzi * cfi);
				//			printf("%f %f %f %f\n",aaxzr,aaxzi,cfr,cfi);
			}

			/* increment time and amplitude */
			taxzr += dtaxzr;
			taxzi += dtaxzi;
			aaxzr += daaxzr;
			aaxzi += daaxzi;
			azg += dazg;
		}

		/* increment times and amplitudes at top and bottom of cell */
		tax0r += dtax0r;
		tax1r += dtax1r;
		tax0i += dtax0i;
		tax1i += dtax1i;
		aax0r += daax0r;
		aax1r += daax1r;
		aax0i += daax0i;
		aax1i += daax1i;
		a0g += da0g;
		a1g += da1g;
	}
}
