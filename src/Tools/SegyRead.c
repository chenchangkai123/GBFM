#include "SegyRead.h"


/* 
SegyRead
Chang Zhimiao 2023.06.07
  */

// int main()
// {
//     int i,j,k;
//     int segyFileNum = 1;
//     char **index;

//     FILE *fp1;

//     fp1 = fopen("Index","r");
//     index = IndexProcess(fp1,segyFileNum);
    
    
//     int *header;
//     float *data;

//     FILE *fp;
//     FILE *headerOut;
//     FILE *dataOut;

//     fp= fopen(index[0],"rb");
//     headerOut = fopen("header.txt","w");
//     dataOut = fopen("data.txt","w");

//     // 跳过3600字节的二进制文件头
//     fseek(fp,3600,SEEK_SET);

//     // 读取240字节的二进制文件头
//     header = (int *)malloc(240);
//     fread(header,240,1,fp);
//     // 字节反序
//     for(int i=0;i<60;i++)
//     {
//         header[i] = swap_int32(header[i]);
//     }

//     // 读取数据
//     data = (float *)malloc(3000*4);
//     fread(data,3000*4,1,fp);

//     // 字节反序
//     for(int i=0;i<3000;i++)
//     {
//         data[i] = swap_float(data[i]);
//     }

//     // 输出文件头
//     for(int i=0;i<60;i++)
//     {
//         fprintf(headerOut,"%d\n",header[i]);
//     }

//     float *dataTo;
//     dataTo = (float *)malloc(3000*4);
    
//     ibm2ieee(data,dataTo,3000,1);

//     // 输出数据
//     for(int i=0;i<3000;i++)
//     {
//         fprintf(dataOut,"%f\n",dataTo[i]);
//     }
//     return 0;
// }


int swap_int32(int i)
{
    int t = 0;
    t |= (i & 0x000000ff) << 24;
    t |= (i & 0x0000ff00) << 8;
    t |= (i & 0x00ff0000) >> 8;
    t |= (i & 0xff000000) >> 24;
    return t;
}

float swap_float(float f)
{
    float t = 0;
    char *p1 = (char *)&f;
    char *p2 = (char *)&t;
    p2[0] = p1[3];
    p2[1] = p1[2];
    p2[2] = p1[1];
    p2[3] = p1[0];
    return t;
}


short swap_int16(short s)
{
    short t = 0;
    t |= (s & 0x00ff) << 8;
    t |= (s & 0xff00) >> 8;
    return t;
}

double swap_double(double d)
{
    double t = 0;
    char *p1 = (char *)&d;
    char *p2 = (char *)&t;
    p2[0] = p1[7];
    p2[1] = p1[6];
    p2[2] = p1[5];
    p2[3] = p1[4];
    p2[4] = p1[3];
    p2[5] = p1[2];
    p2[6] = p1[1];
    p2[7] = p1[0];
    return t;
}

long long swap_int64(long long l)
{
    long long t = 0;
    t |= (l & 0x00000000000000ff) << 56;
    t |= (l & 0x000000000000ff00) << 40;
    t |= (l & 0x0000000000ff0000) << 24;
    t |= (l & 0x00000000ff000000) << 8;
    t |= (l & 0x000000ff00000000) >> 8;
    t |= (l & 0x0000ff0000000000) >> 24;
    t |= (l & 0x00ff000000000000) >> 40;
    t |= (l & 0xff00000000000000) >> 56;
    return t;
}

// 无符号int字节反序
unsigned int swap_uint32(unsigned int i)
{
    unsigned int t = 0;
    t |= (i & 0x000000ff) << 24;
    t |= (i & 0x0000ff00) << 8;
    t |= (i & 0x00ff0000) >> 8;
    t |= (i & 0xff000000) >> 24;
    return t;
}

// 无符号short字节反序
unsigned short swap_uint16(unsigned short s)
{
    unsigned short t = 0;
    t |= (s & 0x00ff) << 8;
    t |= (s & 0xff00) >> 8;
    return t;
}

// 无符号long long字节反序
unsigned long long swap_uint64(unsigned long long l)
{
    unsigned long long t = 0;
    t |= (l & 0x00000000000000ff) << 56;
    t |= (l & 0x000000000000ff00) << 40;
    t |= (l & 0x0000000000ff0000) << 24;
    t |= (l & 0x00000000ff000000) << 8;
    t |= (l & 0x000000ff00000000) >> 8;
    t |= (l & 0x0000ff0000000000) >> 24;
    t |= (l & 0x00ff000000000000) >> 40;
    t |= (l & 0xff00000000000000) >> 56;
    return t;
}



void ibm2ieee(int *from, int *to, int n, int endian)
/***********************************************************************
ibm_to_float - convert between 32 bit IBM and IEEE floating numbers
************************************************************************
Input::
from    input vector
to      output vector, can be same as input vector
endian  byte order =0 little endian (DEC, PC's), =1 other systems
*************************************************************************
Notes:
Up to 3 bits lost on IEEE -> IBM

Assumes sizeof(int) == 4

IBM -> IEEE may overflow or underflow, taken care of by
substituting large number or zero

Only integer shifting and masking are used.
*************************************************************************
Credits: CWP: Brian Sumner,  c.1985
*************************************************************************/
{
  int fconv, fmant, i, t;

  for (i = 0; i < n; ++i) {

    fconv = from[i];

    /* if little endian, i.e. endian=0 do this */
    if (endian == 0)
      fconv = (fconv << 24) | ((fconv >> 24) & 0xff) |
          ((fconv & 0xff00) << 8) | ((fconv & 0xff0000) >> 8);

    if (fconv) {
      fmant = 0x00ffffff & fconv;
      /* The next two lines were added by Toralf Foerster */
      /* to trap non-IBM format data i.e. conv=0 data  */
      if (fmant == 0) {
        printf("mantissa is zero data may not be in IBM FLOAT Format !\n");
      }
      t = (int) ((0x7f000000 & fconv) >> 22) - 130;
      while (!(fmant & 0x00800000)) {
        --t;
        fmant <<= 1;
      }
      if (t > 254) fconv = (0x80000000 & fconv) | 0x7f7fffff;
      else if (t <= 0) fconv = 0;
      else
        fconv = (0x80000000 & fconv) | (t << 23)
            | (0x007fffff & fmant);
    }
    to[i] = fconv;
  }
  return;
}

void ieee2ibm(int *from, int *to, int n, int endian)
/**********************************************************************
 float_to_ibm - convert between 32 bit IBM and IEEE floating numbers
***********************************************************************
Input:
from   input vector
n      number of floats in vectors
endian =0 for little endian machine, =1 for big endian machines
Output:
to     output vector, can be same as input vector
***********************************************************************
Notes:
Up to 3 bits lost on IEEE -> IBM
IBM -> IEEE may overflow or underflow, taken care of by
substituting large number or zero
Only integer shifting and masking are used.
***********************************************************************
Credits:     CWP: Brian Sumner
***********************************************************************/
{
  int fconv, fmant, i, t;

  for (i = 0; i < n; ++i) {
    fconv = from[i];
    if (fconv) {
      fmant = (0x007fffff & fconv) | 0x00800000;
      t = (int) ((0x7f800000 & fconv) >> 23) - 126;
      while (t & 0x3) {
        ++t;
        fmant >>= 1;
      }
      fconv = (0x80000000 & fconv) | (((t >> 2) + 64) << 24) | fmant;
    }
    if (endian == 0)
      fconv = (fconv << 24) | ((fconv >> 24) & 0xff) |
          ((fconv & 0xff00) << 8) | ((fconv & 0xff0000) >> 8);

    to[i] = fconv;
  }
  return;
}


char **IndexProcess(FILE *fp,int segyFileNum)
/***********************************************************************
IndexProcess - Read the index file and store the file name in a char array
Chang Zhimiao, 2023.06.06
************************************************************************/
{
  int i,j,k;
  char **index;
  index = (char **)malloc(segyFileNum*sizeof(char *)); 
  for(i=0;i<segyFileNum;i++)
  {
    index[i] = (char *)malloc(2048*sizeof(char));
  }

  for(i=0;i<segyFileNum;i++)
  {
    fgets(index[i],2048,fp);
    index[i][strlen(index[i])-1]='\0';
  }
  return index;
}

void getShotInfo(FILE *data, int *nt, float *dt, long *ntrtotal, int *nstotal, int *ntrps_max,float *scalel,float *scalco,int endian)
{
  int i,j,k;
  int is, ir;
  int sx,sy,gx,gy;
  segy *tr;
  tr= (segy *)malloc(sizeof(segy));

  if(endian==0)
  {
    fseek(data,3600L,SEEK_SET);
  }
  else if(endian==1){
			fseek(data, 400L, SEEK_SET);
	}
  else{
    fseek(data,240L, SEEK_SET);
  }
 
  fread(tr,240,1,data);
  *scalel = (double)swap_int16(tr->scalel);
  *scalco = (double)swap_int16(tr->scalco);
  //TODO: check if scalel and scalco are positive

if(*scalco!=0)
{
   *scalco=*scalco>0?*scalco:1/fabs(*scalco);
}
else{
  *scalco=1;
}

if(*scalel!=0)
{
   *scalel=*scalel>0?*scalel:1/fabs(*scalel);
}
else{
  *scalel=1;
}

  *nt = swap_uint16(tr->ns);
  if (!getparfloat("dt", dt)) {
		if (tr->dt) { /* is dt field set? */
			*dt = (float)swap_uint16(tr->dt) / 1000000.0;
		} else { /* dt not set, assume 4 ms */
			*dt = 0.004;
			warn("tr.dt not set, assuming dt=0.004\n");
		}
	}//end get dt
  sx = swap_uint32(tr->sx);
  sy = swap_uint32(tr->sy);
  gx = swap_uint32(tr->gx);
  gy = swap_uint32(tr->gy);
  /* get shot info */
	*ntrtotal = 1;
	*nstotal = 1;
	ir = 1;
	is = 1;
	*ntrps_max=1;
  fseek(data,(*nt)*sizeof(float),SEEK_CUR);

  // 读取二进制文件直到文件结束
  while (fread(tr, 240, 1, data) == 1) {
    if (swap_uint32(tr->sx) != sx || swap_uint32(tr->sy) != sy) {
      is++;
      sx = swap_uint32(tr->sx);
      sy = swap_uint32(tr->sy);
      if (ir>(*ntrps_max)) *ntrps_max = ir;
			if ((*ntrps_max)<2) 
				err ("data is not common-shot gathor\n");
			ir = 0;
    }
      ir++;
      *ntrtotal = *ntrtotal + 1;
    
    fseek(data,(*nt)*sizeof(float),SEEK_CUR);
  }
  if (ir>(*ntrps_max)) *ntrps_max = ir;
	*nstotal = is;
}


void SwapHeader(segy *tr)
{
  //TODO:if easier way?
  
  tr->tracl = swap_uint32(tr->tracl);
  tr->tracr = swap_uint32(tr->tracr);
  tr->fldr = swap_uint32(tr->fldr);
  tr->tracf = swap_uint32(tr->tracf);
  tr->ep = swap_uint32(tr->ep);
  tr->cdp = swap_uint32(tr->cdp);
  tr->cdpt = swap_uint32(tr->cdpt);
  tr->trid = swap_uint16(tr->trid);
  tr->nvs = swap_uint16(tr->nvs);
  tr->nhs = swap_uint16(tr->nhs);
  tr->duse = swap_uint16(tr->duse);
  tr->offset = swap_uint32(tr->offset);
  tr->gelev = swap_uint32(tr->gelev);
  tr->selev = swap_uint32(tr->selev);
  tr->sdepth = swap_uint32(tr->sdepth);
  tr->gdel = swap_uint32(tr->gdel);
  tr->sdel = swap_uint32(tr->sdel);
  tr->swdep = swap_uint32(tr->swdep);
  tr->gwdep = swap_uint32(tr->gwdep);
  tr->scalel = swap_uint16(tr->scalel);
  tr->scalco = swap_uint16(tr->scalco);
  tr->sx = swap_uint32(tr->sx);
  tr->sy = swap_uint32(tr->sy);
  tr->gx = swap_uint32(tr->gx);
  tr->gy = swap_uint32(tr->gy);
  tr->counit = swap_uint16(tr->counit);
  tr->wevel = swap_uint16(tr->wevel);
  tr->swevel = swap_uint16(tr->swevel);
  tr->sut = swap_uint16(tr->sut);
  tr->gut = swap_uint16(tr->gut);
  tr->sstat = swap_uint16(tr->sstat);
  tr->gstat = swap_uint16(tr->gstat);
  tr->tstat = swap_uint16(tr->tstat);
  tr->laga = swap_uint16(tr->laga);
  tr->lagb = swap_uint16(tr->lagb);
  tr->delrt = swap_uint16(tr->delrt);
  tr->muts = swap_uint16(tr->muts);
  tr->mute = swap_uint16(tr->mute);
  tr->ns = swap_uint16(tr->ns);
  tr->dt = swap_uint16(tr->dt);
  tr->gain = swap_uint16(tr->gain);
  tr->igc = swap_uint16(tr->igc);
  tr->igi = swap_uint16(tr->igi);
  tr->corr = swap_uint16(tr->corr);
  tr->sfs = swap_uint16(tr->sfs);
  tr->sfe = swap_uint16(tr->sfe);
  tr->slen = swap_uint16(tr->slen);
  tr->styp = swap_uint16(tr->styp);
  tr->stas = swap_uint16(tr->stas);
  tr->stae = swap_uint16(tr->stae);
  tr->tatyp = swap_uint16(tr->tatyp);
  tr->afilf = swap_uint16(tr->afilf);
  tr->afils = swap_uint16(tr->afils);
  tr->nofilf = swap_uint16(tr->nofilf);
  tr->nofils = swap_uint16(tr->nofils);
  tr->lcf = swap_uint16(tr->lcf);
  tr->hcf = swap_uint16(tr->hcf);
  tr->lcs = swap_uint16(tr->lcs);
  tr->hcs = swap_uint16(tr->hcs);
  tr->year = swap_uint16(tr->year);
  tr->day = swap_uint16(tr->day);
  tr->hour = swap_uint16(tr->hour);
  tr->minute = swap_uint16(tr->minute);
  tr->sec = swap_uint16(tr->sec);
  tr->timbas = swap_uint16(tr->timbas);
  tr->trwf = swap_uint16(tr->trwf);
  tr->grnors = swap_uint16(tr->grnors);
  tr->grnofr = swap_uint16(tr->grnofr);
  tr->grnlof = swap_uint16(tr->grnlof);
  tr->gaps = swap_uint16(tr->gaps);
  tr->otrav = swap_uint16(tr->otrav);
  tr->d1 = swap_uint16(tr->d1);
  tr->f1 = swap_uint16(tr->f1);
  tr->d2 = swap_uint16(tr->d2);
  tr->f2 = swap_uint16(tr->f2);
  tr->ungpow = swap_uint16(tr->ungpow);
  tr->unscale = swap_uint16(tr->unscale);
  tr->ntr = swap_uint16(tr->ntr);
  tr->mark = swap_uint16(tr->mark);
  tr->shortpad = swap_uint16(tr->shortpad);
  tr->unass[0] = swap_uint16(tr->unass[0]);
  tr->unass[1] = swap_uint16(tr->unass[1]);
  tr->unass[2] = swap_uint16(tr->unass[2]);
  tr->unass[3] = swap_uint16(tr->unass[3]);
  tr->unass[4] = swap_uint16(tr->unass[4]);
  tr->unass[5] = swap_uint16(tr->unass[5]);
  tr->unass[6] = swap_uint16(tr->unass[6]);
  tr->unass[7] = swap_uint16(tr->unass[7]);
  tr->unass[8] = swap_uint16(tr->unass[8]);
  tr->unass[9] = swap_uint16(tr->unass[9]);
  tr->unass[10] = swap_uint16(tr->unass[10]);
  tr->unass[11] = swap_uint16(tr->unass[11]);
  tr->unass[12] = swap_uint16(tr->unass[12]);
  tr->unass[13] = swap_uint16(tr->unass[13]);
}










