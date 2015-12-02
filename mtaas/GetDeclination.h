#include "stdafx.h"

#define M_PI PI

#define NaN log(-1.0)

#ifndef SEEK_SET
#define SEEK_SET 0
#define SEEK_CUR 1
#define SEEK_END 2
#endif

#define IEXT 0
#define FALSE 0
#define TRUE 1                  /* constants */
#define RECL 81

#define MAXINBUFF RECL+14

/** Max size of in buffer **/

#define MAXREAD MAXINBUFF-2
/** Max to read 2 less than total size (just to be safe) **/

#define MAXMOD 30
/** Max number of models in a file **/

#define PATH MAXREAD
/** Max path and filename length **/

#define EXT_COEFF1 (float)0
#define EXT_COEFF2 (float)0
#define EXT_COEFF3 (float)0

#define MAXDEG 13
#define MAXCOEFF (MAXDEG*(MAXDEG+2)+1) /* index starts with 1!, (from old Fortran?) */

class GetDeclination
{
public:
	GetDeclination():
		declination_input(0),
		x(0),
		C(0),
		z(0),
		modelfile("IGRF11.cof.txt") {}
	static float declination(float latitude, float longitude, float elev, int year, int month, int day);
	const double declination_input;// = 0.1;

private:
	char *modelfile;
	float gh1[MAXCOEFF];
	float gh2[MAXCOEFF];
	float gha[MAXCOEFF];              /* Geomag global variables */
	float ghb[MAXCOEFF];
	float x, C, z;//x=0,C=0,z=0;
	float xtemp, ytemp, ztemp;
	FILE *stream; // = NULL;                /* Pointer to specified model data file */

	static float julday(int i_month, int i_day, int i_year);
	static int calc_geomagXYZ(char* mdfile,float latitude, float longitude,float alt, float sdate, float *X, float *Y, float *Z);
	static int getshc(char* file,int iflag,long strec,int nmax_of_gh,int gh);
	static int extrapsh(float date,float dte1,int nmax1,int nmax2,int gh);
	static int interpsh(float date,float dte1,int nmax1,float dte2,int nmax2,int gh);
	static int shval3(int igdgc,float flat,float flon,float elev,int nmax,int gh,int iext,float ext1,float ext2,float ext3);
};