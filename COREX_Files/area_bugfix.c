/*----------------------------------------*
 * routines from ACCESS to be reimplemented
 * as asa_lib
 *----------------------------------------*/

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "libasa.h"

#define CX	x
#define CY	y
#define CZ	z
#define RAD	radius
#define NRAD	area

static int findneighbors();
static double findarea();
static double calcarea();
static int sortag();

/*
 * The following are dimensioned to the max no of atoms(nsatm)
 */
#define MAXATOMS 10000

static int *CUBE;

/*
 * The following are dimensioned to the max no of intersections of
 * neighbouring spheres (ICT)
 */

#define ICT	2400

static int	*INOV;
static double   *DX, *DY, *DSQ;
static double	*ARCI, *ARCF;
static int	*TAG;

/* OC:
 * NCUBE =  max no of cubes
 * NAC   =  max no of atoms in cube
 */

#define NCUBE	8000
#define	NAC	36

int	*ITAB;
int	**NATM;
int	CURRENT_NAC, OLD_NAC;
int	IDIM, JIDIM, KJIDIM;

/*
 * Inital areas are stored in the occupancy field of each atom in a PDB
 * record, the final answer is stored in the temp_factor field. The number
 * of records is passed to the functions that need it.
 */

extern asa_atom_rec_p a_pdb;
extern char *err_string;

/*----------------------------------------*/

int surfinit(ANO, RH2O)
int ANO;
double RH2O;
{

	int I, J, KJI, N;

	int S, T, U;

    double XMIN, YMIN, ZMIN, XMAX, YMAX, ZMAX, RMAX;

    XMIN =  9999.0,
    YMIN =  9999.0,
    ZMIN =  9999.0,
    XMAX = -9999.0,
    YMAX = -9999.0,
    ZMAX = -9999.0;

    RMAX = 0.0;

		if (!(
		(CUBE = (int *) malloc((MAXATOMS + 1) * sizeof(int))) &&
		(INOV = (int *) malloc((ICT + 1) * sizeof(int))) &&
		(DX = (double *) malloc((ICT + 1) * sizeof(double))) &&
		(DY = (double *) malloc((ICT + 1) * sizeof(double))) &&
		(DSQ = (double *) malloc((ICT + 1) * sizeof(double))) &&
		(ARCI = (double *) malloc((ICT + 1) * sizeof(double))) &&
		(ARCF = (double *) malloc((ICT + 1) * sizeof(double))) &&
		(TAG = (int *) malloc((ICT + 1) * sizeof(int))) 
		)) {
			strcpy(err_string, "asa_lib: malloc() failure.");
			return 0;
		}

	  for (I = 0; I < ANO; I++) {

		  a_pdb[I].RAD = a_pdb[I].RAD + RH2O;

		  if (a_pdb[I].RAD > RMAX)	RMAX = a_pdb[I].RAD;

		  if (XMIN > a_pdb[I].CX)	XMIN = a_pdb[I].CX;

		  if (YMIN > a_pdb[I].CY)	YMIN = a_pdb[I].CY;

		  if (ZMIN > a_pdb[I].CZ)	ZMIN = a_pdb[I].CZ;

		  if (XMAX < a_pdb[I].CX)	XMAX = a_pdb[I].CX;

		  if (YMAX < a_pdb[I].CY)	YMAX = a_pdb[I].CY;

		  if (ZMAX < a_pdb[I].CZ)	ZMAX = a_pdb[I].CZ;
	  }

	  RMAX = RMAX * 2.000;

	  /*
	   * Cubicals containing the atoms are setup. The dimension
	   * of an edge equals the radius of the largest atom sphere;
	   * the cubes have a single index.
	   */

	  IDIM = (XMAX - XMIN) / RMAX + 1.000;
	  if (IDIM < 3)
		  IDIM = 3;

	  JIDIM =  (YMAX - YMIN) / RMAX + 1.000;

	  if (JIDIM < 3)
		  JIDIM = 3;

	  JIDIM = IDIM * JIDIM;

	  KJIDIM = (ZMAX - ZMIN) / RMAX + 1.000;

	if (KJIDIM < 3)
		KJIDIM = 3;

	KJIDIM = JIDIM * KJIDIM;

	if (!(ITAB = (int *) malloc((KJIDIM + 1) * sizeof(int)))) {
	   	strcpy(err_string, "asa_lib: malloc() failure.");
	    return 0;
	}
    
	for (I = 0; I < KJIDIM + 1; I++)
		ITAB[I] = 0;
#ifdef DEBUG
printf("RMAX=%f  IDIM=%d  JIDIM=%d  KJIDIM=%d\n", RMAX, IDIM, JIDIM, KJIDIM);
#endif

	if (!(NATM = (int **) malloc((NAC + 1) * sizeof(int *)))) {
	    strcpy(err_string, "asa_lib: malloc() failure.");
	    return 0;
	}

/* does the zero element get used?  I doubt it. */

	for (I = 0; I < NAC + 1; I++) {
		if(!(NATM[I] = (int *) malloc((KJIDIM + 1) * sizeof(int)))) {
		    strcpy(err_string, "asa_lib: malloc() failure.");
		    return 0;
                }
		for (J = 0; J < KJIDIM + 1; J++)
			NATM[I][J] = 0;
	}
	

	CURRENT_NAC = NAC;

	/*
	 * Prepare upto NCUBE cubes each containing upto ANO atoms.
	 * The cube index is KJI. The atom index for each cube is in ITAB.
	 */

	for (I = 0; I < ANO; I++) {
		S = (a_pdb[I].CX - XMIN) / RMAX + 1.000;
		T = (a_pdb[I].CY - YMIN) / RMAX;
		U = (a_pdb[I].CZ - ZMIN) / RMAX;

		KJI = U * JIDIM + T * IDIM + S;
		N = ITAB[KJI] + 1;

		if (N > CURRENT_NAC) {
#ifdef DEBUG
printf("Had to double NAC. N=%d, I=%d\n", N, I);
#endif
			OLD_NAC = CURRENT_NAC;
			CURRENT_NAC *= 2;
			if (!(NATM = (int **) realloc(NATM, (CURRENT_NAC + 1) * sizeof(int *)))) {
		  		  strcpy(err_string, "asa_lib: realloc() failure.");
		  		  return 0;
            }

			for (I = OLD_NAC+1 ; I < CURRENT_NAC + 1 ; I++) {
				if (!(NATM[I] = (int *) malloc((KJIDIM + 1) * sizeof(int)))) {
		  		  strcpy(err_string, "asa_lib: malloc() failure.");
		  		  return 0;
            	}
				for (J = 0; J < KJIDIM + 1; J++)
					NATM[I][J] = 0;
			}
		}
			

		ITAB[KJI] = N;
		NATM[N][KJI] = I+1;
		CUBE[I+1] = KJI;

	}

	return 1;
}

/*----------------------------------------*/

void surfcleanup()
{
    int I;

    free(CUBE);
    free(INOV);
    free(DX);
    free(DY);
    free(DSQ);
    free(ARCI);
    free(ARCF);
    free(TAG);
    free(ITAB);
    for (I = 0; I < CURRENT_NAC + 1; I++)
		free(NATM[I]);
    free(NATM);
}

/*----------------------------------------*/

/* I can't imagine that this is anything but: 
      double surfarea(num_atoms, slice_width, solvent_radius); */

double surfarea(ANO, p, RH2O)
int ANO;
double p, RH2O;
{

	int IR, IO, NZP;

	double RR, RRX2, RRSQ, XR, YR, ZR;

	double ZRES, ZGRID, AREA;

	double findarea();
	
	double totalarea = 0.0;

	/* Process each atom. */

	for (IR = 0; IR < ANO; IR++) {

		XR = a_pdb[IR].CX;
		YR = a_pdb[IR].CY;
		ZR = a_pdb[IR].CZ;
		RR = a_pdb[IR].RAD;
		RRX2 = a_pdb[IR].RAD * 2.000;
		RRSQ = a_pdb[IR].RAD * a_pdb[IR].RAD;


		IO = findneighbors(IR+1, XR, YR, ZR);

		if (IO != 0) {

			/* OC: Z resolution determined */
			NZP = 1.0 / p + 0.500;
	
			ZRES = RRX2 / NZP;
			ZGRID = a_pdb[IR].CZ - RR - ZRES / 2.000;

			AREA = findarea(IR+1, IO, NZP, ZRES, ZGRID, ZR, RRSQ);
			
		} else

			AREA = TWOPI * RRX2;
		
		/* This would Scale area to VDW shell of atoms
		a_pdb[IR].NRAD = AREA * SQ(RR - RH2O) / RR;
		   while this is _accessible_ surface */
		a_pdb[IR].NRAD = AREA * RR;

		totalarea += a_pdb[IR].NRAD;
	}

	return totalarea;

}

/*----------------------------------------*/

static int findneighbors(IR, XR, YR, ZR)
  int IR;
  double XR, YR, ZR;
{
	
	int	I, J, K, M, MKJI, IN, NM;
	
	int IO = 0;

	/* Find the 'MKJI' cubes neighboring the KJI cube */
	for (K = 1; K <= 3; K++) {
		for (J = 1; J <= 3; J++) {
			for (I = 1; I <= 3; I++) {

				MKJI = CUBE[IR] + (K - 2) *
					JIDIM +	(J - 2) * IDIM + (I - 2);

				if (MKJI < 1)
					continue;

				if (MKJI > KJIDIM)
					return IO;

				NM = ITAB[MKJI];
				
				if (NM < 1)
					continue; 
			        /* OC:
				 * Record the atoms in INOV that neighbor
				 * atom IR.
				 */

				for (M = 1; M <= NM; M++) {
					IN = NATM[M][MKJI];

					if (IN == IR)
						continue;

					IO += 1;										
					if (IO >= ICT)
						asa_error(EX_SOFTWARE,
						      "Overflow in IO");
					DX[IO] = XR - a_pdb[IN-1].CX;
					DY[IO] = YR - a_pdb[IN-1].CY;
					DSQ[IO] = SQ(DX[IO]) + SQ(DY[IO]);
					INOV[IO] = IN;
					


				}
				
			}
		}
	}
	return IO;
}

/*----------------------------------------*/

static double findarea(IR, IO, NZP, ZRES, ZGRID, ZR, RRSQ)
  int IR, IO, NZP;
  double ZRES, ZGRID, ZR, RRSQ;
{
	int I, IN, J, KARC;

	double RSEC2R, RSECR, RSECN, RSEC2N;

	double ALPHA, BETA, TI, TF, B, D_TEMP;

	double AREA = 0.0;

	double calcarea();

	/* OC: Section atom spheres perpendicular to the Z axis */
	for (I = 1; I <= NZP; I++) {
		ZGRID = ZGRID + ZRES;
		
		/*
		 * Find the radius of the circle of intersection of
		 * the IR sphere on the current Z-plane
		 */
		
		RSEC2R = RRSQ - SQ(ZGRID - ZR);

		RSECR = sqrt(RSEC2R);
		
		for (J = 1; J <= ICT; J++)
			ARCI[J] = 0.0;
		/* 34 CONTINUE */

		KARC = 0;

		for (J = 1; J <= IO; J++) {	/* DO 10 J=1,IO */
			IN = INOV[J];
			
			/* OC: Find radius of circle locus */
			RSEC2N = SQ(a_pdb[IN-1].RAD) 
				- SQ(ZGRID - a_pdb[IN-1].CZ);
			
			if (RSEC2N <= 0.0)
				continue;
			else
				RSECN = sqrt(RSEC2N);

			/*
			 * Find intersections of N circles with IR
			 * circles in section.
			 */

			D_TEMP = sqrt(DSQ[J]);

			if (D_TEMP >= RSECR + RSECN)
				continue;
			
			/*
			 * Do the circles intersect, or is one circle
			 * completely inside the other?
			 */

			B = RSECR - RSECN;

			if (D_TEMP > ABS(B))  {
				KARC = KARC + 1;
			} else if (B <= 0.0) {
				goto nine;
			} else
				continue;
			
			/*
			 * If the circles intersect, find the points
			 * of intersection.
			 */
			
			if (KARC >= ICT)
				asa_error(EX_SOFTWARE, "overflow in ICT");
			
			/*
			 * Initial and final arc endpoints are found
			 * for the IR circle intersected by a
			 * neighboring circle contained in the same
			 * plane. The initial endpoint of the
			 * enclosed arc is stored in ARCI, and the
			 * final arc in ARCF (law of cosines).
			 */

			ALPHA = acos((DSQ[J] + RSEC2R - RSEC2N)
				     / (2.000 * D_TEMP * RSECR));

			/*
			 * Alpha is the angle between a line
			 * containing a point of intersection and the
			 * reference circle center and the line
			 * containing both circle centers
			 */
			BETA = atan2(DY[J], DX[J]) + PI;
			
			/*
			 * Beta is the angle between the line
			 * containing both circle centers and the
			 * x-axis
			 */
			TI = BETA - ALPHA;
			TF = BETA + ALPHA;
			if (TI < 0.0)
				TI = TI + TWOPI;
			if (TF > TWOPI)
				TF = TF - TWOPI;
			ARCI[KARC] = TI;
			if (TF < TI) {
				/*
				 * If the arc crosses zero, then it
				 * is broken into two segments. the
				 * first ends at 2PI and the second
				 * begins at zero.
				 */
				ARCF[KARC] = TWOPI;
				KARC = KARC + 1;
			}
			ARCF[KARC] = TF;
			
			
		}
		AREA += calcarea(KARC, ZRES);
  nine:
			;

	}
	
	return AREA;
	
}

/*----------------------------------------*/

static double calcarea(KARC, ZRES)
int KARC;
double ZRES;
{
	int J;
	
	double ARCSUM, T, TT;

        if (KARC == 0) return TWOPI*ZRES;               /***** BUG FIX *****/

	/* OC:
	 * Find the accessible contact surface area for the
	 * sphere IR on this section
	 */

	/* OC:
	 * The arc endpoints are sorted on the value
	 * of the initial arc endpoint.
	 */
	
	/* Set up the tag array */
	for ( J = 1 ; J <= KARC ; J++)
		TAG[J] = J;
	
	sortag(ARCI, KARC, TAG);
	
	/* OC: Calculate the accessible area */
	ARCSUM = ARCI[1];
	
	T = ARCF[TAG[1]];
	
	
	for (J = 2; J <= KARC ; J++) {	/* DO 27 K=2,KARC */
		if (T < ARCI[J])
			ARCSUM = ARCSUM + ARCI[J] - T;
		TT = ARCF[TAG[J]];
		if (TT > T)
			T = TT;
	}	/* 27 CONTINUE */
	
	ARCSUM = ARCSUM + TWOPI - T;

	/* OC:
	 * The area/radius is equal to the accessible
	 * arc length x the section thickness.
	 */
	
	return (ARCSUM * ZRES);
	
}

/*----------------------------------------*/

static int sortag(A,N,TAG)
      double A[];
      int N, TAG[];
{
      int I, J, L, M, K, IJ, TG;
      double T = 0.0, TT = 0.0;
      static int IU[16],IL[16];

      M=1; I=1; J= N;
five:
      if(I >= J) goto seventy;

ten:
      K=I;
      IJ=(J+I)/2;
      T=A[IJ];
      if(A[I] <= T) goto twenty;
      A[IJ]= A[I];
      A[I]=T;
      T=A[IJ];
      TG=TAG[IJ];
      TAG[IJ]=TAG[I];
      TAG[I]=TG;
twenty:
      L=J;
      if(A[J] >= T) goto fourty;
      A[IJ]=A[J];
      A[J]=T;
      T=A[IJ];
      TG=TAG[IJ];
      TAG[IJ]=TAG[J];
      TAG[J]=TG;
      if(A[I] <= T) goto fourty;
      A[IJ]=A[I];
      A[I]=T;
      T=A[IJ];
      TG=TAG[IJ];
      TAG[IJ]=TAG[I];
      TAG[I]=TG;
      goto fourty;
thirty:
      A[L]=A[K];
      A[K]=TT;
      TG=TAG[L];
      TAG[L]=TAG[K];
      TAG[K]=TG;
fourty:
      L=L-1;
      if(A[L] > T) goto fourty;
      TT=A[L];
fifty:
      K=K+1;
      if(A[K] < T) goto fifty;
      if(K <= L) goto thirty;
      if(L-I <= J-K) goto sixty;
      IL[M]=I;
      IU[M]=L;
      I=K;
      M=M+1;
      goto eighty;
sixty:
      IL[M]=K;
      IU[M]=J;
      J=L;
      M=M+1;
      goto eighty;
seventy:
      M=M-1;
      if(M == 0) return 1;
      I=IL[M];
      J=IU[M];
eighty:
      if(J-I >= 1) goto ten;
      if(I == 1) goto five;
      I=I-1;
ninety:
      I=I+1;
      if(I == J) goto seventy;
      T=A[I+1];
      if(A[I] <= T) goto ninety;
      TG=TAG[I+1];
      K=I;
hundred:
      A[K+1]=A[K];
      TAG[K+1]=TAG[K];
      K=K-1;
      if(T < A[K]) goto hundred;
      A[K+1]=T;
      TAG[K+1]=TG;
      goto ninety;

}

/*----------------------------------------*/

