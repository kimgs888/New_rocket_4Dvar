/**************************************************************************************
File:    trace.h
Author:  Paul S.J. Spencer
Purpose: Ray tracing through a distorted cartesian vertex grid
---------------------------------------------------------------------------------------

Method:
   Voxels in a grid of points are each decomposed into 6 tetrahedra. The ray trace
   is solved analytically in each tetrahedron assuming tri-linear gradients in the
   field variable f. The relationship between the field variable and the refractive
   index, RI, at each point is implemented for the following cases;


   Ray path      RI functional   Example cases
   descriptor    form
   ----------    --------------  ----------------------------------------
   LINEAR:       RI   = 1        GHz propagation, straight line ray paths
   CIRCULAR:     RI   = A/f      Seismic/acoustic rays, f represents velocity
   PARABOLIC:    RI^2 = 1 - A*f  HF/VHF, f represents electron density

   A = User defined constant,

Voxel decomposition into six tetrahedra:
       7 ________6
        |\        \       Forms tetrahedra: 5671, 6217, 1237
        | \        \                        4571, 0147, 0137
        | 4\_______5\
        |  |        |     y   z
       3\  |     2  |      \  |
         \ |        |       \ |
          \|________|        \|_______x
           0        1

**************************************************************************************/
#ifndef TRACE_H
#define TRACE_H

#include <math.h>
#include <string.h>

#include "matrices.h"
#include "types.h"
#include "error.h"

/* Limiting values for floating point arithmetic */
#define SMALL 1E-12
#define LARGE 1E+12

/* Maximum path intersections per ray 1200*/
#define MAXPATH 7200

/* Maximum malloc block, bytes */
/* changed by nick from 100000000: this bcz for smaller height res the memo was to small*/
/* so only a few rays could be out put: this new value is the max I could get when changing*/
/* beyond that value the code crushes  320910000*/
#define MAXMEM 320910000


/* Maximum itterations for shooting given Rx, Rv */
#define TIMEOUT 2

/* Vector dot product */
#define DOT(A,B) ((A).x*(B).x + (A).y*(B).y + (A).z*(B).z)
/* Vector magnitude   */
#define MAG(A) (sqrt( fabs((A).x*(A).x + (A).y*(A).y + (A).z*(A).z) +1E-30))
/* Length of vector between points */
#define LEN(A,B) (sqrt(pow((A).x-(B).x,2)+pow((A).y-(B).y,2)+pow((A).z-(B).z,2)))
/* Maximum of two numbers */
#define MAX(A,B) ((A)>(B)?(A):(B))

/* Sign */
#define SIGN(A) ((A)>=0?1:-1)
/* Vector cross product components */
#define CROSSX(A,B) ((A).y*(B).z - (A).z*(B).y)
#define CROSSY(A,B) ((A).z*(B).x - (A).x*(B).z)
#define CROSSZ(A,B) ((A).x*(B).y - (A).y*(B).x)
/* Total number of vertices */
#define SIZE(Grid) ((Grid)->Dim[0]*(Grid)->Dim[1]*(Grid)->Dim[2])
/* Ray exits grid in x, radial */
#define EXITXA(Tx,Fx) (((Tx)==4&&(Fx)==1)||((Tx)==5&&(Fx)==1))
#define EXITXB(Tx,Fx) (((Tx)==0&&(Fx)==1)||((Tx)==1&&(Fx)==0))
#define EXITX( Tx,Fx) ( EXITXA((Tx),(Fx))||EXITXB((Tx),(Fx)) )

/* Copy ray from RAY S to RAY D */
#define COPYRAY(D,S) (memcpy(&(D),&(S),sizeof(RAY)))


typedef struct CPOINT     /* Point in 3D space                               */
{
   double x;
   double y;
   double z;
} CPOINT;

typedef struct VECTOR     /* Vector in 3D space                              */
{
   double x;
   double y;
   double z;
} VECTOR;

typedef struct GRID       /* Grid definition structure                       */
{
   double *X, *Y, *Z;     /* 3D Cartesian vertex coordinates, as from ndgrid */
   long   Dim[3];         /* Size of dimensions of X,Y,Z vertex arrays       */
   double *F;             /* Scalar field at grid vertices                   */
   int    Wrap;           /* Third dimension wrapping, set internally        */
} GRID;

typedef struct VOXEL      /* 3D voxel vertices ordered as in figure above    */
{
   CPOINT Vertex[8];      /* Ordered vertices                                */
   double Field[8];       /* Scalar field value at each vertex               */
   double OPath[8];       /* Contribution to optical path for each vertex    */
/*   CPOINT Sect[8];        /* Ray intersection coordinates with tetra edges
     UINT   NumSect;        /* Number of ray intersections in Sect             */
} VOXEL;

typedef enum PROPAGATOR   /* Raytrace propagator                             */
{
   LINEAR,                /* RI = 1                                          */
   PARABOLIC,             /* Parabolic phase path, RI^2 = 1 + A*f            */
   PARABOLIC_GROUP,       /* Parabolic group path, RI^2 = 1 + A*f            */
   CIRCULAR               /* RI   = A / f       */
} PROPAGATOR;


typedef struct RAY        /* Ray propogation parameters for single ray       */
{
   PROPAGATOR Prop;       /* Raytrace propagator                             */
   double     A;          /* Refractive constant (see PROPAGATOR)            */
   CPOINT     p;          /* Ray position                                    */
   VECTOR     v;          /* Ray vector                                      */
   long       Ix[3];      /* Entry point voxel index                         */
   double     Path[8];    /* Optical path contributions per vertex for voxel */
   double*    Op;         /* Pointer to per vertex optical path              */
/* Itterative shooting parameters                                            */
   CPOINT     Rx;         /* Receiver position                               */
   VECTOR     Rv;         /* Receiver unit vector                            */
   double     Rm;         /* Magnitude of Rv before normalisation            */
} RAY;

typedef struct TRACE      /* Ray tracing parameters passed from mex function */
{
   PROPAGATOR Prop;       /* Raytrace propagator                             */
   PMATRIX A;             /* Propagation constant                            */
   PMATRIX Tx;            /* Transmitter cartesian coords                    */
   PMATRIX Tv;            /* Transmitter transmission vector                 */
   PMATRIX Rx;            /* Receiver cartesian coords                       */
   PMATRIX Rv;            /* Receiver reception vector                       */
   PMATRIX Ix;            /* Transmitter entry point voxel indices           */
} TRACE;


typedef struct PATH       /* Store decomposed path integral matrix           */
{
   int     rows;          /* Number of rows                                  */
   int     cols;          /* Number of columns                               */
   int     size;          /* Maximum entries                                 */
   int     i;             /* Current entry                                   */
   double *row;           /* Row indices                                     */
   double *col;           /* Col indices                                     */
   double *val;           /* Values                                          */
} PATH;

typedef struct SECT       /* Store ray voxel intersections                   */
{
   int     size;          /* Maximum entries                                 */
   int     i;             /* Current entry                                   */
   double *val;           /* Values                                          */
} SECT;



/* Function prototypes */
UINT   InTetra( VOXEL *Voxel, RAY *Ray );
long   Vox(GRID *,long,long,long);
int    InGrid(GRID *,long,long,long);
ERRORS Trace(VOXEL *Voxel,RAY *Ray,UINT *Tx, UINT *Face, double *Dx);
int    Gradient( VOXEL *Voxel, UINT Tx, VECTOR *G );
double LinePlane( CPOINT *A, CPOINT *B, VECTOR *n, CPOINT *P );
double LineTriangle(CPOINT *,VECTOR *,CPOINT *,CPOINT *,CPOINT *,double *, CPOINT *);
int    TrianglePlane(CPOINT *,CPOINT *,CPOINT *,VECTOR *,CPOINT *,CPOINT *,CPOINT *);
int    Parametric(CPOINT  *,CPOINT  *,VECTOR  *,VECTOR  *,double  *,double  *);
int    Roots( double, double, double, double Root[2] );
int    Inverse3x3( double (*In)[3], double (*Out)[3] );
double Determinant3x3( double *A, double *B, double *C );
int    PlanePlane(VECTOR *,VECTOR *,CPOINT *,CPOINT *,CPOINT *,CPOINT *);
ERRORS RayTrace(GRID *,TRACE *,PATH *,SECT *);
ERRORS SetVoxel(GRID *,long Vx[3],VOXEL *);
ERRORS Shoot(GRID *,RAY  *,PATH *,SECT *);
int    ParametricB(CPOINT *,CPOINT *,CPOINT *,CPOINT *,double *,double *);

#endif /* TRACE_H */
