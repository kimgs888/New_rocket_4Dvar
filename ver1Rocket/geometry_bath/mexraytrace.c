/*****************************************************************************
Author:     PSS
File:       MEXRayTrace.c
Purpose:    MEX interface for ray tracing
Notes:      See RAYTRACE.m for more details.
*****************************************************************************/

#include <string.h>
#include <math.h>
#include "mex.h"
#include "matrix.h"
/* #include "stdlib.h"  added this after mex integer cast error */
#include "trace.h"
/* #include "matrices.h" */



void mexFunction
(
   int nlhs,               /* Number of LHS arguments    */
   mxArray *plhs[],        /* Pointer to LHS arguments   */
   int nrhs,               /* Number of RHS arguments    */
   const mxArray *prhs[]   /* Pointer to RHS arguments   */
)
{
   GRID     Grid;          /* Grid definition structure                  */
   TRACE    Trace;         /* Ray trace structure                        */
   SECT    *Sect = NULL;   /* Intersections with grid                    */
   PATH    *Path = NULL;   /* Path matrix                                */
   char    *Info;
   long     Size;
   char     Str[100];
   int     *pInt, i;
   mxArray *Tmp;
   mxArray *pTmp[5];
   UINT     Err, F=0;

/* Check the number of input and output of arguments, print usage instructions */
   if( nrhs != 2 )
   {
      if( (Info = mxCalloc(200,1)) == NULL ) mexErrMsgTxt("Calloc"); Info[0]=0;
      strcat(Info,"Error using MEX function raytrace\n");
      strcat(Info,"Incorrect number of input arguments\n");
      mexErrMsgTxt(Info);
   }


/* Extract raytrace structure fields */
   if( (Tmp = mxGetField( prhs[1], 0, "Prop")) == NULL )
      mexErrMsgTxt("raytrace.c: Ray.Prop field is empty");
   if( !mxIsChar(Tmp))
      mexErrMsgTxt("raytrace.c: Ray.Prop is not a character array");
   if( mxGetString(Tmp,Str,99) != 0 )
      mexErrMsgTxt("raytrace.c: Ray.Prop string is too long");
   if(      strcmp( Str, "CIRCULAR" )        == 0 ) Trace.Prop = CIRCULAR;
   else if( strcmp( Str, "PARABOLIC" )       == 0 ) Trace.Prop = PARABOLIC;
   else if( strcmp( Str, "PARABOLIC_GROUP" ) == 0 ) Trace.Prop = PARABOLIC_GROUP;
   else if( strcmp( Str, "LINEAR" )          == 0 ) Trace.Prop = LINEAR;
   else mexErrMsgTxt("raytrace.c: Ray.Prop unrecognised");

   if( (Tmp = mxGetField(prhs[1],0,"A")) == NULL || !mxIsDouble(Tmp) )
      mexErrMsgTxt("raytrace.c: No Ray.A field or wrong type");
   LoadPacked( &Trace.A, mxGetM(Tmp) ,mxGetN(Tmp) , mxGetPr(Tmp) );

   if( (Tmp = mxGetField(prhs[1],0,"Tx")) == NULL || !mxIsDouble(Tmp) )
      mexErrMsgTxt("raytrace.c: No Ray.Tx field or wrong type");
   LoadPacked( &Trace.Tx, mxGetM(Tmp) ,mxGetN(Tmp) , mxGetPr(Tmp) );

   if( (Tmp = mxGetField(prhs[1],0,"Tv")) == NULL || !mxIsDouble(Tmp) )
      mexErrMsgTxt("raytrace.c: No Ray.Tv field or wrong type");
   LoadPacked( &Trace.Tv, mxGetM(Tmp) ,mxGetN(Tmp) , mxGetPr(Tmp) );

   if( (Tmp = mxGetField(prhs[1],0,"Ix")) == NULL || !mxIsDouble(Tmp) )
      mexErrMsgTxt("raytrace.c: No Ray.Ix field or wrong type");
   LoadPacked( &Trace.Ix, mxGetM(Tmp) ,mxGetN(Tmp) , mxGetPr(Tmp) );

   if( (Tmp = mxGetField(prhs[1],0,"Rx")) != NULL  )
   {
      if( !mxIsDouble(Tmp) ) mexErrMsgTxt("raytrace.c: Ray.Rx wrong type");
      LoadPacked( &Trace.Rx, mxGetM(Tmp) ,mxGetN(Tmp) , mxGetPr(Tmp) );
   }
   else InitPacked( &Trace.Rx,0,0,NULL);

   if( (Tmp = mxGetField(prhs[1],0,"Rv")) != NULL )
   {
      if( !mxIsDouble(Tmp) ) mexErrMsgTxt("raytrace.c: Ray.Rv wrong type");
      LoadPacked( &Trace.Rv, mxGetM(Tmp) ,mxGetN(Tmp) , mxGetPr(Tmp) );
   }
   else InitPacked( &Trace.Rv,0,0,NULL);

   if( (Trace.Tx.Rows != Trace.Ix.Rows) ||
       (Trace.A.Rows  > 1 && Trace.A.Rows  != MAX(Trace.Tx.Rows,Trace.Tv.Rows)) ||
       MAX(Trace.Tx.Rows,Trace.Tv.Rows) < MAX(Trace.Rx.Rows,Trace.Rv.Rows) )
      mexErrMsgTxt("raytrace.c: Incorrect number of rows for field in Ray structure");


/* Extract Grid structure fields */
   if( (Tmp = mxGetField(prhs[0],0,"X")) == NULL || !mxIsDouble(Tmp) )
      mexErrMsgTxt("raytrace.c: No Grid.X field or wrong type");
   if( mxGetNumberOfDimensions(Tmp) != 3 )
      mexErrMsgTxt("raytrace.c: Grid.X must be 3 dimensional");
   pInt = (int *)mxGetDimensions(Tmp);
   Grid.Dim[0] = pInt[0]; Grid.Dim[1] = pInt[1]; Grid.Dim[2] = pInt[2];
   Size = SIZE(&Grid);
   if( Grid.Dim[0]<2 || Grid.Dim[1]<2 || Grid.Dim[2]<2 )
      mexErrMsgTxt("raytrace.c: Grid dimension(s) too small");
   if( (Tmp = mxGetField(prhs[0],0,"Y")) == NULL || !mxIsDouble(Tmp))
      mexErrMsgTxt("raytrace.c: No Grid.Y field or wrong type");
   if( mxGetNumberOfDimensions(Tmp) != 3 )
      mexErrMsgTxt("raytrace.c: Grid.Y must be 3 dimensional");
   pInt = (int *)mxGetDimensions(Tmp);
   Grid.Dim[0] = pInt[0]; Grid.Dim[1] = pInt[1]; Grid.Dim[2] = pInt[2];
   Size = SIZE(&Grid);
   if( Size != SIZE(&Grid) ) mexErrMsgTxt("raytrace.c: X,Y,Z,F must be same size");
   if( (Tmp = mxGetField(prhs[0],0,"Z")) == NULL || !mxIsDouble(Tmp))
      mexErrMsgTxt("raytrace.c: No Grid.Z field or wrong type");
   if( mxGetNumberOfDimensions(Tmp) != 3 )
      mexErrMsgTxt("raytrace.c: Grid.Z must be 3 dimensional");
   pInt = (int *)mxGetDimensions(Tmp);
   Grid.Dim[0] = pInt[0]; Grid.Dim[1] = pInt[1]; Grid.Dim[2] = pInt[2];
   Size = SIZE(&Grid);
   if( Size != SIZE(&Grid) )
      mexErrMsgTxt("raytrace.c: X,Y,Z,F must be same size");
   if( (Tmp = mxGetField(prhs[0],0,"F")) != NULL && Trace.Prop != LINEAR )
   {
      if( !mxIsDouble(Tmp) ) mexErrMsgTxt("raytrace.c: Grid.F wrong type");
      if( mxGetNumberOfDimensions(Tmp) != 3 )
         mexErrMsgTxt("Grid.F must be 3 dimensional");
      pInt = (int *)mxGetDimensions(Tmp);
      Grid.Dim[0] = pInt[0]; Grid.Dim[1] = pInt[1]; Grid.Dim[2] = pInt[2];
      Size = SIZE(&Grid);
      if( Size != SIZE(&Grid) )
         mexErrMsgTxt("raytrace.c: X,Y,Z,F must be same size");
      Grid.F = mxGetPr(mxGetField(prhs[0],0,"F"));
   }
   else
   {
      if( Trace.Prop != LINEAR ) mexErrMsgTxt("raytrace.c: No Grid.F field");
      Grid.F = NULL;
   }
   Grid.X = mxGetPr(mxGetField(prhs[0],0,"X"));
   Grid.Y = mxGetPr(mxGetField(prhs[0],0,"Y"));
   Grid.Z = mxGetPr(mxGetField(prhs[0],0,"Z"));

   if( (Path = (PATH *)malloc(sizeof(PATH)))==NULL )
      mexErrMsgTxt("Memory allocation failure");

   if( nlhs == 2 && (Sect = (SECT *)malloc(sizeof(SECT)))==NULL )
      mexErrMsgTxt("Memory allocation failure");

/* Solve the problem */
   if( (Err = RayTrace( &Grid, &Trace, Path, Sect )) != OK )
   {
      sprintf(Str,"raytrace.c: Failed with code %d", Err );
      mexErrMsgTxt( Str );
   }

/* Create MatLAB arrays for the output */
   pTmp[0] = mxCreateDoubleMatrix( 1, Path->i, mxREAL );
   pTmp[1] = mxCreateDoubleMatrix( 1, Path->i, mxREAL );
   pTmp[2] = mxCreateDoubleMatrix( 1, Path->i, mxREAL );
   pTmp[3] = mxCreateDoubleScalar( Path->rows ); /* was  pTmp[3] = mxCreateScalarDouble( Path->rows ); */
   pTmp[4] = mxCreateDoubleScalar( Path->cols ); /* was  pTmp[4] = mxCreateScalarDouble( Path->cols ); */

   memcpy( mxGetPr(pTmp[0]),Path->row,Path->i*sizeof(double) );
   memcpy( mxGetPr(pTmp[1]),Path->col,Path->i*sizeof(double) );
   memcpy( mxGetPr(pTmp[2]),Path->val,Path->i*sizeof(double) );

   mexCallMATLAB(1,plhs,5,pTmp,"sparse");

   if( nlhs == 2 && Sect != NULL )
   {
      plhs[1] = mxCreateDoubleMatrix( 4, Sect->i/4, mxREAL );
      memcpy( mxGetPr(plhs[1]), Sect->val, Sect->i*sizeof(double) );
   }

   mxDestroyArray(pTmp[0]);
   mxDestroyArray(pTmp[1]);
   mxDestroyArray(pTmp[2]);
   mxDestroyArray(pTmp[3]);
   mxDestroyArray(pTmp[4]);

   if( Path != NULL && Path->row != NULL ) free(Path->row);
   if( Path != NULL && Path->col != NULL ) free(Path->col);
   if( Path != NULL && Path->val != NULL ) free(Path->val);
   if( Sect != NULL && Sect->val != NULL ) free(Sect->val);
   return;
}




