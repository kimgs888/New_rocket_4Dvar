/**************************************************************************************
File:    trace.c
Author:  Paul S.J. Spencer
Purpose: Ray tracing through a warped cartesian vertex grid
**************************************************************************************/

#include "trace.h"

/* Include for debugging */
 #define VERBOSE

#ifdef VERBOSE
   extern int mexPrintf(const char *format, ...);
#else
   int mexPrintf(const char *format, ...) { return 0;}
#endif


/* Tetrahedral decomposition of a cube, provides indices to vertices in voxel.
   Faces of index 0 and 1 are external to voxel, faces 2 and 3 are internal
   Tetra[TetraIndex][FaceIndex][VertexIndex] */
UINT Tetra[6][4][3] = {
                         { {6,7,5},{5,6,1},{6,7,1},{5,7,1} },
                         { {2,1,6},{6,7,2},{2,1,7},{6,7,1} },
                         { {2,3,1},{3,7,2},{3,7,1},{2,1,7} },
                         { {4,7,5},{5,4,1},{5,7,1},{4,7,1} },
                         { {0,4,1},{4,7,0},{4,7,1},{0,7,1} },
                         { {0,3,1},{3,0,7},{3,7,1},{0,7,1} },
                      };


/* Vertex origin for gradient */
UINT TetraOrigin[6] = {1,7,7,1,7,7};

/* Index of adjoining tetrahedral to a particular face, Face < 2 implies external
   Face[TetraIndex][FaceIndex] */
UINT AdjoinTetra[6][4] = {{2,4,1,3},{5,3,2,0},{0,4,5,1},{5,1,0,4},{2,0,3,5},{3,1,2,4}};

/* Index of adjoining face to a particular face, Face < 2 implies external
   Face[TetraIndex][FaceIndex] */
UINT AdjoinFace[6][4]  = {{0,1,3,2},{1,1,3,2},{0,0,2,2},{0,1,3,2},{1,1,3,3},{0,0,2,3}};




/**************************************************************************************
Function: RayTrace
Purpose:  Entry point function for ray tracing multiple rays
Returns:  Non-zero on error
Notes:    The matrices Trace.Path and Trace.Sect must be free'd externally.  If Sect
          is NULL no terms are created.
**************************************************************************************/
ERRORS RayTrace
(
   GRID    *Grid,  /* Grid definition structure                            */
   TRACE   *Trace, /* Ray trace parameter structure                        */
   PATH    *Path,  /* Output sparse optical path matrix, H                 */
   SECT    *Sect   /* Output intersection matrix, [RayIndex..;x..;y..;z..] */
)
{
   SINT i, j, k;
   long Vx[3];
   double de;
   RAY  Ray;
   UINT Err;
   Ray.Rm = 0;

/* Initialise Path Matrix */
   if( Path != NULL )
   {
      Path->size = MAX(Trace->Tx.Rows,Trace->Tv.Rows) * MAXPATH * 4;
      Path->rows = 0; /* Use as ray counter in shoot */ 
      Path->cols = SIZE(Grid);
      Path->i    = 0;

      if( Path->size*sizeof(double) > MAXMEM ) Path->size = MAXMEM/sizeof(double);  

      if( (Path->row=(double *)malloc(sizeof(double)*Path->size))==NULL ) return MEMORY;
      if( (Path->col=(double *)malloc(sizeof(double)*Path->size))==NULL ) return MEMORY;
      if( (Path->val=(double *)malloc(sizeof(double)*Path->size))==NULL ) return MEMORY;
   }

/* Initialise intersection matrix */
   if( Sect != NULL )
   {
      Sect->size = MAX(Trace->Tx.Rows,Trace->Tv.Rows) * MAXPATH * 4;
      Sect->i    = 0;

      if( Sect->size*sizeof(double) > MAXMEM ) Sect->size = MAXMEM/sizeof(double);  
      if( (Sect->val=(double *)malloc(sizeof(double)*Sect->size))==NULL ) return MEMORY;
   }

/* Define wrapping condition for third dimension */
   if( (pow(Grid->X[Vox(Grid,0,0,0)]-Grid->X[Vox(Grid,0,0,Grid->Dim[2]-1)],2)+
        pow(Grid->Y[Vox(Grid,0,0,0)]-Grid->Y[Vox(Grid,0,0,Grid->Dim[2]-1)],2)+
        pow(Grid->Z[Vox(Grid,0,0,0)]-Grid->Z[Vox(Grid,0,0,Grid->Dim[2]-1)],2)) <
     2*(pow(Grid->X[Vox(Grid,0,0,0)]-Grid->X[Vox(Grid,0,0,1)],2)+
        pow(Grid->Y[Vox(Grid,0,0,0)]-Grid->Y[Vox(Grid,0,0,1)],2)+
        pow(Grid->Z[Vox(Grid,0,0,0)]-Grid->Z[Vox(Grid,0,0,1)],2)) && Grid->Dim[2]>2 )
        Grid->Wrap = 1;
   else Grid->Wrap = 0;


/* Offset to prevent propagation along boundary */
   de = (fabs(Grid->X[1]-Grid->X[0])+fabs(Grid->Y[1]-Grid->Y[0]))*1E-8;
   for( i = 0; i < MAX(Trace->Tx.Rows,Trace->Rx.Rows); i++ )
   {
      for( j = 0; j < 3; j++ )
      {
         if( !IsEmpty(Trace->Tx) && i < Trace->Tx.Rows ) 
            Elem(Trace->Tx,i,j)=Elem(Trace->Tx,i,j)-de*(j+1);
         if( !IsEmpty(Trace->Rx) && i < Trace->Rx.Rows ) 
            Elem(Trace->Rx,i,j)=Elem(Trace->Rx,i,j)+de*(j+1);
      }
   }

/* Trace all rays */
   for( i = 0; i < MAX(Trace->Tx.Rows,Trace->Tv.Rows); i++ )
   {
      double m;

      Ray.A     = Elem(Trace->A,Trace->A.Rows==1?0:i,0);
      Ray.Prop  = Trace->Prop;
      Ray.Ix[0] = (long)Elem(Trace->Ix,Trace->Ix.Rows==1?0:i,0);
      Ray.Ix[1] = (long)Elem(Trace->Ix,Trace->Ix.Rows==1?0:i,1);
      Ray.Ix[2] = (long)Elem(Trace->Ix,Trace->Ix.Rows==1?0:i,2);

      Ray.p.x   = Elem(Trace->Tx,Trace->Tx.Rows==1?0:i,0);
      Ray.p.y   = Elem(Trace->Tx,Trace->Tx.Rows==1?0:i,1);
      Ray.p.z   = Elem(Trace->Tx,Trace->Tx.Rows==1?0:i,2);
      Ray.v.x   = Elem(Trace->Tv,Trace->Tv.Rows==1?0:i,0);
      Ray.v.y   = Elem(Trace->Tv,Trace->Tv.Rows==1?0:i,1);
      Ray.v.z   = Elem(Trace->Tv,Trace->Tv.Rows==1?0:i,2);
      m         = MAG(Ray.v); Ray.v.x/=m; Ray.v.y/=m; Ray.v.z/=m;
      
      if( !IsEmpty(Trace->Rx) )
      {
         Ray.Rx.x = Elem(Trace->Rx,Trace->Rx.Rows==1?0:i,0);
         Ray.Rx.y = Elem(Trace->Rx,Trace->Rx.Rows==1?0:i,1);
         Ray.Rx.z = Elem(Trace->Rx,Trace->Rx.Rows==1?0:i,2);
         Ray.Rv.x = Elem(Trace->Rv,Trace->Rv.Rows==1?0:i,0);
         Ray.Rv.y = Elem(Trace->Rv,Trace->Rv.Rows==1?0:i,1);
         Ray.Rv.z = Elem(Trace->Rv,Trace->Rv.Rows==1?0:i,2);
         Ray.Rm   = MAG(Ray.Rv); Ray.Rv.x/=Ray.Rm; Ray.Rv.y/=Ray.Rm; Ray.Rv.z/=Ray.Rm;
      }

      if( (Err=Shoot( Grid, &Ray, Path, Sect )) != OK && Err != PLANE_SECT ) break;
      Path->rows++; 
   }

   return Err;
}





/**************************************************************************************
Function: Shoot
Purpose:  Shoot a ray, returns when ray exits grid
Returns:  Non-zero on error
Notes:    The entry point of the ray may be within a voxel.
          Ensure that the Path and Sect matrices are initialised to empty with
          Sect.Rows = 4
**************************************************************************************/
ERRORS Shoot
(
   GRID    *Grid,   /* Grid definition                                            */
   RAY     *Ray,    /* Entering ray returned as exit ray                          */
   PATH    *Path,   /* Output appended optical path matrix with ray data, or NULL */
   SECT    *Sect    /* Append intersection matrix, or NULL                        */
)
{
   VOXEL  Voxel;         /* Current voxel                       */
   UINT   Tx;            /* Current tetrahedral index           */
   UINT   Fx = 4;        /* Current face index                  */
   double Dx[2]={.5,.5}; /* Parametric coords of entry/exit     */
   long   Vx[3];         /* Voxel index                         */
   long   i,j,k;         /* Looping variables                   */
   int    Count=0;       /* Voxel counter                       */
   ERRORS Err;           /* Error state                         */

   int SectIo=Sect!=NULL?Sect->i:0; /* Starting Sect point      */
   int PathIo=Path!=NULL?Path->i:0; /* Starting Path point      */

/* Compute ray entry point */
   Vx[0]=Ray->Ix[0]; Vx[1]=Ray->Ix[1]; Vx[2]=Ray->Ix[2];

   if( (Err=SetVoxel(Grid,Vx,&Voxel))!=OK || (Tx=InTetra(&Voxel,Ray)) > 5 )
   {
      for( i = -2; i <= 2; i++ )
      {
         for( j = -2; j <= 2; j++ )
         {
            for( k = -2; k <= 2; k++ )
            {
               Vx[0]=Ray->Ix[0]+i; Vx[1]=Ray->Ix[1]+j; Vx[2]=Ray->Ix[2]+k;
               if( (Err=SetVoxel(Grid,Vx,&Voxel))==OK && (Tx=InTetra(&Voxel,Ray)) <= 5 )
                  break;
            }
            if( Tx <= 5 ) break;
         }
         if( Tx <= 5 ) break;
      }
   }
   if( Err != OK || Tx > 5 ) Err = ENTRY_RAY;

   if( Sect != NULL ) /* Starting point */
   {
      if( Sect->i+4>=Sect->size ) Err = DOVERFLOW;
      else
      {
         Sect->val[Sect->i++] = Path->rows;
         Sect->val[Sect->i++] = Ray->p.x;
         Sect->val[Sect->i++] = Ray->p.y;
         Sect->val[Sect->i++] = Ray->p.z;
      }
   }

   while( Err==OK && (Err = Trace( &Voxel, Ray, &Tx, &Fx, Dx ))==OK && Count++<MAXPATH )
   {
      if( Sect != NULL ) /* Current ray position */
      {
         if( Sect->i+4>=Sect->size ) {Err = DOVERFLOW; break; };
         Sect->val[Sect->i++] = Path->rows;
         Sect->val[Sect->i++] = Ray->p.x;
         Sect->val[Sect->i++] = Ray->p.y;
         Sect->val[Sect->i++] = Ray->p.z;
      }

      if( Path != NULL && Err == OK )
      {
         for( i = 0; i < 8; i++ ) /* Find the voxel index for each vertex */
         {
            if( Ray->Path[i] != 0.0 )
            {
               int Vi;
               switch( i )
               {
                  case 0: { Vi = Vox(Grid,Vx[0]  ,Vx[1],  Vx[2]);   }; break;
                  case 1: { Vi = Vox(Grid,Vx[0]+1,Vx[1],  Vx[2]);   }; break;
                  case 2: { Vi = Vox(Grid,Vx[0]+1,Vx[1]+1,Vx[2]);   }; break;
                  case 3: { Vi = Vox(Grid,Vx[0]  ,Vx[1]+1,Vx[2]);   }; break;
                  case 4: { Vi = Vox(Grid,Vx[0]  ,Vx[1],  Vx[2]+1); }; break;
                  case 5: { Vi = Vox(Grid,Vx[0]+1,Vx[1],  Vx[2]+1); }; break;
                  case 6: { Vi = Vox(Grid,Vx[0]+1,Vx[1]+1,Vx[2]+1); }; break;
                  case 7: { Vi = Vox(Grid,Vx[0],  Vx[1]+1,Vx[2]+1); }; break;
               }

               if( Path->i >= Path->size ) {Err = DOVERFLOW; break;} 
               Path->row[Path->i]   = Path->rows+1;
               Path->col[Path->i]   = Vi+1;
               Path->val[Path->i++] = Ray->Path[i];
            }
         }
      }
      if( Fx == 4 || Err != OK ) break; 

/* Compute new entry point voxel, Fx and Tx are already correct */
      if(      Tx==0 && Fx==0 ) Vx[2] -= 1;
      else if( Tx==0 && Fx==1 ) Vx[0] -= 1;
      else if( Tx==1 && Fx==0 ) Vx[0] -= 1;
      else if( Tx==1 && Fx==1 ) Vx[1] -= 1;
      else if( Tx==2 && Fx==0 ) Vx[2] += 1;
      else if( Tx==2 && Fx==1 ) Vx[1] -= 1;
      else if( Tx==3 && Fx==0 ) Vx[2] -= 1;
      else if( Tx==3 && Fx==1 ) Vx[1] += 1;
      else if( Tx==4 && Fx==0 ) Vx[1] += 1;
      else if( Tx==4 && Fx==1 ) Vx[0] += 1;
      else if( Tx==5 && Fx==0 ) Vx[2] += 1;
      else if( Tx==5 && Fx==1 ) Vx[0] += 1;

      if( Err == OK ) Err = SetVoxel( Grid, Vx, &Voxel );
   }

   if( Ray->Rm != 0 && Fx != 4 ) Err = PLANE_SECT; 
   if( Err == BOUND ) Err = OK;
   

/* On error remove Path and Sect data */
   if( Sect != NULL && (Err != OK || Count>=MAXPATH) ) Sect->i = SectIo;
   if( Path != NULL && (Err != OK || Count>=MAXPATH) ) Path->i = PathIo;

/*   if( Err != OK ) mexPrintf("Error %d\n",Err); */

   return OK;
}

/**************************************************************************************
Function: Vox
Purpose:  Vertex index from component indices in each dimension
Returns:  Vertex index
Notes:    This function will always return a valid voxel index within range
**************************************************************************************/
long Vox
(
   GRID *Grid,   /* Grid definition structure */
   long  i,
   long  j,
   long  k
)
{
   return (long)((i % Grid->Dim[0])+(i<0?Grid->Dim[0]:0)) +
                ((j % Grid->Dim[1])+(j<0?Grid->Dim[1]:0))*Grid->Dim[0] +
                ((k % Grid->Dim[2])+(k<0?Grid->Dim[2]:0))*Grid->Dim[0]*Grid->Dim[1];
}

/**************************************************************************************
Function: InGrid
Purpose:  Test if a vertex is within the grid
Returns:  Booolean
Notes:    To test a voxel referenced by its lower left corner test with
**************************************************************************************/
int InGrid
(
   GRID *Grid,   /* Grid definition structure */
   long  i,
   long  j,
   long  k
)
{
   return (i>=0&&i<Grid->Dim[0])&&(j>=0&&j<Grid->Dim[1])&&
          (Grid->Wrap==1?1:(k>=0&&k<Grid->Dim[2]));
}

/**************************************************************************************
Function: SetVoxel
Purpose:  Set voxel coordinates and field from voxel index
Returns:  Non-zero on error
Notes:
**************************************************************************************/
ERRORS SetVoxel
(
   GRID *Grid,   /* Grid definition structure                          */
   long  Vx[3],  /* Voxel index                                        */
   VOXEL *V      /* Output voxel with coordinates and field values set */
)
{
   long  m;      /* Grid voxel index of current vertex */

   if( !InGrid(Grid,Vx[0],Vx[1],Vx[2]) || !InGrid(Grid,Vx[0]+1,Vx[1]+1,Vx[2]+1) )
      return BOUND;

   m = Vox(Grid,Vx[0],Vx[1],Vx[2]);
   if( m < 0 || m >= SIZE(Grid) ) return BOUND;

   if( Grid->F != NULL ) V->Field[0] = Grid->F[m]; 
   V->Vertex[0].x=Grid->X[m]; V->Vertex[0].y=Grid->Y[m]; V->Vertex[0].z=Grid->Z[m];
   m = Vox(Grid,Vx[0]+1,Vx[1],Vx[2]);     if( Grid->F != NULL ) V->Field[1] = Grid->F[m];
   V->Vertex[1].x=Grid->X[m]; V->Vertex[1].y=Grid->Y[m]; V->Vertex[1].z=Grid->Z[m];
   m = Vox(Grid,Vx[0]+1,Vx[1]+1,Vx[2]);   if( Grid->F != NULL ) V->Field[2] = Grid->F[m];
   V->Vertex[2].x=Grid->X[m]; V->Vertex[2].y=Grid->Y[m]; V->Vertex[2].z=Grid->Z[m];
   m = Vox(Grid,Vx[0],Vx[1]+1,Vx[2]);     if( Grid->F != NULL ) V->Field[3] = Grid->F[m];
   V->Vertex[3].x=Grid->X[m]; V->Vertex[3].y=Grid->Y[m]; V->Vertex[3].z=Grid->Z[m];
   m = Vox(Grid,Vx[0],Vx[1],Vx[2]+1);     if( Grid->F != NULL ) V->Field[4] = Grid->F[m];
   V->Vertex[4].x=Grid->X[m]; V->Vertex[4].y=Grid->Y[m]; V->Vertex[4].z=Grid->Z[m];
   m = Vox(Grid,Vx[0]+1,Vx[1],Vx[2]+1);   if( Grid->F != NULL ) V->Field[5] = Grid->F[m];
   V->Vertex[5].x=Grid->X[m]; V->Vertex[5].y=Grid->Y[m]; V->Vertex[5].z=Grid->Z[m];
   m = Vox(Grid,Vx[0]+1,Vx[1]+1,Vx[2]+1); if( Grid->F != NULL ) V->Field[6] = Grid->F[m];
   V->Vertex[6].x=Grid->X[m]; V->Vertex[6].y=Grid->Y[m]; V->Vertex[6].z=Grid->Z[m];
   m = Vox(Grid,Vx[0],Vx[1]+1,Vx[2]+1);   if( Grid->F != NULL ) V->Field[7] = Grid->F[m];
   V->Vertex[7].x=Grid->X[m]; V->Vertex[7].y=Grid->Y[m]; V->Vertex[7].z=Grid->Z[m];

   return OK;
}



/**************************************************************************************
Function: InTetra
Purpose:  Determine which tetrahedral a ray is located within
Returns:  Index of tetrahedral, 99 on error
Notes:
**************************************************************************************/
UINT InTetra
(
   VOXEL  *Voxel,     /* Voxel vertices ordered as in diagram in trace.h */
   RAY    *Ray        /* Ray definition, uses Ray.p,Ray.v                */
)
{
   double Dx[2];
   UINT   i,j,k;
   CPOINT P;

   for( i = 0; i < 6; i++ )        /* Loop through tetrahedra */
   {
      k = 0;
      for( j = 0; j < 4; j++ )     /* Loop through all faces */
      {
         if( LineTriangle( &Ray->p, &Ray->v, &Voxel->Vertex[Tetra[i][j][0]],
             &Voxel->Vertex[Tetra[i][j][1]],
             &Voxel->Vertex[Tetra[i][j][2]],Dx,&P)>0.0 ) k++;
      }
      if( k == 1 ) return i;
   }
   return 99;
}


/**************************************************************************************
Function: Trace
Purpose:  Compute ray path through a cartesian voxel
Returns:  Non-zero on error
Notes:    Path.NumSect and Path.NumPath should be initialised externally
**************************************************************************************/
ERRORS Trace
(
   VOXEL  *Voxel,     /* Voxel vertices ordered as in diagram above */
   RAY    *Ray,       /* Ray definition                             */
   UINT   *Tx,        /* Entry point tetrahedral, returned updated  */
   UINT   *Face,      /* Entry point face index, returned updated   */
   double *De         /* Parametic coords of entry                  */
)
{
   VECTOR t, r, s;         /* Gradient in field t, propagation plane vectors s,n */
   double Dx[2];           /* Parametric coords of exit                          */
   int    i;               /* Looping variable                                   */
   int    Loops = 0;       /* Trapped ray loop counter                           */

   for(i=0;i<8;i++) Ray->Path[i] = 0;

   do /* While within this voxel */
   {
      double OPath = LARGE;   /* Current minimum optical path within tetrahedron */
      double GPath;           /* Group path                                      */
      UINT   OFace;           /* Face of outgoing ray                            */
      UINT   F;               /* Current face                                    */
      double fd;              /* Fractional distance along ray                   */
      double tm, rm;          /* Magnitude of field gradient and n (|n|=sin(i))  */
      double e;               /* Interpolated field value at ray entry point     */
      CPOINT Ex, EX;          /* Ray exit point                                  */
      VECTOR Vx, VX;          /* Ray exit vector                                 */
      CPOINT O;               /* Origin of s,t coord system in cartesians        */
      CPOINT I[2];            /* Cartesian face intersection points              */
      double sc[2],tc[2];     /* Parametric face s,t intersection points         */
      double Root[2];         /* Intersection roots                              */
      int    NRoot;           /* Number of roots                                 */

/* Compute normalised propagation vectors s,t,r */
      if( Ray->Prop != LINEAR )
      {
         int O = Gradient( Voxel, *Tx, &t );

         e = Voxel->Field[O] + (Ray->p.x - Voxel->Vertex[O].x)*t.x +
                               (Ray->p.y - Voxel->Vertex[O].y)*t.y +
                               (Ray->p.z - Voxel->Vertex[O].z)*t.z;

         if( (tm = MAG(t)) > 1E-11*Ray->A ) {t.x=t.x/tm; t.y=t.y/tm; t.z=t.z/tm;} /* was 1E-8 */
         else tm = 0;

         r.x = CROSSX(Ray->v,t); r.y = CROSSY(Ray->v,t); r.z = CROSSZ(Ray->v,t);
         if( (rm=MAG(r)) > 1E-11 ) {r.x=r.x/rm; r.y=r.y/rm; r.z=r.z/rm;} /* was 1E-9 */
         else
         {
            r.x=t.y; r.y=t.z; r.z=t.x; rm = 0;
            r.x = CROSSX(r,Ray->v); r.y = CROSSY(r,Ray->v); r.z = CROSSZ(r,Ray->v);
            fd=MAG(r); r.x=r.x/fd; r.y=r.y/fd; r.z=r.z/fd;
         }
         s.x = CROSSX(t,r); s.y = CROSSY(t,r); s.z = CROSSZ(t,r);

         if( Ray->Prop==PARABOLIC || Ray->Prop==PARABOLIC_GROUP )
         {
            if( Ray->A*e >= 1.0-1E-6 ) return 9991;                /* No path      */
            if( Ray->A*e < 1E-6 ) tm = 0; 
            t.x=-t.x; t.y=-t.y; t.z=-t.z;                          /* Reverse t    */
         }
         if( Ray->Prop==CIRCULAR )
         {
            if( e <= 0.0 ) return 9992;                            /* No path      */
         }
      }

/* -------------------------------- Linear Rays ----------------------------------------*/
      if( Ray->Prop == LINEAR || tm == 0 )
      {
/* mexPrintf("LINEAR "); */
         VX.x = Ray->v.x; VX.y = Ray->v.y; VX.z = Ray->v.z;

         for( OFace = 0; OFace < 4; OFace++ )
         {
            if( OFace != *Face && ((fd = LineTriangle( &Ray->p, &Ray->v,
                &Voxel->Vertex[Tetra[*Tx][OFace][0]], &Voxel->Vertex[Tetra[*Tx][OFace][1]],
                &Voxel->Vertex[Tetra[*Tx][OFace][2]], Dx, &EX )) > 0 )  )
            {
               if( Ray->Prop == PARABOLIC || Ray->Prop == PARABOLIC_GROUP )
               { double ne=sqrt(1-Ray->A*e); OPath = fd*ne; GPath = fd/ne; }
               else if( Ray->Prop == CIRCULAR )
               { double ne=Ray->A/e; OPath = fd*ne; GPath = fd/ne; }
               else OPath = fd;
               break;
            }
         }

         if( OFace == 4 ) {return TRACE_FAIL;}

         if( Ray->Rm != 0 ) /* Test for Rx plane intersection */
         {
            double d;
            CPOINT C;

            C.x = Ray->p.x+Ray->v.x; C.y = Ray->p.y+Ray->v.y; C.z = Ray->p.z+Ray->v.z;

            if( (d = LinePlane(&Ray->p,&C,&Ray->Rv,&Ray->Rx)) > 0.0 && d < fd )
            {
               if( Ray->Prop == PARABOLIC || Ray->Prop == PARABOLIC_GROUP )
               { double ne=sqrt(1-Ray->A*e); OPath = d*ne; GPath = d/ne; }
               else if( Ray->Prop == CIRCULAR )
               { double ne=Ray->A/e; OPath = d*ne; GPath = fd/ne; }
               else OPath = d;
               EX.x = Ray->p.x+d*VX.x; EX.y = Ray->p.y+d*VX.y; EX.z = Ray->p.z+d*VX.z;
               Dx[0]=0.5; Dx[1]=0.5; OFace = 4;
            }
         }
      }
/* -------------------------------- Parabolic Rays -------------------------------------*/
      else if( Ray->Prop == PARABOLIC || Ray->Prop == PARABOLIC_GROUP )
      {
         double nn = 1.0-Ray->A*e;               /* RI^2 at entry ray point     */
         double sini = rm;                       /* sin(i)                      */
         double cosi = DOT(t,Ray->v);            /* cos(i)                      */
         double at = tm==0?0:sqrt(Ray->A*tm);
         double kk = nn*sini*sini;
         double ep = 4*kk/(Ray->A*tm);           /* Parabola constant, s^2=ep*t */
         double es = ep>0?sqrt(ep):0;            /* sqrt(ep)                    */
         double to = (nn-kk)/(Ray->A*tm);        /* Ray entry in t              */
         double so = sini==0?0:0.5*ep*cosi/sini; /* Ray entry in s              */
         double st = to<=0?0:SIGN(cosi)*sqrt(to);

         O.x = Ray->p.x-t.x*to-s.x*so;           /* Parabola origin             */
         O.y = Ray->p.y-t.y*to-s.y*so;
         O.z = Ray->p.z-t.z*to-s.z*so;

         for( F = 0; F < 4+(UINT)(Ray->Rm!=0); F++ ) /* Test faces for intersections */
         {
            if( (F<4  && TrianglePlane( &Voxel->Vertex[Tetra[*Tx][F][0]],
                                        &Voxel->Vertex[Tetra[*Tx][F][1]],
                                        &Voxel->Vertex[Tetra[*Tx][F][2]],
                                        &r, &Ray->p, &I[0], &I[1] ) == OK) ||
                (F==4 && PlanePlane(&r,&Ray->Rv,&Ray->p,&Ray->Rx,&I[0],&I[1]) == OK) )
            {
               double ds,dt;  /* Slope of intersecting line = dt/ds */
               double OP[2];  /* Optical path                       */
               double m,c,tx,sx;

               sc[0] = (I[0].x-O.x)*s.x+(I[0].y-O.y)*s.y+(I[0].z-O.z)*s.z;
               tc[0] = (I[0].x-O.x)*t.x+(I[0].y-O.y)*t.y+(I[0].z-O.z)*t.z;
               sc[1] = (I[1].x-O.x)*s.x+(I[1].y-O.y)*s.y+(I[1].z-O.z)*s.z;
               tc[1] = (I[1].x-O.x)*t.x+(I[1].y-O.y)*t.y+(I[1].z-O.z)*t.z;
               ds    = sc[1]-sc[0];
               dt    = tc[1]-tc[0];
               m     = ds==0.0 ? 0.0 : dt/ds;
               c     = tc[0] - m*sc[0];

/* mexPrintf("Face %d ds %.10f kk %f Root[0] %f to %f\n",F,ds,kk,pow(Root[0],2),to); */
               if( fabs(ds) < 1E-9 ) /* Single root, problem could be here */
               {
                  if( ep == 0 || (F==*Face & F!=4) ) NRoot = 0;
                  else if( sini > 0.5 )
                  {
                     Root[0] = sc[0]/es; NRoot=1;
                     OP[0] = sqrt(kk)*( sc[0]-so + 4.0/(3.0*ep*ep)*(pow(sc[0],3)-pow(so,3)) );
                  } 
                  else
                  {
                     Root[0] = sc[0]/es; NRoot=1;
                     OP[0]   = at*.5*( (ep*Root[0]+4.0/3.0*pow(Root[0],3))-st*(ep+4.0/3.0*to) ); 
                  }
               }
               else
               {
                  Root[0] = m*m*ep*.25+c; NRoot = Root[0]>0?2:0; 
                  Root[0] = Root[0]>0?sqrt(Root[0]):0;
                  Root[1] = m*es*.5 + Root[0];
                  Root[0] = m*es*.5 - Root[0];
                  OP[0]   = at*.5*( (ep*Root[0]+4.0/3.0*pow(Root[0],3)) - st*(ep+4.0/3.0*to) );
                  OP[1]   = at*.5*( (ep*Root[1]+4.0/3.0*pow(Root[1],3)) - st*(ep+4.0/3.0*to) );
               }

               if( F==*Face && F!=4 && NRoot==2 )
               { if( fabs(OP[1])>fabs(OP[0]) ) {Root[0]=Root[1];OP[0]=OP[1];} NRoot=1;};

               for(i = 0; i < NRoot; i++ )
               {
                  double tx = Root[i]*Root[i];         /* s,t coords at exit */
                  double sx = es*Root[i];
                  double Valid = fabs(ds)>fabs(dt)?(sx-sc[0])/ds:(dt==0?2:(tx-tc[0])/dt);

                  if( OP[i]>0 && OP[i]<OPath && Valid>=0 && Valid<=1 )  /* Test optical path */
                  {
                     EX.x  = O.x+sx*s.x+tx*t.x;
                     EX.y  = O.y+sx*s.y+tx*t.y;
                     EX.z  = O.z+sx*s.z+tx*t.z;
                     VX.x  = es*s.x+2*Root[i]*t.x;
                     VX.y  = es*s.y+2*Root[i]*t.y;
                     VX.z  = es*s.z+2*Root[i]*t.z;
                     OPath = OP[i]; OFace = F;
                     GPath = 2/at*(Root[i]-st);
                  }
               }
            }
         }
      }
/* -------------------------------- Circular Rays --------------------------------------*/
      else if( Ray->Prop == CIRCULAR )
      {
         double cosi = DOT(Ray->v,t);      /* cos(i)                      */
         double R  = e/(tm*rm);            /* Ray arc radius of curvature */
/*         double to = e/tm;  */               /* Origin of ray arc on t axis */
/*         double so = -to*cosi/rm;  */        /* Origin of ray arc on s axis */
         double to =  R*rm;                /* Origin of ray arc on t axis */
         double so = -R*sqrt(1.0-rm*rm)*SIGN(cosi); /* Origin of ray arc on s axis */

         double ds, dt;

/* mexPrintf("CIRCULAR %7.2E %.10f ",R,pow(to,2.0)+pow(so,2.0)-pow(R,2.0)); */

         O.x = Ray->p.x-t.x*to-s.x*so;     /* Ray arc origin */
         O.y = Ray->p.y-t.y*to-s.y*so;
         O.z = Ray->p.z-t.z*to-s.z*so;

         for( F = 0; F < 4+(UINT)(Ray->Rm!=0); F++ )
         {
            if( (F<4  && TrianglePlane( &Voxel->Vertex[Tetra[*Tx][F][0]],
                                        &Voxel->Vertex[Tetra[*Tx][F][1]],
                                        &Voxel->Vertex[Tetra[*Tx][F][2]],
                                        &r, &Ray->p, &I[0], &I[1] ) == OK) ||
                (F==4 && PlanePlane(&r,&Ray->Rv,&Ray->p,&Ray->Rx,&I[0],&I[1]) == OK) )
            {
               Parametric( &I[0], &O, &s, &t, &sc[0], &tc[0] );
               Parametric( &I[1], &O, &s, &t, &sc[1], &tc[1] );

               ds = sc[1]-sc[0]; dt = tc[1]-tc[0];

               NRoot = Roots( pow(ds,2.0)+pow(dt,2.0), sc[0]*ds+tc[0]*dt,
                              sc[0]*sc[0]+tc[0]*tc[0]-R*R, Root );

               if( F==*Face && NRoot==1 ) NRoot = 0;
               if( F==*Face && NRoot==2 ) /* Test for reflecting ray */
               {
                  if( fabs(tc[0] + Root[0]*dt-to)+fabs(sc[0] + Root[0]*ds-so) <
                      fabs(tc[0] + Root[1]*dt-to)+fabs(sc[0] + Root[1]*ds-so) )
                     Root[0] = Root[1];
                  NRoot = 1;
               }

               for( i = 0; i < NRoot; i++ ) /* Loop for each root */
               {
                  VECTOR Rx;                /* Vector from origin to ray exit point */
                  double Rm, p;
                  double tx = tc[0] + Root[i]*dt; /* s,t coords at exit */
                  double sx = sc[0] + Root[i]*ds;

                  Ex.x = I[0].x+Root[i]*(I[1].x-I[0].x);
                  Ex.y = I[0].y+Root[i]*(I[1].y-I[0].y);
                  Ex.z = I[0].z+Root[i]*(I[1].z-I[0].z);

                  Rx.x = Ex.x-O.x; Rx.y = Ex.y-O.y; Rx.z = Ex.z-O.z;
                  Rm   = MAG(Rx);  Rx.x/=Rm; Rx.y/=Rm; Rx.z/=Rm;
                  Vx.x = CROSSX(Rx,r); Vx.y = CROSSY(Rx,r); Vx.z = CROSSZ(Rx,r);


                  p = (R-so)*tx/(to*(R-sx));
                  p = p>0?Ray->A/tm*log(p):LARGE;

                  if( p < OPath && p > 0.0 )            /* Min optical path? */
                  {
                     OPath = p; OFace = F;
                     EX.x = Ex.x; EX.y = Ex.y; EX.z = Ex.z;
                     VX.x = Vx.x; VX.y = Vx.y; VX.z = Vx.z;
                  }
               }
            }
         }
      }
      else return 99;
/* -------------------------------------------------------------------------------------*/
/* mexPrintf("OFace %d OPath %.4f\n",OFace,OPath); */

      if( OPath == LARGE ) return(ERRORS)999;

/* Phase path contributions per vertex in tetrahedron */
      if( Ray->Prop == PARABOLIC_GROUP ) OPath = GPath;

      ParametricB(&EX,&Voxel->Vertex[Tetra[*Tx][OFace%4][0]],
         &Voxel->Vertex[Tetra[*Tx][OFace%4][1]],&Voxel->Vertex[Tetra[*Tx][OFace%4][2]],
         &Dx[0],&Dx[1]);

      if( *Face == 4 ) {De[0] = Dx[0]; De[1] = Dx[1];}
      else
      {
         ParametricB(&Ray->p,&Voxel->Vertex[Tetra[*Tx][(*Face)%4][0]],
            &Voxel->Vertex[Tetra[*Tx][(*Face)%4][1]],&Voxel->Vertex[Tetra[*Tx][(*Face)%4][2]],
            &De[0],&De[1]);
      }

      Ray->Path[Tetra[*Tx][OFace%4][0]]   += OPath*(1-(Dx[0]+Dx[1]))/2.0; 
      Ray->Path[Tetra[*Tx][OFace%4][1]]   += OPath*Dx[0]/2.0;
      Ray->Path[Tetra[*Tx][OFace%4][2]]   += OPath*Dx[1]/2.0;
      Ray->Path[Tetra[*Tx][(*Face)%4][0]] += OPath*(1-(De[0]+De[1]))/2.0; 
      Ray->Path[Tetra[*Tx][(*Face)%4][1]] += OPath*De[0]/2.0;
      Ray->Path[Tetra[*Tx][(*Face)%4][2]] += OPath*De[1]/2.0;

      tm = MAG(VX);
      Ray->p.x = EX.x; Ray->p.y = EX.y; Ray->p.z = EX.z;    /* New ray position    */
      Ray->v.x=VX.x/tm; Ray->v.y=VX.y/tm; Ray->v.z=VX.z/tm; /* New ray vector      */

      if( Ray->Rm!=0 && OFace==4 )                          /* Intersect Rx plane  */
      {
         *Face=4;
         if( LEN(Ray->p,Ray->Rx) > Ray->Rm ) return PLANE_SECT;
         else return OK;
      }

      *Face = AdjoinFace[*Tx][OFace];                       /* Set new face        */
      *Tx   = AdjoinTetra[*Tx][OFace];                      /* Set new tetrahedron */
      De[0] = Dx[0]; De[1] = Dx[1];                         /* New parametric      */
   }
   while( *Face >= 2 && Loops++ < 12 );                     /* While within voxel  */

   if( Loops == 12 ) return TIMEOUT;

   return OK;
};


/**************************************************************************************
Function: Roots
Purpose:  Compute real roots of a quadratic ax^2+2bx+c=0
Returns:  Number of roots in range 0 to 1
Notes:    Only roots in the range 0 to 1 are returned in Root, a single valid root is
          retured in Root[0]
**************************************************************************************/
int Roots( double a, double b, double c, double Root[2] )
{
   double R = b*b - a*c;
   double *pRoot = Root;

   if( R < 0.0 || b == 0.0 ) return 0;

   if( fabs(a) < 1E-15 )
   {
      *pRoot = -c/(2*b); return (*pRoot>=0.0 && *pRoot<=1.0);
   }

   R = sqrt(R);

   *pRoot = (-b + R)/a; pRoot += (*pRoot>=0.0 && *pRoot<=1.0);
   *pRoot = (-b - R)/a; pRoot += (*pRoot>=0.0 && *pRoot<=1.0);

   return pRoot-Root;
}


/**************************************************************************************
Function: Gradient
Purpose:  Compute gradient within tetrahedron, Go=const, G=[d/dx,d/dy,d/dz]
Returns:  Index of origin vertex
Notes:    Gradient in tetrahedron is computed from an origin at the last vertex of the
          second face in Tetra, ie Tetra[TetIx][1][2].
**************************************************************************************/
int Gradient( VOXEL *Voxel, UINT Tx, VECTOR *G )
{
   double M[3][3], I[3][3];
   int O;

   O = (int)TetraOrigin[Tx];
   M[0][0] = Voxel->Vertex[Tetra[Tx][0][0]].x - Voxel->Vertex[O].x;
   M[1][0] = Voxel->Vertex[Tetra[Tx][0][0]].y - Voxel->Vertex[O].y;
   M[2][0] = Voxel->Vertex[Tetra[Tx][0][0]].z - Voxel->Vertex[O].z;
   M[0][1] = Voxel->Vertex[Tetra[Tx][0][1]].x - Voxel->Vertex[O].x;
   M[1][1] = Voxel->Vertex[Tetra[Tx][0][1]].y - Voxel->Vertex[O].y;
   M[2][1] = Voxel->Vertex[Tetra[Tx][0][1]].z - Voxel->Vertex[O].z;
   M[0][2] = Voxel->Vertex[Tetra[Tx][0][2]].x - Voxel->Vertex[O].x;
   M[1][2] = Voxel->Vertex[Tetra[Tx][0][2]].y - Voxel->Vertex[O].y;
   M[2][2] = Voxel->Vertex[Tetra[Tx][0][2]].z - Voxel->Vertex[O].z;
   if( Inverse3x3( M, I ) != OK ) {G->x=0.0; G->y=0.0; G->z=0.0; return O; }

   G->x = I[0][0]*(Voxel->Field[Tetra[Tx][0][0]]-Voxel->Field[O]) +
          I[1][0]*(Voxel->Field[Tetra[Tx][0][1]]-Voxel->Field[O]) +
          I[2][0]*(Voxel->Field[Tetra[Tx][0][2]]-Voxel->Field[O]);
   G->y = I[0][1]*(Voxel->Field[Tetra[Tx][0][0]]-Voxel->Field[O]) +
          I[1][1]*(Voxel->Field[Tetra[Tx][0][1]]-Voxel->Field[O]) +
          I[2][1]*(Voxel->Field[Tetra[Tx][0][2]]-Voxel->Field[O]);
   G->z = I[0][2]*(Voxel->Field[Tetra[Tx][0][0]]-Voxel->Field[O]) +
          I[1][2]*(Voxel->Field[Tetra[Tx][0][1]]-Voxel->Field[O]) +
          I[2][2]*(Voxel->Field[Tetra[Tx][0][2]]-Voxel->Field[O]);

  return O;
}

/**************************************************************************************
Function: LinePlane
Purpose:  Compute intersection of a line from points A to B with a plane defined by its
          normal n and a point on the plane P
Returns:  Fractional distance from A along AB to intersection, -99 if none
**************************************************************************************/
double LinePlane
(
   CPOINT  *A,      /* Line origin    */
   CPOINT  *B,      /* Line end point */
   VECTOR  *N,      /* Plane normal   */
   CPOINT  *P       /* Point in plane */
)
{
   double Denominator;
   VECTOR AB, PA;

   AB.x = B->x-A->x;  AB.y = B->y-A->y;  AB.z = B->z-A->z;
   PA.x = P->x-A->x;  PA.y = P->y-A->y;  PA.z = P->z-A->z;
   if( (Denominator = DOT(AB,*N)) == 0.0 ) return -99.0;
   return DOT(PA,*N)/Denominator;
}

/**************************************************************************************
Function: TrianglePlane
Purpose:  Compute intersection of a triangle with a plane
Returns:  Zero for valid intersection, non-zero otherwise
**************************************************************************************/
int TrianglePlane
(
   CPOINT  *A,      /* Triangle vertex A         */
   CPOINT  *B,      /* Triangle vertex B         */
   CPOINT  *C,      /* Triangle vertex C         */
   VECTOR  *N,      /* Plane normal              */
   CPOINT  *P,      /* Point in plane            */
   CPOINT  *Ia,     /* First intersection point  */
   CPOINT  *Ib      /* Second intersection point */
)
{
   VECTOR a, b, c, D;
   double sa,sb,sc, Den;
   double F;             /* Fractional intersection */

   a.x = P->x-A->x; a.y = P->y-A->y; a.z = P->z-A->z;
   b.x = P->x-B->x; b.y = P->y-B->y; b.z = P->z-B->z;
   c.x = P->x-C->x; c.y = P->y-C->y; c.z = P->z-C->z;

/* Do sign test to find unique vertex on one side of plane */
   sa = DOT(a,*N); sb = DOT(b,*N); sc = DOT(c,*N);

   if( (sa>=0.0 && sb<0.0 && sc<0.0) || (sa<0.0 && sb>=0.0 && sc>=0.0) )      /* A odd */
   {
      D.x = B->x-A->x; D.y = B->y-A->y; D.z = B->z-A->z; ;
      if( (Den = DOT(D,*N)) == 0.0 ) return 1; F = sa / Den;
      Ia->x = A->x+F*(B->x-A->x); Ia->y = A->y+F*(B->y-A->y); Ia->z = A->z+F*(B->z-A->z);
      D.x = C->x-A->x; D.y = C->y-A->y; D.z = C->z-A->z; ;
      if( (Den = DOT(D,*N)) == 0.0 ) return 1; F = sa / Den;
      Ib->x = A->x+F*(C->x-A->x); Ib->y = A->y+F*(C->y-A->y); Ib->z = A->z+F*(C->z-A->z);
      return 0;
   }
   else if( (sb>=0.0 && sa<0.0 && sc<0.0) || (sb<0.0 && sa>=0.0 && sc>=0.0) ) /* B odd */
   {
      D.x = A->x-B->x; D.y = A->y-B->y; D.z = A->z-B->z; ;
      if( (Den = DOT(D,*N)) == 0.0 ) return 1; F = sb / Den;
      Ia->x = B->x+F*(A->x-B->x); Ia->y = B->y+F*(A->y-B->y); Ia->z = B->z+F*(A->z-B->z);
      D.x = C->x-B->x; D.y = C->y-B->y; D.z = C->z-B->z; ;
      if( (Den = DOT(D,*N)) == 0.0 ) return 1; F = sb / Den;
      Ib->x = B->x+F*(C->x-B->x); Ib->y = B->y+F*(C->y-B->y); Ib->z = B->z+F*(C->z-B->z);
      return 0;
   }
   else if( (sc>=0.0 && sa<0.0 && sb<0.0) || (sc<0.0 && sa>=0.0 && sb>=0.0) ) /* C odd */
   {
      D.x = A->x-C->x; D.y = A->y-C->y; D.z = A->z-C->z; ;
      if( (Den = DOT(D,*N)) == 0.0 ) return 1; F = sc / Den;
      Ia->x = C->x+F*(A->x-C->x); Ia->y = C->y+F*(A->y-C->y); Ia->z = C->z+F*(A->z-C->z);
      D.x = B->x-C->x; D.y = B->y-C->y; D.z = B->z-C->z; ;
      if( (Den = DOT(D,*N)) == 0.0 ) return 1; F = sc / Den;
      Ib->x = C->x+F*(B->x-C->x); Ib->y = C->y+F*(B->y-C->y); Ib->z = C->z+F*(B->z-C->z);
      return 0;
   }

   return 1;
}


/**************************************************************************************
Function: LineTriangle
Purpose:  Determines if a vector intersects a triangle
Returns:  Fractional distance along vector of intersection, -99.0 if no intersection
Notes:    Only intersections in the direction of the vector from its origin are
          considered valid
**************************************************************************************/
double LineTriangle
(
   CPOINT  *O,   /* Origin of intersecting vector                   */
   VECTOR  *V,   /* Intersecting vector                             */
   CPOINT  *Va,  /* Plane origin point                              */
   CPOINT  *Vb,  /* Point in plane defining s axis                  */
   CPOINT  *Vc,  /* Point in plane defining t axis                  */
   double  *Dx,  /* Parametric coords of intersection, Dx[2]        */
   CPOINT  *I    /* Output intersection point of line with triangle */
)
{
   CPOINT B;           /* Line end, A + V                               */
   VECTOR N;           /* Plane normal                                  */
   VECTOR w;           /* Vector from origin (Va) to intersection point */
   VECTOR u;           /* Vector VaVb                                   */
   VECTOR v;           /* Vector VaVc                                   */
   double Den, uv;     /* Denominator and dot(u,v)                      */
   double F;           /* Fractional intersection along AB              */

/* Determine vector end point and the triangle plane axis */
   B.x = O->x + V->x; B.y = O->y + V->y; B.z = O->z + V->z;
   u.x = Vb->x-Va->x; u.y = Vb->y-Va->y; u.z = Vb->z-Va->z;
   v.x = Vc->x-Va->x; v.y = Vc->y-Va->y; v.z = Vc->z-Va->z;
/* Plane normal */
   N.x = CROSSX(u,v); N.y = CROSSY(u,v); N.z = CROSSZ(u,v);

   if( (F = LinePlane(O,&B,&N,Va)) < 0.0 ) {return -99.0;}  /* No intersection */

/* Intersection point on plane */
   I->x = O->x + (B.x-O->x)*F; I->y = O->y + (B.y-O->y)*F; I->z = O->z + (B.z-O->z)*F;
/* Compute vector from origin to intersection point in plane */
   w.x =  I->x-Va->x;  w.y =  I->y-Va->y;  w.z =  I->z-Va->z;

/* Compute denominator */
   uv = DOT(u,v); Den = uv*uv - DOT(u,u)*DOT(v,v);
   if( Den == 0.0 ) return 1;

/* Parametric coords of point in plane */
   Dx[0] = (uv * DOT(w,v) - DOT(v,v)*DOT(w,u))/Den;
   Dx[1] = (uv * DOT(w,u) - DOT(u,u)*DOT(w,v))/Den;

   if( Dx[0]>=0.0 && Dx[1]>=0.0 && Dx[0]+Dx[1]<=1.0 )
   return F;      /* Intersection inside triangle  */
   return -99.0;  /* Intersection outside triangle */
}


/**************************************************************************************
Function: Parametric
Purpose:  Determine the parametric coords of a point in plane
Returns:  Zero if coords are valid, non-zero otherwise
Notes:    Parametric coords range from 0 to 1 for points located on VaVb and VaVc axis
**************************************************************************************/
int Parametric
(
   CPOINT  *A,   /* Point in plane                    */
   CPOINT  *Va,  /* Plane origin point                */
   VECTOR  *u,   /* Unit vector defining s axis       */
   VECTOR  *v,   /* Unit vector defining t axis       */
   double  *s,   /* Output parametric point on s axis */
   double  *t    /* Output parametric point on t axis */
)
{
   VECTOR w;           /* Vector from origin (Va) to intersection point */
   double Den, uv;     /* Denominator and dot(u,v)                      */

   w.x =  A->x-Va->x;  w.y =  A->y-Va->y;  w.z =  A->z-Va->z;

/* New assumes s and t are orthogonal and A is in plane */
   *s = DOT(w,*u);
   *t = DOT(w,*v);
   return 0;
}


/**************************************************************************************
Function: ParametricB
Purpose:  Determine the parametric coords of a point in plane
Returns:  Zero if coords are valid, non-zero otherwise
Notes:    Parametric coords range from 0 to 1 for points located on VaVb and VaVc axis
**************************************************************************************/
int ParametricB
(
   CPOINT  *Po,  /* Point in plane                    */
   CPOINT  *Va,  /* Origin point                      */
   CPOINT  *Vb,  /* Point defining s axis             */
   CPOINT  *Vc,  /* Point defining t axis             */
   double  *s,   /* Output parametric point on s axis */
   double  *t    /* Output parametric point on t axis */
)
{
   VECTOR w,u,v;       /* Vector from origin (Va) to intersection point */
   double Den, uv;     /* Denominator and dot(u,v)                      */

   w.x = Po->x-Va->x; w.y = Po->y-Va->y;  w.z = Po->z-Va->z;
   u.x = Vb->x-Va->x; u.y = Vb->y-Va->y;  u.z = Vb->z-Va->z;
   v.x = Vc->x-Va->x; v.y = Vc->y-Va->y;  v.z = Vc->z-Va->z;

/*   *s = DOT(u,w)/pow(MAG(u),2); Assumes orthogonal vectors
     *t = DOT(v,w)/pow(MAG(v),2); */

/* Compute denominator */
   uv = DOT(u,v); Den = uv*uv - DOT(u,u)*DOT(v,v);
   if( Den == 0.0 ) return 1;

/* Parametric coords of point in plane */
   *s = (uv * DOT(w,v) - DOT(v,v)*DOT(w,u))/Den;
   *t = (uv * DOT(w,u) - DOT(u,u)*DOT(w,v))/Den;


   return 0;
}


/**************************************************************************************
Function: PlanePlane
Purpose:  Determine intersection line of two planes
Returns:  Non-zero if planes are parallel
Notes:    The returned points are seperated by a distance of 2E6 and located
          either side of the point on the line of intersection nearest to Pb
**************************************************************************************/
int PlanePlane
(
   VECTOR  *na,  /* Normal to plane a                 */
   VECTOR  *nb,  /* Normal to plane b                 */
   CPOINT  *Pa,  /* point in plane a                  */
   CPOINT  *Pb,  /* point in plane b                  */
   CPOINT  *A,   /* Return point on intersection line */
   CPOINT  *B    /* Return point on intersection line */
)
{
   VECTOR n;
   double m, da, db, f;

   n.x = CROSSX(*na,*nb); n.y = CROSSY(*na,*nb); n.z = CROSSZ(*na,*nb);
   if( fabs(m = MAG(n)) < 1E-10 ) return 1; /* Parallel */
   n.x /= m; n.y /= m; n.z /= m;

   da = -DOT(*na,*Pa);
   db = -DOT(*nb,*Pb);

   if( fabs(n.x)>=fabs(n.y)&&fabs(n.x)>=fabs(n.z) )
   {
      if( (m = na->z*nb->y-na->y*nb->z) == 0 ) return 1;
      A->x = 0; A->y = (nb->z*da-na->z*db)/m; A->z = (na->y*db-nb->y*da)/m;
   }
   else if( fabs(n.y)>=fabs(n.x)&&fabs(n.y)>=fabs(n.z) )
   {
      if( (m = na->z*nb->x-na->x*nb->z) == 0 ) return 1;
      A->x = (nb->z*da-na->z*db)/m; A->y = 0; A->z = (na->x*db-nb->x*da)/m;
   }
   else
   {
      if( (m = na->x*nb->y-nb->x*na->y) == 0 ) return 1;
      A->x = (na->y*db-nb->y*da)/m; A->y = (nb->x*da-na->x*db)/m;  A->z = 0;
   }
   f = DOT(*Pb,n)-DOT(*A,n);
   A->x = A->x + n.x*(f-1E6);
   A->y = A->y + n.y*(f-1E6);
   A->z = A->z + n.z*(f-1E6);
   B->x = A->x + n.x*2E6;
   B->y = A->y + n.y*2E6;
   B->z = A->z + n.z*2E6;

   return OK;
}


/**************************************************************************************
Function: Inverse3x3
Purpose:  Compute inverse of a 3x3 matrix
Returns:  Non-zero on error
Notes:
  | a11 a12 a13 |-1             |   a33a22-a32a23  -(a33a12-a32a13)   a23a12-a22a13  |
  | a21 a22 a23 |    =  1/DET * | -(a33a21-a31a23)   a33a11-a31a13  -(a23a11-a21a13) |
  | a31 a32 a33 |               |   a32a21-a31a22  -(a32a11-a31a12)   a22a11-a21a12  |
**************************************************************************************/
int Inverse3x3( double (*In)[3], double (*Out)[3] )
{
   double D = Determinant3x3( In[0], In[1], In[2] );
   if( fabs(D) < SMALL ) return 1;

   Out[0][0] =  (In[2][2]*In[1][1]-In[1][2]*In[2][1])/D;
   Out[1][0] = -(In[2][2]*In[1][0]-In[1][2]*In[2][0])/D;
   Out[2][0] =  (In[2][1]*In[1][0]-In[1][1]*In[2][0])/D;

   Out[0][1] = -(In[2][2]*In[0][1]-In[0][2]*In[2][1])/D;
   Out[1][1] =  (In[2][2]*In[0][0]-In[0][2]*In[2][0])/D;
   Out[2][1] = -(In[2][1]*In[0][0]-In[0][1]*In[2][0])/D;

   Out[0][2] =  (In[1][2]*In[0][1]-In[0][2]*In[1][1])/D;
   Out[1][2] = -(In[1][2]*In[0][0]-In[0][2]*In[1][0])/D;
   Out[2][2] =  (In[1][1]*In[0][0]-In[0][1]*In[1][0])/D;

   return 0;
};


/**************************************************************************************
Function: Determinant3x3
Purpose:  Compute determinant of a 3x3 matrix
Returns:  Determinant
**************************************************************************************/
double Determinant3x3( double *A, double *B, double *C )
{
   return A[0] * (B[1]*C[2]-C[1]*B[2]) -
          B[0] * (A[1]*C[2]-A[2]*C[1]) +
          C[0] * (A[1]*B[2]-B[1]*A[2]);
}


