/*****************************************************************************
Author:     Paul Spencer
File:       matrix.c
Purpose:    Functions for manipulating packed and sparse two-dimensional
            matrices. The matrices are stored to enable easy integration
            with MatLAB via its MEX interface.
*****************************************************************************/

#include <string.h>

#include "matrices.h"
#include "types.h"
#include "error.h"

#define VERBOSE

#ifdef VERBOSE /* Include for debugging */
   extern int mexPrintf(const char *format, ...);
#endif


/*---------------------------------------------------------------------------
   Sparse matrix functions
-----------------------------------------------------------------------------*/


/*****************************************************************************
Function:   InitSparse
Purpose:    Initialise the sparse matrix structure
Returns:    Non-zero on error
Notes:      NB. This function must be called prior to calling AppendSparse()
            with NumElem set to the number of non-zero elements
*****************************************************************************/
int InitSparse
(
   SMATRIX *M,          /* Input sparse matrix                     */
   SINT     Rows,       /* Number of rows                          */
   SINT     Cols,       /* Number of cols                          */
   SINT     NumElem     /* Number of non-zero elements             */
)
{
   SINT i;

   M->MajorIx  = 0L;
   M->NumElem  = NumElem;
   M->NumCols  = Cols;
   M->NumRows  = Rows;

   if( Rows == 0 && Cols == 0 )
   {
      M->Major = NULL;
      M->Minor = NULL;
      M->Elem  = NULL;
      return SPARSE_MALLOC;
   }

   M->Major = (SINT   *)malloc( (M->NumCols+1)  * sizeof(SINT) );
   M->Minor = (SINT   *)malloc( (M->NumElem+1)  * sizeof(SINT) );
   M->Elem  = (double *)malloc( (M->NumElem+1)  * sizeof(double));

   if( M->Major == NULL || M->Minor == NULL || M->Elem == NULL )
      return SPARSE_MALLOC;

   for( i = 0; i < M->NumCols+1; i++ ) M->Major[i] = 0L;
   M->NumElem = 0;

   return OK;
}


/*****************************************************************************
Function:   LoadSparse
Purpose:    Load a sparse matrix structure from matlab pointers
Returns:    Non-zero on error
Notes:
*****************************************************************************/
void LoadSparse
(
   SMATRIX *M,          /* Input sparse matrix                     */
   SINT     Rows,       /* Number of rows                          */
   SINT     Cols,       /* Number of cols                          */
   SINT     NumElem,    /* Number of non-zero elements             */
   SINT    *pMajor,
   SINT    *pMinor,
   double  *pElem
)
{
   M->MajorIx  = 0;
   M->NumElem  = NumElem;
   M->NumCols  = Cols;
   M->NumRows  = Rows;
   M->Major    = pMajor;
   M->Minor    = pMinor;
   M->Elem     = pElem;
}


/*****************************************************************************
Function:   FreeSparse
Purpose:    Free allocated memory for a sparse matrix
Notes:
*****************************************************************************/
void FreeSparse
(
   SMATRIX *M           /* Input sparse matrix                     */
)
{
   if( M->Major != NULL ) free( M->Major );
   if( M->Minor != NULL ) free( M->Minor );
   if( M->Elem  != NULL ) free( M->Elem  );
   M->NumRows = 0; M->Major = NULL;
   M->NumCols = 0; M->Minor = NULL;
   M->NumElem = 0; M->Elem  = NULL;
}


/*****************************************************************************
Function:   AppendSparse
Purpose:    Sequentially append data to a sparse matrix, the major index will
            be incremented after a call to this function
Returns:    Non-zero on error
Notes:      This function assumes that all data are appended sequentially, even
            all zero's (using Num = 0).
*****************************************************************************/
int AppendSparse
(
   SMATRIX *M,          /* Input sparse matrix to be updated      */
   double  *Data,       /* Input data to be appended              */
   SINT    *Minor,      /* Input minor indices of Data elements   */
   SINT     Num         /* Input number of terms to be appended   */
)
{
   M->MajorIx += 1;
   if( M->MajorIx > M->NumCols || Num > M->NumRows ) return SPARSE_BOUND;
   M->Major[M->MajorIx] = M->NumElem + Num;

   if( Num == 0 ) return( OK );
/*
   if( (M->Elem = (double *)realloc(M->Elem,(sizeof(double)*(M->NumElem+Num)))) == NULL )
      return SPARSE_MALLOC;
   if( (M->Minor = (SINT *)realloc(M->Minor,(sizeof(SINT)*(M->NumElem+Num)))) == NULL )
      return SPARSE_MALLOC;
*/
   memcpy( &M->Elem[M->NumElem],  Data, sizeof(double)*Num );
   memcpy( &M->Minor[M->NumElem], Minor, sizeof(SINT)*Num );

   M->NumElem += Num;

   return( OK );
}


/*****************************************************************************
Function:   MultiplySparse
Purpose:    Multiply two sparse matrices given sparsity of result
Returns:    Non-zero on error
Notes:      R must be allocated and release externally and is the same as S
            in terms of space requirements
*****************************************************************************/
int MultiplySparse
(
   SMATRIX *S,          /* Matrix defining sparsity structure     */
   SMATRIX *A,          /* Matrix A transposed                    */
   SMATRIX *B,          /* Matrix B                               */
   SMATRIX *R           /* Result, R = (A.' * B ).*(S~=0)         */
)
{
   int Err, j;
   double *pData;       /* Temp storage for column */
   SINT   *pMinor;      /* Temp storage for column */

   if( (pData  = (double *)malloc( S->NumRows * sizeof(double))) == NULL )
      return SPARSE_MALLOC;
   if( (pMinor = (SINT   *)malloc( S->NumRows * sizeof(SINT  ))) == NULL )
      return SPARSE_MALLOC;

   for( j = 0; j < S->NumCols; j++ )      /* Loop through cols in S                */
   {
      int e = S->Major[j];                /* First element in S for col j          */
      int k = 0;                          /* Index of elements appended            */
      memset(pData,0,S->NumRows*sizeof(double));

      while( e < S->Major[j+1] )          /* Loop till last element in S for col j */
      {
         int i  = S->Minor[e];            /* Row in S                              */
         int eB = B->Major[j];            /* First element in B for col j          */
         int eA = A->Major[i];            /* First element in A for col j          */

         while( eB < B->Major[j+1] && eA < A->Major[i+1] )
         {
            while( eA < A->Major[i+1] && A->Minor[eA] < B->Minor[eB] ) eA++;

            if( eA != A->Major[i+1] && A->Minor[eA] == B->Minor[eB] )
               pData[k] += A->Elem[eA]*B->Elem[eB];
            eB++;
         }

         pMinor[k++] = i; e++;
      }
      if( (Err = AppendSparse( R, pData, pMinor, k )) != OK ) return Err;
   }

   return( OK );
}



/*---------------------------------------------------------------------------
   Packed matrix functions
-----------------------------------------------------------------------------*/


/*****************************************************************************
Function:   InitPacked
Purpose:    Initialise a packed matrix
Returns:    Non-zero on error
Notes:      The elements of the matrix may be optionally provided in Elem,
            otherwise they are set to zero
*****************************************************************************/
int InitPacked
(
   PMATRIX *M,          /* Input packed matrix                           */
   SINT     Rows,       /* Input number of rows                          */
   SINT     Cols,       /* Input number of cols                          */
   double  *Elem        /* Input pointer to column major data or NULL    */
)
{
   M->Rows = Rows;
   M->Cols = Cols;

   if( Rows == 0 || Cols == 0 )
   {
      M->Elem = NULL;
      return( OK );
   }

   if( (M->Elem = (double *)calloc( (Rows*Cols+1),sizeof(double) )) == NULL )
      return PACKED_MALLOC;

   if( Elem != NULL ) memcpy( M->Elem, Elem, Rows*Cols*sizeof(double) );
   else               memset( M->Elem, 0, (Rows*Cols+1)*sizeof(double ) );

   return( OK );
}

/*****************************************************************************
Function:   LoadPacked
Purpose:    Load a packed matrix structure from matlab pointers
Returns:    Non-zero on error
Notes:
*****************************************************************************/
void LoadPacked
(
   PMATRIX *M,          /* Input packed matrix                     */
   SINT     Rows,       /* Number of rows                          */
   SINT     Cols,       /* Number of cols                          */
   double  *pElem       /* Pointer to data                         */
)
{
   M->Cols  = Cols;
   M->Rows  = Rows;
   M->Elem  = pElem;
}

/*****************************************************************************
Function:   ReallocPacked
Purpose:    Reallocate a packed matrix
Returns:    Non-zero on error
Notes:      The size of the matrix can be made smaller or larger
*****************************************************************************/
int ReallocPacked
(
   PMATRIX *A,          /* Input packed matrix returned with B appended */
   SINT Cols            /* New number of columns for the matrix         */
)
{
   if( Cols == 0 || A->Rows == 0 ) {A->Cols = 0; return OK;}
   A->Elem = (double *)realloc( A->Elem, A->Rows*(Cols)*sizeof(double) );

   if( A->Elem == NULL ) return PACKED_MALLOC;
   A->Cols = Cols;
   return OK;
}

/*****************************************************************************
Function:   CopyRows
Purpose:    Copy a number of rows from one matrix to another
*****************************************************************************/
void CopyRows
(
   PMATRIX  A,       /* Input matrix                     */
   PMATRIX *B,       /* Matrix to change                 */
   SINT     RowA,    /* Start row in Ma to copy from     */
   SINT     RowB,    /* Start row in Mb to copy to       */
   SINT     NumRows  /* Number of rows to copy           */
)
{
   SINT  i;

   for( i = 0; i < A.Cols; i++ )
   {
      memcpy( &Elem(*B,RowB,i), &Elem(A,RowA,i), sizeof(double)*NumRows );
   }
}


/*****************************************************************************
Function:   FreePacked
Purpose:    Free a packed matrix
*****************************************************************************/
void FreePacked
(
   PMATRIX *M           /* Input packed matrix                     */
)
{
   if( M->Elem != NULL ) free( M->Elem );
   M->Elem = NULL;
}


