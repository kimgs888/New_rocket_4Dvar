/***************************************************************************
Author:  Paul Spencer
File:    matrices.h
Notes:   Definitions for sparse and packed two-dimensional matrices.
         Matrices are stored so as to enable easy integration with matlab.
****************************************************************************/
#ifndef MATRICES_H
#define MATRICES_H

#include <math.h>
#include <stdlib.h>

#include "types.h"
#include "error.h"

/* Type for element indexes */
#define SINT int

/* Macro to return the element of a packed matrix, PMatrix         */
/* Usage:  x = Elem(PMatrix,row,col)                               */
#define Elem(M,i,j) ((M).Elem[(j)*(M).Rows+(i)])

/* Macro to return TRUE if packed or sparse matrix is empty        */
#define IsEmpty(M) ((M).Elem == NULL||(M).Rows==0||(M).Cols==0 )

/* Macro to return TRUE if sparse matrix is fully built            */
#define SpDone(M) ((M).NumCols == (M).MajorIx)

typedef struct PMATRIX  /* Packed matrix definition                */
{
   SINT    Rows;        /* Number of rows in the matrix            */
   SINT    Cols;        /* Number of columns in the matrix         */
   double *Elem;        /* Vector of column major elements         */
} PMATRIX;

typedef struct CELL     /* Cell array of basis maps                */
{
   SINT    Size;        /* Number of cells                         */
   SINT    Rows;        /* Number of expanded rows                 */
   SINT    Cols;        /* Number of expanded cols                 */
   PMATRIX M[4];        /* Cell matrices                           */
} CELL;

typedef struct SMATRIX  /* Sparse matrix definition                */
{
   SINT    NumCols;     /* Number of columns                       */
   SINT    NumRows;     /* Number of rows                          */
   SINT    NumElem;     /* Number of non-zero elements             */
   double *Elem;        /* Vector of non-zero elements             */
   SINT   *Major;       /* End indices for major axis              */
   SINT   *Minor;       /* Minor indices for all elements          */
   SINT    MajorIx;     /* Current index of major for append       */
} SMATRIX;


int InitSparse
(
   SMATRIX *M,          /* Input sparse matrix                     */
   SINT     Rows,       /* Number of rows                          */
   SINT     Cols,       /* Number of cols                          */
   SINT     NumElem     /* Number of non-zero elements             */
);

void LoadSparse
(
   SMATRIX *M,          /* Input sparse matrix                     */
   SINT     Rows,       /* Number of rows                          */
   SINT     Cols,       /* Number of cols                          */
   SINT     NumElem,    /* Number of non-zero elements             */
   SINT    *pMajor,
   SINT    *pMinor,
   double  *pElem
);


int MultiplySparse
(
   SMATRIX *S,          /* Matrix defining sparsity structure     */
   SMATRIX *A,          /* Matrix A                               */
   SMATRIX *B,          /* Matrix B                               */
   SMATRIX *R           /* Result, R = (A * B ).*(S~=0)           */
);


void FreeSparse
(
   SMATRIX *M           /* Input sparse matrix                     */
);

int AppendSparse
(
   SMATRIX *M,          /* Sparse matrix to be updated             */
   double  *Data,       /* Data to be appended                     */
   SINT    *Minor,      /* Minor indices of Data elements          */
   SINT     Num         /* Number of terms to be appended          */
);

int InitPacked
(
   PMATRIX *M,          /* Input packed matrix                     */
   SINT     Rows,       /* Number of rows                          */
   SINT     Cols,       /* Number of cols                          */
   double  *Elem        /* Pointer to column major data or NULL    */
);

void LoadPacked
(
   PMATRIX *M,          /* Input packed matrix                     */
   SINT     Rows,       /* Number of rows                          */
   SINT     Cols,       /* Number of cols                          */
   double  *pElem       /* Pointer to data                         */
);

int ReallocPacked
(
   PMATRIX *A,          /* Input packed matrix returned with B appended */
   SINT Cols             /* Extra columns to append to the end           */
);

void CopyRows
(
   PMATRIX  A,          /* Input matrix                            */
   PMATRIX *B,          /* Matrix to change                        */
   SINT     RowA,       /* Start row in Ma to copy from            */
   SINT     RowB,       /* Start row in Mb to copy to              */
   SINT     NumRows     /* Number of rows to copy                  */
);

void FreePacked
(
   PMATRIX *M           /* Input packed matrix                     */
);


int CellMat
(
   CELL    *A,          /* Input cell array of matrices            */
   PMATRIX *B,          /* Input packed matrix                     */
   PMATRIX *C           /* Input packed matrix                     */
);

void CellCol
(
   CELL   *A,            /* Input cell array of matrices            */
   int     Col,          /* Input required column                   */
   double *pCol          /* Output column vector                    */
);


#endif
