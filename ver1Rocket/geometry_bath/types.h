/********************************************************************
Author:  Paul Spencer
File:    types.h
Notes:   General type definitions
*********************************************************************/
#ifndef TYPES_H
#define TYPES_H

/* #ifndef _WINDEF_ */
typedef signed short SWORD;
typedef unsigned short UWORD;
typedef unsigned long ULONG;
typedef signed long SLONG;
typedef unsigned int UINT;
/* typedef unsigned int BOOL; */
/*#endif */
typedef double CALC;

#define TRUE  1
#define FALSE 0
#define ZERO  1E-12

#include <math.h>
#ifndef M_PI
   #define M_PI 3.14159265358979
#endif


typedef struct RANGE      /* Range definition                          */
{
   double Min;            /* Minimum                                   */
   double Max;            /* Maximum                                   */
} RANGE;


#endif
