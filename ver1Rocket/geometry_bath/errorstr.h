/**********************************************************************
Author:  Paul Spencer
File:    errorstr.h
Notes:   Error strings returned by each function.  Include this file
         in the entry point function.
***********************************************************************/
#ifndef ERRORSTR_H
#define ERRORSTR_H

#include "error.h"




char ErrString[END_ERR][40] =
{
   "                                      |",
   "Insufficient memory for sparse matrix"
   "Insufficient memory for packed matrix",
   "Ray dosn't intersect grid"

};


#endif
