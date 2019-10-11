/**********************************************************************
Author:  Paul Spencer
File:    error.h
Notes:   Errors returned by all software modules.
***********************************************************************/
#ifndef ERROR_H
#define ERROR_H


typedef enum ERRORS
{
   OK  = 0,         /* No error                            */
   SPARSE_MALLOC,
   SPARSE_BOUND,
   ENTRY_RAY,
   PACKED_MALLOC,
   MEMORY,          /* Memory allocation (fatal)           */
   BOUND,           /* Bounds on variable exceeded         */
   DOVERFLOW,       /* Overflow of predefined buffer       */
   TIMEOUT,         /* Timeout in loop                     */
   TRACE_FAIL,      /* Failed to trace ray                 */
   PLANE_SECT,      /* Intersect Rx plane                  */


   END_ERR          /* Last error for string               */
} ERRORS;


#endif
