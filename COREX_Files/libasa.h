/*----------------------------------------*
 * This is the asa_lib header file
 *----------------------------------------*/

#define ASA_VERSION "27-Mar-96"

#ifdef M_PI
  #define PI	         M_PI
#else
  #define PI             3.1415926535897932
#endif
#define TWOPI            (2*PI)
#define EX_SOFTWARE      1
#define ABS(x)           ((x) < 0.0 ? -(x) : (x))
#define SQ(x)            ((x) * (x))

int surfinit();
double surfarea();
void surfcleanup();

#include "asa_types.h"
