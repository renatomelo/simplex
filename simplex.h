#ifndef _SIMPLEX_
#define _SIMPLEX_

#define MAX 1123 /* ~10^3      */
#define FEA   0  /* FEAsible   */
#define NFEA -1  /* iNFEAsible */
#define UNBD  1  /* UNBounDed  */

int simplex(int m, int n, double z[], double A[][MAX], double b[],
            double *z0, double x[]);

#endif
