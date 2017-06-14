#include<bits/stdc++.h>
#include"simplex.h"

#define EPS 1e-9 /* 10^-9 */
#define correct_double(var) \
  (fabs(var) < EPS ? 0. : var) /* Avoid printing -0.000 */
#define FOR(var, a, b) \
  for (var = (a); var <= (b); ++var) /* Avoid silly mistakes ;) */

double A[MAX][MAX];

int main(void) {
  int n, m, i, j, ans;
  double z[MAX], b[MAX], z0, x[MAX];
  scanf("%d %d", &m, &n);
  FOR(j, 1, n) scanf("%lf", &z[j]);
  FOR(i, 1, m) {
    FOR(j, 1, n) scanf("%lf", &A[i][j]);
    scanf("%lf", &b[i]);
  }
  ans = simplex(m, n, z, A, b, &z0, x);
  switch (ans) {
  case NFEA: printf("infeasible\n"); break;
  case UNBD: printf("unbounded\n"); break;
  case FEA:
    printf("%.5lf\n", correct_double(z0));
    FOR(j, 1, n)
      printf("%s%.8lf", j > 1 ? " " : "", correct_double(x[j]));
    printf("\n");
  }
  return 0;
}


