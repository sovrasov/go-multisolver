#ifndef MAP_HPP
#define MAP_HPP


enum MapType {
  Simple = 1, Linear = 2, Noninjective = 3
};

void mapd(double x, int m, double* y, int n, int key = 1);    /* map x to y         */
void invmad(int, double *, int, int *, double *, int, int);  /* map y to x         */
void xyd(double *xx, int m, double* y, int n);        /* get preimage       */

#endif
