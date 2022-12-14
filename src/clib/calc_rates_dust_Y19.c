#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "grackle_macros.h"
#include "grackle_types.h"
#include "grackle_chemistry_data.h"
#include "phys_constants.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#define tiny 1.0e-20
#define huge 1.0e+20
#define tevk 1.1605e+4

extern int grackle_verbose;

int calc_rates_dust_Y19(int iSN, chemistry_data *my_chemistry, chemistry_data_storage *my_rates)
{

  int NTd, Nmom;
  int iTd, imom, itab0, itab;

  my_chemistry->SN0_XC [iSN] =   2.50000e-01;
  my_chemistry->SN0_XO [iSN] =   2.93867e-01;
  my_chemistry->SN0_XMg[iSN] =   6.00000e-02;
  my_chemistry->SN0_XAl[iSN] =   2.85361e-03;
  my_chemistry->SN0_XSi[iSN] =   7.00000e-02;
  my_chemistry->SN0_XS [iSN] =   1.58191e-02;
  my_chemistry->SN0_XFe[iSN] =   6.64078e-02;

  my_chemistry->SN0_fC [iSN] =   0.00000e+00;
  my_chemistry->SN0_fO [iSN] =   1.73867e-01;
  my_chemistry->SN0_fMg[iSN] =   0.00000e+00;
  my_chemistry->SN0_fAl[iSN] =   2.85361e-03;
  my_chemistry->SN0_fSi[iSN] =   0.00000e+00;
  my_chemistry->SN0_fS [iSN] =   1.58191e-02;
  my_chemistry->SN0_fFe[iSN] =   6.64078e-02;

  my_chemistry->SN0_fMgSiO3  [iSN] =   2.50000e-01;
  my_chemistry->SN0_fAC      [iSN] =   2.50000e-01;

  itab0 = 3 * iSN;
  my_chemistry->SN0_r0MgSiO3  [itab0 + 0] =   1.00000e-05;
  my_chemistry->SN0_r0AC      [itab0 + 0] =   1.00000e-05;

  my_chemistry->SN0_r0MgSiO3  [itab0 + 1] =   1.00000e-10;
  my_chemistry->SN0_r0AC      [itab0 + 1] =   1.00000e-10;

  my_chemistry->SN0_r0MgSiO3  [itab0 + 2] =   1.00000e-15;
  my_chemistry->SN0_r0AC      [itab0 + 2] =   1.00000e-15;

  NTd =            35;
 Nmom =             4;

  double Y19_kpMgSiO3[] = 
  {  2.19890e-02,   2.19890e-07,   2.19890e-12,   2.19890e-17,
     3.90612e-02,   3.90612e-07,   3.90612e-12,   3.90612e-17,
     6.05539e-02,   6.05539e-07,   6.05539e-12,   6.05539e-17,
     8.76116e-02,   8.76116e-07,   8.76116e-12,   8.76116e-17,
     1.43288e-01,   1.43288e-06,   1.43288e-11,   1.43288e-16,
     2.19266e-01,   2.19266e-06,   2.19266e-11,   2.19266e-16,
     3.36256e-01,   3.36256e-06,   3.36256e-11,   3.36256e-16,
     5.14336e-01,   5.14336e-06,   5.14336e-11,   5.14336e-16,
     7.97216e-01,   7.97216e-06,   7.97216e-11,   7.97216e-16,
     1.25414e+00,   1.25414e-05,   1.25414e-10,   1.25414e-15,
     2.03450e+00,   2.03450e-05,   2.03450e-10,   2.03450e-15,
     3.34654e+00,   3.34654e-05,   3.34654e-10,   3.34654e-15,
     5.45913e+00,   5.45913e-05,   5.45913e-10,   5.45913e-15,
     8.82166e+00,   8.82166e-05,   8.82166e-10,   8.82166e-15,
     1.41836e+01,   1.41836e-04,   1.41836e-09,   1.41836e-14,
     2.28449e+01,   2.28449e-04,   2.28449e-09,   2.28449e-14,
     3.71258e+01,   3.71258e-04,   3.71258e-09,   3.71258e-14,
     6.14485e+01,   6.14485e-04,   6.14485e-09,   6.14485e-14,
     1.03898e+02,   1.03898e-03,   1.03898e-08,   1.03898e-13,
     1.75627e+02,   1.75627e-03,   1.75627e-08,   1.75627e-13,
     2.82290e+02,   2.82290e-03,   2.82290e-08,   2.82290e-13,
     4.14908e+02,   4.14908e-03,   4.14908e-08,   4.14908e-13,
     5.60606e+02,   5.60606e-03,   5.60606e-08,   5.60606e-13,
     7.12020e+02,   7.12020e-03,   7.12020e-08,   7.12020e-13,
     8.42130e+02,   8.42130e-03,   8.42130e-08,   8.42130e-13,
     8.96812e+02,   8.96812e-03,   8.96812e-08,   8.96812e-13,
     8.41845e+02,   8.41845e-03,   8.41845e-08,   8.41845e-13,
     6.97883e+02,   6.97883e-03,   6.97883e-08,   6.97883e-13,
     5.19082e+02,   5.19082e-03,   5.19082e-08,   5.19082e-13,
     3.53464e+02,   3.53464e-03,   3.53464e-08,   3.53464e-13,
     2.24610e+02,   2.24610e-03,   2.24610e-08,   2.24610e-13,
     1.35389e+02,   1.35389e-03,   1.35389e-08,   1.35389e-13,
     7.84898e+01,   7.84898e-04,   7.84898e-09,   7.84898e-14,
     4.43113e+01,   4.43113e-04,   4.43113e-09,   4.43113e-14,
     2.49396e+01,   2.49396e-04,   2.49396e-09,   2.49396e-14  };

  double Y19_kpAC[] = 
  {  6.76020e-02,   6.76020e-07,   6.76020e-12,   6.76020e-17,
     1.20181e-01,   1.20181e-06,   1.20181e-11,   1.20181e-16,
     1.86375e-01,   1.86375e-06,   1.86375e-11,   1.86375e-16,
     2.69708e-01,   2.69708e-06,   2.69708e-11,   2.69708e-16,
     4.44368e-01,   4.44368e-06,   4.44368e-11,   4.44368e-16,
     6.87406e-01,   6.87406e-06,   6.87406e-11,   6.87406e-16,
     1.07797e+00,   1.07797e-05,   1.07797e-10,   1.07797e-15,
     1.71241e+00,   1.71241e-05,   1.71241e-10,   1.71241e-15,
     2.74163e+00,   2.74163e-05,   2.74163e-10,   2.74163e-15,
     4.35812e+00,   4.35812e-05,   4.35812e-10,   4.35812e-15,
     6.98720e+00,   6.98720e-05,   6.98720e-10,   6.98720e-15,
     1.13206e+01,   1.13206e-04,   1.13206e-09,   1.13206e-14,
     1.85159e+01,   1.85159e-04,   1.85159e-09,   1.85159e-14,
     3.09414e+01,   3.09414e-04,   3.09414e-09,   3.09414e-14,
     5.34575e+01,   5.34575e-04,   5.34575e-09,   5.34575e-14,
     9.60912e+01,   9.60912e-04,   9.60912e-09,   9.60912e-14,
     1.76000e+02,   1.76000e-03,   1.76000e-08,   1.76000e-13,
     3.10598e+02,   3.10598e-03,   3.10598e-08,   3.10598e-13,
     4.95502e+02,   4.95502e-03,   4.95502e-08,   4.95502e-13,
     6.87389e+02,   6.87389e-03,   6.87389e-08,   6.87389e-13,
     8.21760e+02,   8.21760e-03,   8.21760e-08,   8.21760e-13,
     8.55209e+02,   8.55209e-03,   8.55209e-08,   8.55209e-13,
     7.90315e+02,   7.90315e-03,   7.90315e-08,   7.90315e-13,
     6.63814e+02,   6.63814e-03,   6.63814e-08,   6.63814e-13,
     5.19410e+02,   5.19410e-03,   5.19410e-08,   5.19410e-13,
     3.88956e+02,   3.88956e-03,   3.88956e-08,   3.88956e-13,
     2.88141e+02,   2.88141e-03,   2.88141e-08,   2.88141e-13,
     2.20698e+02,   2.20698e-03,   2.20698e-08,   2.20698e-13,
     1.84716e+02,   1.84716e-03,   1.84716e-08,   1.84716e-13,
     1.78316e+02,   1.78316e-03,   1.78316e-08,   1.78316e-13,
     2.05010e+02,   2.05010e-03,   2.05010e-08,   2.05010e-13,
     2.82760e+02,   2.82760e-03,   2.82760e-08,   2.82760e-13,
     4.70437e+02,   4.70437e-03,   4.70437e-08,   4.70437e-13,
     9.49808e+02,   9.49808e-03,   9.49808e-08,   9.49808e-13,
     2.18634e+03,   2.18634e-02,   2.18634e-07,   2.18634e-12  };


  itab0 = Nmom * NTd * iSN;
  itab  = 0;
  for(imom = 0; imom < Nmom; imom++) {
    for(iTd = 0; iTd < NTd; iTd++) {
      my_rates->SN0_kpMgSiO3  [itab0] = Y19_kpMgSiO3  [itab];
      my_rates->SN0_kpAC      [itab0] = Y19_kpAC      [itab];
      itab0++;
      itab ++;
    }
  }

  return SUCCESS;
}
