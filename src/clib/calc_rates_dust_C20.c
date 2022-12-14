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

int calc_rates_dust_C20(int iSN, chemistry_data *my_chemistry, chemistry_data_storage *my_rates)
{

  int NTd, Nmom;
  int iTd, imom, itab0, itab;

  my_chemistry->SN0_XC [iSN] =   1.00183e-01;
  my_chemistry->SN0_XO [iSN] =   6.06515e-01;
  my_chemistry->SN0_XMg[iSN] =   2.75968e-02;
  my_chemistry->SN0_XAl[iSN] =   1.87118e-04;
  my_chemistry->SN0_XSi[iSN] =   1.00051e-01;
  my_chemistry->SN0_XS [iSN] =   6.02208e-02;
  my_chemistry->SN0_XFe[iSN] =   3.07560e-02;

  my_chemistry->SN0_fC [iSN] =   8.74563e-02;
  my_chemistry->SN0_fO [iSN] =   6.04383e-01;
  my_chemistry->SN0_fMg[iSN] =   2.63753e-02;
  my_chemistry->SN0_fAl[iSN] =   1.87118e-04;
  my_chemistry->SN0_fSi[iSN] =   6.44592e-02;
  my_chemistry->SN0_fS [iSN] =   6.02018e-02;
  my_chemistry->SN0_fFe[iSN] =   2.69505e-02;

  my_chemistry->SN0_fSiM     [iSN] =   3.44388e-02;
  my_chemistry->SN0_fFeM     [iSN] =   3.77223e-03;
  my_chemistry->SN0_fMg2SiO4 [iSN] =   1.90086e-03;
  my_chemistry->SN0_fMgSiO3  [iSN] =   2.57266e-06;
  my_chemistry->SN0_fAC      [iSN] =   1.27270e-02;
  my_chemistry->SN0_fSiO2D   [iSN] =   1.65484e-03;
  my_chemistry->SN0_fMgO     [iSN] =   9.48713e-04;
  my_chemistry->SN0_fFeS     [iSN] =   5.23050e-05;
  my_chemistry->SN0_fAl2O3   [iSN] =   1.31693e-29;

  itab0 = 3 * iSN;
  my_chemistry->SN0_r0SiM     [itab0 + 0] =   1.24861e-05;
  my_chemistry->SN0_r0FeM     [itab0 + 0] =   6.67024e-06;
  my_chemistry->SN0_r0Mg2SiO4 [itab0 + 0] =   1.41253e-06;
  my_chemistry->SN0_r0MgSiO3  [itab0 + 0] =   1.01138e-06;
  my_chemistry->SN0_r0AC      [itab0 + 0] =   7.95099e-07;
  my_chemistry->SN0_r0SiO2D   [itab0 + 0] =   1.40285e-06;
  my_chemistry->SN0_r0MgO     [itab0 + 0] =   1.29303e-06;
  my_chemistry->SN0_r0FeS     [itab0 + 0] =   1.68897e-06;
  my_chemistry->SN0_r0Al2O3   [itab0 + 0] =   9.21063e-08;

  my_chemistry->SN0_r0SiM     [itab0 + 1] =   2.86508e-10;
  my_chemistry->SN0_r0FeM     [itab0 + 1] =   7.50596e-11;
  my_chemistry->SN0_r0Mg2SiO4 [itab0 + 1] =   4.77566e-12;
  my_chemistry->SN0_r0MgSiO3  [itab0 + 1] =   1.31688e-12;
  my_chemistry->SN0_r0AC      [itab0 + 1] =   2.51133e-12;
  my_chemistry->SN0_r0SiO2D   [itab0 + 1] =   3.98828e-12;
  my_chemistry->SN0_r0MgO     [itab0 + 1] =   1.06240e-11;
  my_chemistry->SN0_r0FeS     [itab0 + 1] =   3.16618e-12;
  my_chemistry->SN0_r0Al2O3   [itab0 + 1] =   9.03508e-15;

  my_chemistry->SN0_r0SiM     [itab0 + 2] =   1.01028e-14;
  my_chemistry->SN0_r0FeM     [itab0 + 2] =   1.22752e-15;
  my_chemistry->SN0_r0Mg2SiO4 [itab0 + 2] =   3.08016e-17;
  my_chemistry->SN0_r0MgSiO3  [itab0 + 2] =   2.89696e-18;
  my_chemistry->SN0_r0AC      [itab0 + 2] =   4.21640e-17;
  my_chemistry->SN0_r0SiO2D   [itab0 + 2] =   1.93974e-17;
  my_chemistry->SN0_r0MgO     [itab0 + 2] =   1.57687e-16;
  my_chemistry->SN0_r0FeS     [itab0 + 2] =   6.72598e-18;
  my_chemistry->SN0_r0Al2O3   [itab0 + 2] =   9.36936e-22;

  NTd =            35;
 Nmom =             4;

  double C20_kpSiM[] = 
  {  1.53894e-01,   1.90916e-06,   4.34900e-11,   1.52207e-15,
     1.93844e-01,   2.40648e-06,   5.48634e-11,   1.92178e-15,
     2.44138e-01,   3.03256e-06,   6.91797e-11,   2.42474e-15,
     3.07454e-01,   3.82073e-06,   8.72020e-11,   3.05783e-15,
     3.87243e-01,   4.81526e-06,   1.09978e-10,   3.85936e-15,
     4.87709e-01,   6.06778e-06,   1.38669e-10,   4.86916e-15,
     6.14251e-01,   7.64642e-06,   1.74856e-10,   6.14383e-15,
     7.73625e-01,   9.63590e-06,   2.20495e-10,   7.75268e-15,
     9.74036e-01,   1.21392e-05,   2.77959e-10,   9.78002e-15,
     1.22376e+00,   1.52600e-05,   3.49642e-10,   1.23105e-14,
     1.52029e+00,   1.89682e-05,   4.34881e-10,   1.53223e-14,
     1.83650e+00,   2.29254e-05,   5.25941e-10,   1.85447e-14,
     2.15714e+00,   2.69443e-05,   6.18616e-10,   2.18354e-14,
     2.55729e+00,   3.19712e-05,   7.34846e-10,   2.59753e-14,
     3.23398e+00,   4.04866e-05,   9.32032e-10,   3.30022e-14,
     4.34499e+00,   5.44911e-05,   1.25683e-09,   4.45846e-14,
     5.84259e+00,   7.34292e-05,   1.69748e-09,   6.03345e-14,
     7.53509e+00,   9.49825e-05,   2.20255e-09,   7.84868e-14,
     9.31292e+00,   1.17996e-04,   2.75076e-09,   9.84428e-14,
     1.14053e+01,   1.46050e-04,   3.44207e-09,   1.24265e-13,
     1.44699e+01,   1.88865e-04,   4.53802e-09,   1.66375e-13,
     1.88525e+01,   2.52409e-04,   6.22006e-09,   2.32633e-13,
     2.35897e+01,   3.27106e-04,   8.33967e-09,   3.20293e-13,
     2.71065e+01,   3.99935e-04,   1.08043e-08,   4.33797e-13,
     2.87408e+01,   4.69399e-04,   1.38602e-08,   5.94648e-13,
     2.88428e+01,   5.34587e-04,   1.74785e-08,   8.07129e-13,
     2.80514e+01,   5.86377e-04,   2.09655e-08,   1.03092e-12,
     2.67893e+01,   6.13530e-04,   2.33810e-08,   1.20137e-12,
     2.53851e+01,   6.14601e-04,   2.43011e-08,   1.28034e-12,
     2.53848e+01,   6.23633e-04,   2.48207e-08,   1.31287e-12,
     3.52560e+01,   8.08661e-04,   3.01974e-08,   1.52398e-12,
     9.64731e+01,   1.84326e-03,   5.73265e-08,   2.51143e-12,
     3.56496e+02,   5.70937e-03,   1.44867e-07,   5.33930e-12,
     1.16716e+03,   1.65654e-02,   3.60906e-07,   1.15456e-11,
     3.02792e+03,   3.93985e-02,   7.71226e-07,   2.22345e-11  };

  double C20_kpFeM[] = 
  {  1.10506e-02,   1.69983e-07,   3.50089e-12,   9.16404e-17,
     1.85837e-02,   2.76971e-07,   5.53036e-12,   1.40642e-16,
     2.79449e-02,   4.08949e-07,   8.01809e-12,   2.00458e-16,
     3.96605e-02,   5.73574e-07,   1.11125e-11,   2.74730e-16,
     6.19392e-02,   8.69904e-07,   1.63567e-11,   3.93140e-16,
     9.16817e-02,   1.25867e-06,   2.31128e-11,   5.43193e-16,
     1.36080e-01,   1.82251e-06,   3.26005e-11,   7.47226e-16,
     2.02056e-01,   2.63811e-06,   4.59129e-11,   1.02492e-15,
     2.99593e-01,   3.81390e-06,   6.45602e-11,   1.40281e-15,
     4.39758e-01,   5.46636e-06,   9.01024e-11,   1.90718e-15,
     6.41903e-01,   7.79581e-06,   1.25163e-10,   2.58082e-15,
     9.24985e-01,   1.09870e-05,   1.71974e-10,   3.45658e-15,
     1.30585e+00,   1.51941e-05,   2.32250e-10,   4.55695e-15,
     1.80164e+00,   2.05666e-05,   3.07564e-10,   5.90101e-15,
     2.42648e+00,   2.72192e-05,   3.99035e-10,   7.50126e-15,
     3.19207e+00,   3.52384e-05,   5.07449e-10,   9.36601e-15,
     4.10565e+00,   4.46618e-05,   6.32978e-10,   1.14947e-14,
     5.17172e+00,   5.54917e-05,   7.75385e-10,   1.38820e-14,
     6.40478e+00,   6.78244e-05,   9.35760e-10,   1.65481e-14,
     7.85614e+00,   8.21045e-05,   1.11980e-09,   1.95943e-14,
     9.64785e+00,   9.94326e-05,   1.34166e-09,   2.32660e-14,
     1.20159e+01,   1.21941e-04,   1.62875e-09,   2.80374e-14,
     1.53778e+01,   1.53371e-04,   2.02922e-09,   3.47445e-14,
     2.04423e+01,   1.99938e-04,   2.62238e-09,   4.47607e-14,
     2.83887e+01,   2.71546e-04,   3.53030e-09,   6.01372e-14,
     4.11773e+01,   3.83718e-04,   4.93303e-09,   8.37180e-14,
     6.21107e+01,   5.61284e-04,   7.10254e-09,   1.19520e-13,
     9.68132e+01,   8.45164e-04,   1.04719e-08,   1.73689e-13,
     1.54804e+02,   1.30324e-03,   1.57446e-08,   2.55980e-13,
     2.51780e+02,   2.04562e-03,   2.40381e-08,   3.81486e-13,
     4.12838e+02,   3.24595e-03,   3.70714e-08,   5.72650e-13,
     6.76592e+02,   5.16747e-03,   5.73792e-08,   8.61246e-13,
     1.09863e+03,   8.17937e-03,   8.83910e-08,   1.28821e-12,
     1.74646e+03,   1.26954e-02,   1.33666e-07,   1.89228e-12,
     2.67882e+03,   1.89754e-02,   1.94753e-07,   2.68226e-12  };

  double C20_kpMg2SiO4[] = 
  {  1.05240e-01,   1.48654e-07,   5.02591e-13,   3.24156e-18,
     1.32588e-01,   1.87284e-07,   6.33194e-13,   4.08391e-18,
     1.67016e-01,   2.35915e-07,   7.97614e-13,   5.14437e-18,
     2.10360e-01,   2.97139e-07,   1.00461e-12,   6.47941e-18,
     2.71887e-01,   3.84048e-07,   1.29844e-12,   8.37459e-18,
     3.55694e-01,   5.02428e-07,   1.69868e-12,   1.09560e-17,
     4.84932e-01,   6.84981e-07,   2.31589e-12,   1.49369e-17,
     6.99767e-01,   9.88442e-07,   3.34188e-12,   2.15544e-17,
     1.05860e+00,   1.49530e-06,   5.05559e-12,   3.26077e-17,
     1.62902e+00,   2.30104e-06,   7.77981e-12,   5.01791e-17,
     2.54260e+00,   3.59152e-06,   1.21430e-11,   7.83229e-17,
     3.96490e+00,   5.60062e-06,   1.89363e-11,   1.22144e-16,
     6.10634e+00,   8.62558e-06,   2.91646e-11,   1.88128e-16,
     9.28774e+00,   1.31196e-05,   4.43615e-11,   2.86175e-16,
     1.39267e+01,   1.96729e-05,   6.65240e-11,   4.29191e-16,
     2.05387e+01,   2.90138e-05,   9.81191e-11,   6.33137e-16,
     3.00660e+01,   4.24748e-05,   1.43665e-10,   9.27312e-16,
     4.55128e+01,   6.43028e-05,   2.17560e-10,   1.40505e-15,
     7.47912e+01,   1.05684e-04,   3.57721e-10,   2.31202e-15,
     1.29637e+02,   1.83206e-04,   6.20346e-10,   4.01209e-15,
     2.14853e+02,   3.03661e-04,   1.02847e-09,   6.65466e-15,
     3.20060e+02,   4.52395e-04,   1.53258e-09,   9.92043e-15,
     4.29862e+02,   6.07661e-04,   2.05913e-09,   1.33350e-14,
     5.30972e+02,   7.50675e-04,   2.54443e-09,   1.64851e-14,
     5.99878e+02,   8.48166e-04,   2.87540e-09,   1.86348e-14,
     6.06737e+02,   8.57905e-04,   2.90867e-09,   1.88528e-14,
     5.43451e+02,   7.68444e-04,   2.60547e-09,   1.68884e-14,
     4.33701e+02,   6.13270e-04,   2.07941e-09,   1.34791e-14,
     3.13436e+02,   4.43222e-04,   1.50292e-09,   9.74306e-15,
     2.09092e+02,   2.95687e-04,   1.00276e-09,   6.50192e-15,
     1.31208e+02,   1.85574e-04,   6.29536e-10,   4.08396e-15,
     7.91293e+01,   1.11988e-04,   3.80343e-10,   2.47142e-15,
     4.74361e+01,   6.72978e-05,   2.29477e-10,   1.49906e-15,
     2.98121e+01,   4.25812e-05,   1.46690e-10,   9.70392e-16,
     2.07674e+01,   3.01013e-05,   1.05718e-10,   7.13724e-16  };

  double C20_kpMgSiO3[] = 
  {  2.19890e-02,   2.22393e-08,   2.89570e-14,   6.37012e-20,
     3.90612e-02,   3.95059e-08,   5.14391e-14,   1.13159e-19,
     6.05539e-02,   6.12433e-08,   7.97425e-14,   1.75423e-19,
     8.76116e-02,   8.86090e-08,   1.15374e-13,   2.53808e-19,
     1.43288e-01,   1.44919e-07,   1.88693e-13,   4.15103e-19,
     2.19266e-01,   2.21762e-07,   2.88749e-13,   6.35214e-19,
     3.36256e-01,   3.40084e-07,   4.42810e-13,   9.74136e-19,
     5.14336e-01,   5.20191e-07,   6.77322e-13,   1.49005e-18,
     7.97217e-01,   8.06292e-07,   1.04985e-12,   2.30962e-18,
     1.25414e+00,   1.26842e-06,   1.65157e-12,   3.63350e-18,
     2.03450e+00,   2.05766e-06,   2.67923e-12,   5.89467e-18,
     3.34648e+00,   3.38458e-06,   4.40703e-12,   9.69677e-18,
     5.45894e+00,   5.52111e-06,   7.18906e-12,   1.58199e-17,
     8.82120e+00,   8.92167e-06,   1.16172e-11,   2.55686e-17,
     1.41826e+01,   1.43441e-05,   1.86785e-11,   4.11217e-17,
     2.28421e+01,   2.31025e-05,   3.00847e-11,   6.62628e-17,
     3.71183e+01,   3.75417e-05,   4.88915e-11,   1.07761e-16,
     6.14292e+01,   6.21307e-05,   8.09232e-11,   1.78540e-16,
     1.03850e+02,   1.05037e-04,   1.36825e-10,   3.02237e-16,
     1.75513e+02,   1.77524e-04,   2.31276e-10,   5.11420e-16,
     2.82073e+02,   2.85307e-04,   3.71726e-10,   8.22603e-16,
     4.14541e+02,   4.19298e-04,   5.46336e-10,   1.20964e-15,
     5.60007e+02,   5.66440e-04,   7.38095e-10,   1.63478e-15,
     7.11090e+02,   7.19267e-04,   9.37240e-10,   2.07555e-15,
     8.40892e+02,   8.50566e-04,   1.10826e-09,   2.45241e-15,
     8.95414e+02,   9.05715e-04,   1.18001e-09,   2.60842e-15,
     8.40504e+02,   8.50170e-04,   1.10753e-09,   2.44576e-15,
     6.96768e+02,   7.04778e-04,   9.18057e-10,   2.02574e-15,
     5.18260e+02,   5.24216e-04,   6.82816e-10,   1.50584e-15,
     3.52903e+02,   3.56958e-04,   4.64940e-10,   1.02501e-15,
     2.24241e+02,   2.26818e-04,   2.95430e-10,   6.51237e-16,
     1.35153e+02,   1.36707e-04,   1.78065e-10,   3.92585e-16,
     7.83236e+01,   7.92263e-05,   1.03204e-10,   2.27687e-16,
     4.41657e+01,   4.46783e-05,   5.82153e-11,   1.28647e-16,
     2.46133e+01,   2.49140e-05,   3.25047e-11,   7.21799e-17  };

  double C20_kpAC[] = 
  {  3.27960e-01,   2.60760e-07,   8.23594e-13,   1.38274e-17,
     4.38752e-01,   3.48855e-07,   1.10197e-12,   1.85031e-17,
     5.78230e-01,   4.59761e-07,   1.45242e-12,   2.43895e-17,
     7.53824e-01,   5.99382e-07,   1.89361e-12,   3.18000e-17,
     1.04013e+00,   8.27053e-07,   2.61333e-12,   4.38935e-17,
     1.41735e+00,   1.12702e-06,   3.56171e-12,   5.98308e-17,
     1.95293e+00,   1.55292e-06,   4.90855e-12,   8.24705e-17,
     2.71532e+00,   2.15922e-06,   6.82677e-12,   1.14729e-16,
     3.79677e+00,   3.01935e-06,   9.54991e-12,   1.60553e-16,
     5.29747e+00,   4.21303e-06,   1.33318e-11,   2.24238e-16,
     7.37841e+00,   5.86846e-06,   1.85820e-11,   3.12737e-16,
     1.02169e+01,   8.12703e-06,   2.57565e-11,   4.33854e-16,
     1.40423e+01,   1.11717e-05,   3.54480e-11,   5.97791e-16,
     1.92026e+01,   1.52804e-05,   4.85668e-11,   8.20344e-16,
     2.61626e+01,   2.08251e-05,   6.63455e-11,   1.12316e-15,
     3.55324e+01,   2.82955e-05,   9.04435e-11,   1.53595e-15,
     4.81644e+01,   3.83784e-05,   1.23253e-10,   2.10252e-15,
     6.53233e+01,   5.20969e-05,   1.68441e-10,   2.89162e-15,
     8.87789e+01,   7.08925e-05,   2.31406e-10,   4.00785e-15,
     1.20729e+02,   9.65780e-05,   3.19528e-10,   5.60280e-15,
     1.63666e+02,   1.31274e-04,   4.42922e-10,   7.90380e-15,
     2.20664e+02,   1.77696e-04,   6.17150e-10,   1.12922e-14,
     2.96276e+02,   2.39957e-04,   8.67756e-10,   1.64168e-14,
     3.97347e+02,   3.24130e-04,   1.23022e-09,   2.41672e-14,
     5.32072e+02,   4.37186e-04,   1.73755e-09,   3.52998e-14,
     7.07827e+02,   5.84876e-04,   2.40229e-09,   4.99094e-14,
     9.32126e+02,   7.72477e-04,   3.21611e-09,   6.73490e-14,
     1.21808e+03,   1.00957e-03,   4.17439e-09,   8.68091e-14,
     1.58941e+03,   1.31461e-03,   5.30260e-09,   1.07984e-13,
     2.08259e+03,   1.71725e-03,   6.66833e-09,   1.31338e-13,
     2.74876e+03,   2.26047e-03,   8.38813e-09,   1.58137e-13,
     3.65886e+03,   3.00609e-03,   1.06436e-08,   1.90525e-13,
     4.91536e+03,   4.04510e-03,   1.37047e-08,   2.31589e-13,
     6.67260e+03,   5.51150e-03,   1.79301e-08,   2.84849e-13,
     9.16963e+03,   7.59478e-03,   2.36873e-08,   3.52355e-13  };

  double C20_kpSiO2D[] = 
  {  7.60358e-02,   1.06666e-07,   3.03247e-13,   1.47482e-18,
     9.07206e-02,   1.27267e-07,   3.61815e-13,   1.75967e-18,
     1.09208e-01,   1.53201e-07,   4.35546e-13,   2.11827e-18,
     1.32481e-01,   1.85851e-07,   5.28369e-13,   2.56972e-18,
     1.58907e-01,   2.22922e-07,   6.33765e-13,   3.08233e-18,
     1.91565e-01,   2.68735e-07,   7.64012e-13,   3.71581e-18,
     2.30490e-01,   3.23342e-07,   9.19258e-13,   4.47087e-18,
     2.76795e-01,   3.88301e-07,   1.10394e-12,   5.36908e-18,
     3.33074e-01,   4.67253e-07,   1.32840e-12,   6.46082e-18,
     4.05326e-01,   5.68613e-07,   1.61658e-12,   7.86245e-18,
     5.08161e-01,   7.12876e-07,   2.02672e-12,   9.85731e-18,
     6.72474e-01,   9.43381e-07,   2.68206e-12,   1.30448e-17,
     9.48552e-01,   1.33068e-06,   3.78319e-12,   1.84008e-17,
     1.41789e+00,   1.98912e-06,   5.65528e-12,   2.75074e-17,
     2.19503e+00,   3.07934e-06,   8.75500e-12,   4.25864e-17,
     3.46724e+00,   4.86415e-06,   1.38300e-11,   6.72794e-17,
     5.76869e+00,   8.09313e-06,   2.30131e-11,   1.11983e-16,
     1.17202e+01,   1.64442e-05,   4.67706e-11,   2.27727e-16,
     3.16486e+01,   4.44090e-05,   1.26336e-10,   6.15454e-16,
     8.68419e+01,   1.21859e-04,   3.46686e-10,   1.68908e-15,
     1.92333e+02,   2.69887e-04,   7.67805e-10,   3.74058e-15,
     3.36287e+02,   4.71894e-04,   1.34247e-09,   6.53939e-15,
     5.05924e+02,   7.09951e-04,   2.01945e-09,   9.83185e-15,
     7.20787e+02,   1.01144e-03,   2.87594e-09,   1.39829e-14,
     9.77602e+02,   1.37169e-03,   3.89813e-09,   1.89179e-14,
     1.18669e+03,   1.66488e-03,   4.72879e-09,   2.29108e-14,
     1.23874e+03,   1.73774e-03,   4.93372e-09,   2.38742e-14,
     1.11204e+03,   1.55988e-03,   4.42750e-09,   2.14067e-14,
     8.76525e+02,   1.22945e-03,   3.48896e-09,   1.68597e-14,
     6.22287e+02,   8.72816e-04,   2.47658e-09,   1.19634e-14,
     4.07334e+02,   5.71310e-04,   1.62094e-09,   7.82836e-15,
     2.50598e+02,   3.51474e-04,   9.97166e-10,   4.81513e-15,
     1.47088e+02,   2.06293e-04,   5.85254e-10,   2.82583e-15,
     8.33219e+01,   1.16862e-04,   3.31541e-10,   1.60077e-15,
     4.59660e+01,   6.44731e-05,   1.82924e-10,   8.83269e-16  };

  double C20_kpMgO[] = 
  {  2.25389e-04,   2.91423e-10,   2.39426e-15,   3.55346e-20,
     4.04967e-04,   5.23622e-10,   4.30206e-15,   6.38511e-20,
     6.31042e-04,   8.15942e-10,   6.70384e-15,   9.94996e-20,
     9.15653e-04,   1.18395e-09,   9.72751e-15,   1.44378e-19,
     1.52197e-03,   1.96795e-09,   1.61693e-14,   2.39993e-19,
     2.37407e-03,   3.06977e-09,   2.52225e-14,   3.74371e-19,
     3.77209e-03,   4.87751e-09,   4.00763e-14,   5.94853e-19,
     6.14348e-03,   7.94395e-09,   6.52733e-14,   9.68874e-19,
     1.01907e-02,   1.31776e-08,   1.08281e-13,   1.60732e-18,
     1.68897e-02,   2.18408e-08,   1.79477e-13,   2.66430e-18,
     2.96125e-02,   3.82967e-08,   3.14746e-13,   4.67293e-18,
     6.10610e-02,   7.89818e-08,   6.49304e-13,   9.64262e-18,
     1.43400e-01,   1.85524e-07,   1.52569e-12,   2.26652e-17,
     3.27383e-01,   4.23645e-07,   3.48510e-12,   5.17901e-17,
     6.39452e-01,   8.27698e-07,   6.81196e-12,   1.01271e-16,
     1.05149e+00,   1.36261e-06,   1.12352e-11,   1.67335e-16,
     1.55882e+00,   2.03505e-06,   1.69775e-11,   2.55740e-16,
     2.94474e+00,   3.91135e-06,   3.34978e-11,   5.17225e-16,
     1.32877e+01,   1.75173e-05,   1.48167e-10,   2.26035e-15,
     7.10035e+01,   9.20831e-05,   7.58876e-10,   1.12902e-14,
     2.54672e+02,   3.27285e-04,   2.65806e-09,   3.89823e-14,
     5.99878e+02,   7.67225e-04,   6.18279e-09,   8.99815e-14,
     9.93398e+02,   1.26701e-03,   1.01648e-08,   1.47279e-13,
     1.24125e+03,   1.58042e-03,   1.26439e-08,   1.82695e-13,
     1.24435e+03,   1.58257e-03,   1.26385e-08,   1.82292e-13,
     1.05321e+03,   1.33853e-03,   1.06770e-08,   1.53820e-13,
     7.84243e+02,   9.96171e-04,   7.93973e-09,   1.14294e-13,
     5.30434e+02,   6.73540e-04,   5.36536e-09,   7.71941e-14,
     3.33897e+02,   4.23882e-04,   3.37536e-09,   4.85451e-14,
     1.99187e+02,   2.52825e-04,   2.01271e-09,   2.89397e-14,
     1.14133e+02,   1.44851e-04,   1.15293e-09,   1.65742e-14,
     6.34329e+01,   8.04978e-05,   6.40626e-10,   9.20833e-15,
     3.44484e+01,   4.37133e-05,   3.47851e-10,   4.99954e-15,
     1.83799e+01,   2.33224e-05,   1.85579e-10,   2.66709e-15,
     9.67434e+00,   1.22756e-05,   9.76752e-11,   1.40372e-15  };

  double C20_kpFeS[] = 
  {  5.18099e-02,   8.75057e-08,   1.64042e-13,   3.48481e-19,
     9.98914e-02,   1.68714e-07,   3.16278e-13,   6.71882e-19,
     1.60422e-01,   2.70949e-07,   5.07932e-13,   1.07902e-18,
     2.36626e-01,   3.99656e-07,   7.49210e-13,   1.59157e-18,
     3.67294e-01,   6.20351e-07,   1.16293e-12,   2.47046e-18,
     5.36239e-01,   9.05695e-07,   1.69785e-12,   3.60681e-18,
     7.64225e-01,   1.29076e-06,   2.41970e-12,   5.14028e-18,
     1.04974e+00,   1.77298e-06,   3.32371e-12,   7.06072e-18,
     1.38085e+00,   2.33223e-06,   4.37210e-12,   9.28786e-18,
     1.74382e+00,   2.94529e-06,   5.52140e-12,   1.17294e-17,
     2.10325e+00,   3.55236e-06,   6.65945e-12,   1.41471e-17,
     2.42167e+00,   4.09020e-06,   7.66778e-12,   1.62893e-17,
     2.66848e+00,   4.50708e-06,   8.44938e-12,   1.79499e-17,
     2.81672e+00,   4.75754e-06,   8.91905e-12,   1.89480e-17,
     2.89933e+00,   4.89720e-06,   9.18118e-12,   1.95057e-17,
     3.08200e+00,   5.20610e-06,   9.76116e-12,   2.07399e-17,
     3.69875e+00,   6.24871e-06,   1.17179e-11,   2.49020e-17,
     5.06356e+00,   8.55561e-06,   1.60466e-11,   3.41079e-17,
     7.08905e+00,   1.19792e-05,   2.24709e-11,   4.77704e-17,
     9.21663e+00,   1.55761e-05,   2.92218e-11,   6.21311e-17,
     1.08429e+01,   1.83272e-05,   3.43894e-11,   7.31341e-17,
     1.17217e+01,   1.98178e-05,   3.71989e-11,   7.91392e-17,
     1.19809e+01,   2.02656e-05,   3.80618e-11,   8.10292e-17,
     1.18982e+01,   2.01424e-05,   3.78698e-11,   8.07160e-17,
     1.17194e+01,   1.98680e-05,   3.74208e-11,   7.99219e-17,
     1.16101e+01,   1.97358e-05,   3.72971e-11,   7.99631e-17,
     1.19784e+01,   2.05514e-05,   3.92888e-11,   8.53362e-17,
     1.49136e+01,   2.65380e-05,   5.30405e-11,   1.20988e-16,
     2.57980e+01,   4.84122e-05,   1.02966e-10,   2.50769e-16,
     6.12919e+01,   1.18083e-04,   2.58706e-10,   6.49359e-16,
     1.55902e+02,   3.00837e-04,   6.60279e-10,   1.66030e-15,
     3.47036e+02,   6.65409e-04,   1.45009e-09,   3.62025e-15,
     6.57260e+02,   1.25056e-03,   2.70176e-09,   6.68619e-15,
     1.09327e+03,   2.06455e-03,   4.42254e-09,   1.08501e-14,
     1.65394e+03,   3.10062e-03,   6.58692e-09,   1.60229e-14  };

  double C20_kpAl2O3[] = 
  {  9.93250e-04,   9.14846e-11,   8.97410e-18,   9.30612e-25,
     1.81240e-03,   1.66933e-10,   1.63752e-17,   1.69810e-24,
     2.84365e-03,   2.61918e-10,   2.56926e-17,   2.66432e-24,
     4.14191e-03,   3.81496e-10,   3.74225e-17,   3.88071e-24,
     7.18271e-03,   6.61573e-10,   6.48964e-17,   6.72974e-24,
     1.13364e-02,   1.04415e-09,   1.02425e-16,   1.06215e-23,
     1.77361e-02,   1.63360e-09,   1.60247e-16,   1.66176e-23,
     2.59477e-02,   2.38995e-09,   2.34440e-16,   2.43114e-23,
     3.45425e-02,   3.18159e-09,   3.12095e-16,   3.23642e-23,
     4.22006e-02,   3.88695e-09,   3.81286e-16,   3.95393e-23,
     4.71420e-02,   4.34208e-09,   4.25932e-16,   4.41691e-23,
     4.91934e-02,   4.53102e-09,   4.44466e-16,   4.60911e-23,
     5.05162e-02,   4.65286e-09,   4.56418e-16,   4.73304e-23,
     5.78201e-02,   5.32560e-09,   5.22410e-16,   5.41738e-23,
     8.84237e-02,   8.14438e-09,   7.98916e-16,   8.28474e-23,
     1.78786e-01,   1.64673e-08,   1.61535e-15,   1.67511e-22,
     4.36404e-01,   4.01956e-08,   3.94295e-15,   4.08884e-22,
     1.63796e+00,   1.50867e-07,   1.47992e-14,   1.53467e-21,
     8.50819e+00,   7.83659e-07,   7.68723e-14,   7.97165e-21,
     3.92751e+01,   3.61749e-06,   3.54854e-13,   3.67984e-20,
     1.41439e+02,   1.30275e-05,   1.27792e-12,   1.32520e-19,
     3.83709e+02,   3.53420e-05,   3.46684e-12,   3.59511e-19,
     7.70411e+02,   7.09598e-05,   6.96073e-12,   7.21827e-19,
     1.16399e+03,   1.07211e-04,   1.05167e-11,   1.09058e-18,
     1.37566e+03,   1.26707e-04,   1.24292e-11,   1.28891e-18,
     1.33070e+03,   1.22566e-04,   1.20230e-11,   1.24678e-18,
     1.09978e+03,   1.01297e-04,   9.93663e-12,   1.03043e-18,
     8.05638e+02,   7.42044e-05,   7.27901e-12,   7.54832e-19,
     5.38690e+02,   4.96167e-05,   4.86711e-12,   5.04718e-19,
     3.36338e+02,   3.09789e-05,   3.03884e-12,   3.15127e-19,
     1.99460e+02,   1.83715e-05,   1.80214e-12,   1.86881e-19,
     1.13787e+02,   1.04805e-05,   1.02808e-12,   1.06611e-19,
     6.30411e+01,   5.80648e-06,   5.69582e-13,   5.90655e-20,
     3.41529e+01,   3.14570e-06,   3.08575e-13,   3.19991e-20,
     1.81893e+01,   1.67535e-06,   1.64342e-13,   1.70422e-20  };


  itab0 = Nmom * NTd * iSN;
  itab  = 0;
  for(imom = 0; imom < Nmom; imom++) {
    for(iTd = 0; iTd < NTd; iTd++) {
      my_rates->SN0_kpSiM     [itab0] = C20_kpSiM     [itab];
      my_rates->SN0_kpFeM     [itab0] = C20_kpFeM     [itab];
      my_rates->SN0_kpMg2SiO4 [itab0] = C20_kpMg2SiO4 [itab];
      my_rates->SN0_kpMgSiO3  [itab0] = C20_kpMgSiO3  [itab];
      my_rates->SN0_kpAC      [itab0] = C20_kpAC      [itab];
      my_rates->SN0_kpSiO2D   [itab0] = C20_kpSiO2D   [itab];
      my_rates->SN0_kpMgO     [itab0] = C20_kpMgO     [itab];
      my_rates->SN0_kpFeS     [itab0] = C20_kpFeS     [itab];
      my_rates->SN0_kpAl2O3   [itab0] = C20_kpAl2O3   [itab];
      itab0++;
      itab ++;
    }
  }

  return SUCCESS;
}
