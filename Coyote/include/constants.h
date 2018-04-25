/* ============================================================ *
 * constants.h							*
 * Definition of constants for the Coyote emulator of the non-  *
 * linear power spectrum.					*
 * ============================================================ */

#ifndef __CONSTANTS_H
#define __CONSTANTS_H


#define m    37
#define neta 6000
#define p    5
#define peta 5
#define rs   6

const double ksim[nsim];

const double kemu[1000];
const double mean[neta];
const double K[neta][peta];
const double x[m][p], xmin[p], xrange[p], aemu[rs];
const double lamws[peta], lamz[peta];
const double beta[peta][p];
const double w_coyote[peta][m];
const double KrigBasis[peta][m];

/* The above defined quantities are undefined at the end of emu.c */

#endif


