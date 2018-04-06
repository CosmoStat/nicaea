/* ==================================================== *
 * fr_constants.h                                      	*
 * Martin Kilbinger 2013			   	*
 * Constants from .h files in FrankenEmu.           	*
 * ==================================================== */

#ifndef __FR_CONSTANTS_H
#define __FR_CONSTANTS_H

#include "coyote.h"

/* These constants have the same name but different *
 * values compared to coyote v<2.                   */

#define m 137
#define neta 5500
#define p 6
#define peta 9
#define rs fr_rs

/* neta_over_rs = neta / rs */
# define neta_over_rs 500

const double fr_ksim[fr_nsim];

const double fr_kemu[neta_over_rs]; /* MKDEBUG: Size was 1000 in original Coyote emulator, for some reason */
const double fr_mean[neta];
const double fr_K[neta][peta];
const double fr_x[m][p];
const double fr_xmin[p];
const double fr_xrange[p];
const double fr_aemu[rs];
const double fr_lamz[peta];
const double fr_beta[peta][p];
const double fr_KrigBasis[peta][m];

/* The above defined quantities are undefined at the end of emu.c */

#endif

