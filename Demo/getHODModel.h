/* ------------------------------------------------------------ *
 * getHODModel.h                                                *
 * Jean Coupon 2015                                             *
 * ------------------------------------------------------------ */

#include "hod.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_histogram.h>

#define COSMO     0
#define HMF       1
#define DM        2
#define SMF       3
#define WTHETA    4
#define WP        5
#define XI        6
#define GGL       7
#define MSMH      8
#define SHMR      9
#define TESTS     10

/* ---------------------------------------------------------------- *
 * Variables
 * ---------------------------------------------------------------- */

char MYNAME[100];

/* ---------------------------------------------------------------- *
 * Structures
 * ---------------------------------------------------------------- */

typedef struct ConfigGetModel
{
  /* Input and output files */
  char fileOutName[1000], configpmcFileInName[1000], sampleFileInName[1000],
    configmcmcFileInName[1000];
  
  /* Cosmological model */
  cosmo_hm *model;
  
  /* mean redshift, a, range, Halo mass */
  double z, a, min, max, delta, log10Mhalo;
  
  /* output type, nbins and binsize*/
  int type, nbins;
  
} ConfigGetModel;

/* ---------------------------------------------------------------- *
 * Main fonctions
 * ---------------------------------------------------------------- */

/* Dark matter functions */
void writeCosmo(ConfigGetModel para, error **err);
void writeHMF(ConfigGetModel para, error **err);
void writeTPCF_DM(ConfigGetModel para, error **err);

/* Galaxies functions */
void writeSMF(ConfigGetModel para, error **err);
void writeTPCF_GAL(ConfigGetModel para, error **err);
void writeGGL(ConfigGetModel para, error **err);

/* SHMR functions */
void writeMSMH(ConfigGetModel para, error **err);
void writeSHMR(ConfigGetModel para, error **err);

/* Others */
void writeTESTS(ConfigGetModel para, error **err);



/*----------------------------------------------------------------*
 * Utils
 *----------------------------------------------------------------*/

void initPara(int argc, char **argv, ConfigGetModel *para, error **err);
long fileSize(FILE *fileIn);

