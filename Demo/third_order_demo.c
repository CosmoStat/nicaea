/* ============================================================ *
 * third_order_demo.c						*
 * Martin Kilbinger 2010					*
 * ============================================================ */

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <time.h>

#include "errorlist.h"
#include "io.h"
#include "maths.h"
#include "cosmo.h"
#include "cmb_bao.h"
#include "lensing.h"
#include "lensing_3rd.h"
#include "nofz.h"


#define NTHETA 20
#define TH_MIN (1.0 / sqrt(2.0))
#define TH_FAC (sqrt(2.0))
void usage(int ex)
{
   fprintf(stderr, "Usage: lensingdemo [OPTIONS]\n");
   fprintf(stderr, "OPTIONS\n");
   fprintf(stderr, "  -t FNAME      File containing angular scales in arcmin\n");
   fprintf(stderr, "                 (default: %d bins between %g and %g arcmin)\n",
	   NTHETA, TH_MIN, TH_MIN * pow(TH_FAC, NTHETA-1));
   fprintf(stderr, "  -d            Diagonal only (three equal angular scales\n");
   fprintf(stderr, "  -h            This message\n");

   if (ex >= 0) exit(ex);
}

int SLC_test(cosmo_3rd *self, error **err)
{

   double theta, SLC;

   theta = 2.0 * arcmin;

   SLC = map3_SLC_t1(self, theta, 0, err);
   forwardError(*err, __LINE__, 0.0);
   printf("SLC = %g\n", SLC);
   return 0;


   double a1, a2, q, abserr, chisqr, w1, w2, f1, f2;

   a1    = 0.4; 
   a2    = 0.8;

   w1  = w(self->lens->cosmo, a1, 0, err);       forwardError(*err, __LINE__, 0.0);
   f1  = f_K(self->lens->cosmo, w1, err);      	 forwardError(*err, __LINE__, 0.0);
   w2  = w(self->lens->cosmo, a2, 0, err);       forwardError(*err, __LINE__, 0.0);
   f2 = f_K(self->lens->cosmo, w2, err);         forwardError(*err, __LINE__, 0.0);

   q     = Q_mc(self, a1, a2, f1, f2, theta, &abserr, &chisqr, err);
   forwardError(*err, __LINE__, 0.0);

   printf("q, relerr, chisqr = %g %g %g\n", q, abserr/q, chisqr);



   return 0;
}

int main(int argc, char** argv)
{
   error *myerr = NULL, **err;
   FILE *F;
   //const double THETA[NTHETA] = {1.0, 5.4772, 30.0};
   double *THETA;
   int i, j, k, c, diag_only;
   size_t Nrec, Ntheta;
   double theta[3], res;
   cosmo_3rd *model;
   char *theta_fname;
   clock_t c_start;

   err = &myerr;

   theta_fname = NULL;
   diag_only   = 0;
   while (1) {

      static struct option long_option[] = {
	{"", required_argument, 0, 't'},
	{0, 0, 0, 0}
      };

      int option_index = 0;

      c = getopt_long(argc, argv, "dt:h", long_option, &option_index);
      switch (c) {
	 case 't' :
	    theta_fname = optarg;
	    break;
	 case 'd' :
	    diag_only = 1;
	    break;
	 case 'h' :
	    usage(0);
	 case -1 :
	    goto end_arg;
	 default :
	    usage(1);
      }

   }
 end_arg:

   F = fopen_err("cosmo_3rd.par", "r", err);
   quitOnError(*err, __LINE__, stderr);
   read_cosmological_parameters_lens_3rd(&model, F, err);
   quitOnError(*err, __LINE__, stderr);
   fclose(F);

   printf("# Parameters:\n");
   dump_param_3rd(model, stdout, err);
   quitOnError(*err, __LINE__, stderr);

   c_start = clock();


   if (theta_fname != NULL) {
      Nrec = 0;
      THETA = (double*)read_any_list_count(theta_fname, &Nrec, "%lg", sizeof(double), &Ntheta, err);
      quitOnError(*err, __LINE__, stderr);
   } else {
      Ntheta = NTHETA;
      THETA = malloc_err(sizeof(double) * Ntheta, err);
      quitOnError(*err, __LINE__, stderr);
      for (i=0; i<Ntheta; i++) {
	 THETA[i] = TH_MIN * pow(TH_FAC, i);
      }
   }

   for (i=0; i<Ntheta; i++) {
      theta[0] = THETA[i]*arcmin;
      for (j=i;j<Ntheta; j++) {
	 theta[1]= THETA[j]*arcmin;
	 for (k=j; k<Ntheta; k++) {
	    theta[2] = THETA[k]*arcmin;

	    if (diag_only && (i != j || j != k)) continue;

	    res = map3(model, theta, 0, 0, 0, fgauss, err);
	    quitOnError(*err, __LINE__, stderr);

	    printf("%8.3f %8.3f %8.3f  % g\n", theta[0]/arcmin, theta[1]/arcmin, theta[2]/arcmin, res);
	    fflush(stdout);
	 }
      }
   }

   free(THETA);

   end_clock(c_start, stderr);

   return 0;
}
#undef NTHETA
#undef TH_MIN
#undef TH_FAC
