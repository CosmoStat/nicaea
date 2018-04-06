/* ============================================================ *
 * decomp_eb_demo.c						*
 * Martin Kilbinger, Liping Fu 2009				*
 * ============================================================ */

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

#include "errorlist.h"
#include "io.h"
#include "maths.h"
#include "cosmo.h"
#include "lensing.h"
#include "nofz.h"
#include "decomp_eb.h"

void usage(int ex)
{
   fprintf(stderr, "Usage: decomp_eb_decom [OPTIONS]\n");
   fprintf(stderr, "OPTIONS:\n");
   fprintf(stderr, "  -c IND        Coefficients for E-/B-mode estimator, optimised according to:\n");
   fprintf(stderr, "                 0: S/N, eta=1/50\n");
   fprintf(stderr, "                 1: FoM, eta=1/10\n");
   fprintf(stderr, "                 2: FoM, eta=1/50\n");
   fprintf(stderr, "  -x FNAME      2PCF from file FNAME (xipm format; default: 2PCF from theory)\n");
   fprintf(stderr, "  -T 'THMIN THMAX NTH'\n");
   fprintf(stderr, "                Calculate E-/B-modes on NTH logarithmis scales between THMIN\n");
   fprintf(stderr, "                 and THMAX (in arcmin); default: '1 100 15'\n");
   fprintf(stderr, "  -o FNAME      Ouput file name FNAME (default: 'REB')\n");
   fprintf(stderr, "  -f            Calculate and write filter functions T+, T-\n");
   fprintf(stderr, "  -h            This message\n");

   exit(ex);
}

#define N 6
int main(int argc, char** argv)
{
   error *myerr = NULL, **err;
   double eta, R, Psi, rp, rm, thmin, thmax, fth;
   const double *a;
   char choice;
   cosmo_lens *model;
   FILE *F, *Fc;
   char *xi_fname, *out_fname;
   int c, do_filter, Nth;


   err = &myerr;

   printf("decomp_eb_demo (M. Kilbinger, L. Fu 2009)\n");


   xi_fname  = NULL;
   out_fname = NULL;
   choice    = '?';
   thmin     = 1.0;
   thmax     = 100.0;
   Nth       = 15;
   do_filter = 0;
   while (1) {

      static struct option long_option[] = {
	{"", required_argument, 0, 'c'},
	{"", required_argument, 0, 'x'},
	{"", required_argument, 0, 'T'},
	{"", required_argument, 0, 'o'},
	{0, 0, 0, 0}
      };

      int option_index = 0;

      c = getopt_long(argc, argv, "c:x:T:o:fh", long_option, &option_index);
      switch (c) {
	 case 'x' :
	    xi_fname = optarg;
	    break;
	 case 'c' :
	    choice = optarg[0];
	    break;
	 case 'T' :
	    sscanf(optarg, "%lf %lf %d\n", &thmin, &thmax, &Nth);
	    break;
	 case 'o' :
	    out_fname = optarg;
	    break;
	 case 'f' :
	    do_filter = 1;
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

   fth = pow(thmax/thmin, 1.0/((double)Nth - 0.999));

   if (choice == '?') {
      printf("Choice of coefficients for filter function:\n");
      printf("  Covariance from CFHTLS Wide 3rd data release (Fu et al. 2008)\n");
      printf("     S/N, eta=1/50          [0]\n");
      printf("     FoM, eta=1/10          [1]\n");
      printf("     FoM, eta=1/50          [2]\n");
      printf("Choice? ");
      choice = getchar();
   }

   switch (choice) {
      case '0' : a = a_FK10_SN;        eta = eta_FK10_SN; break;
      case '1' : a = a_FK10_FoM_eta10; eta = eta_FK10_FoM_eta10; break;
      case '2' : a = a_FK10_FoM_eta50; eta = eta_FK10_FoM_eta50; break;
      default  : assert(0);
   }

   R = (1.0+eta)/(1.0-eta);    /* FK09 (10) */

   printf("(%s)\n", out_fname);
   if (out_fname == NULL) {
      out_fname = "REB";
   }
   printf("(%s)\n", out_fname);


   if (do_filter == 1) {
      double tp, tm, x;
      printf("Writing filter functions T+, T- to file 'Tpm'\n");
      F = fopen_err("Tpm", "w", err);
      quitOnError(*err, __LINE__, stderr);
      for (x=-1.0; x<=1.0005; x+=0.005) {
	 tp = Tp(x, a, N, cheby2, err);
	 quitOnError(*err, __LINE__, stderr);
	 tm = Tm(x, a, N, cheby2, R, err);
	 quitOnError(*err, __LINE__, stderr);
	 fprintf(F, "% f % f % f\n", x, tp, tm);
      }
      fclose(F);

      return 0;
   }


   printf("Writing shear functions R_E, R_B to file '%s'\n", out_fname);
   F = fopen_err(out_fname, "w", err);
   quitOnError(*err, __LINE__, stderr);


   if (xi_fname == NULL) {

      printf("Optimised E-/B-mode function from theory\n");

      Fc = fopen_err("cosmo_lens.par", "r", err);
      quitOnError(*err, __LINE__, stderr);
      read_cosmological_parameters_lens(&model, Fc, err);
      quitOnError(*err, __LINE__, stderr);
      fclose(Fc);

      for (Psi=thmin*arcmin; Psi<=thmax*arcmin; Psi*=fth) {
	 rp = RR(model, Psi*eta, Psi, a, N, cheby2, +1, err);
	 quitOnError(*err, __LINE__, stderr);
	 rm = RR(model, Psi*eta, Psi, a, N, cheby2, -1, err);
	 quitOnError(*err, __LINE__, stderr);
	 fprintf(F, "% f % g % g\n", Psi/arcmin, 0.5*fabs(rp+rm), 0.5*fabs(rp-rm));
	 fflush(F);
      }

   } else {

      datcov *dc;
      int i_bin=0, j_bin=0, Nxi, Nzbin;
      double *xip, *xim, *theta, *theta2;

      printf("Optimised E-/B-mode function from data (file %s)\n", xi_fname);

      dc = malloc_err(sizeof(datcov), err);
      quitOnError(*err, __LINE__, stderr);
      dc->format = angle_center;
      read_data_tomo(dc, xi_fname, 0, second_order, err);
      quitOnError(*err, __LINE__, stderr);

      Nzbin    = dc->Nzbin;
      dc->type = xipm;

      /* Header */
      fprintf(F, "# Psi[arcmin]");	    fprintf(F, "(RE RB)^{mn}\n  ");
      for (i_bin=0; i_bin<Nzbin; i_bin++) {
	 for (j_bin=i_bin; j_bin<Nzbin; j_bin++) {
	    fprintf(F, "          %d%d ", i_bin, j_bin);
	    fprintf(F, "              ");
	 }
      }
      fprintf(F, "\n");

      for (Psi=thmin*arcmin; Psi<=thmax*arcmin; Psi*=fth) {
	 fprintf(F, "% f", Psi/arcmin);
	 for (i_bin=0; i_bin<Nzbin; i_bin++) {
	    for (j_bin=i_bin; j_bin<Nzbin; j_bin++) {
	       datcov2xipm(dc,  i_bin,  j_bin,  &xip,  &xim,  &theta, &theta2,  &Nxi,  err);
	       quitOnError(*err, __LINE__, stderr);

	       rp = RR_data(xip, xim, theta, Nxi, Psi*eta, Psi, a, N, cheby2, +1, err);
	       quitOnError(*err, __LINE__, stderr);
	       rm = RR_data(xip, xim, theta, Nxi, Psi*eta, Psi, a, N, cheby2, -1, err);
	       quitOnError(*err, __LINE__, stderr);

	       free(xip); free(xim); free(theta);

	       fprintf(F, "  % g % g", 0.5*fabs(rp+rm), 0.5*fabs(rp-rm));
	    }
	 }
	 fprintf(F, "\n");
      }

      del_data_cov(&dc);
   }

   fclose(F);

   return 0;
}
