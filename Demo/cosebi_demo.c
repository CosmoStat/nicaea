/* ============================================================ *
 * cosebi_demo.c						*
 * Martin Kilbinger, Liping Fu 2012.				*
 *  2013. Tomography.						*
 * Reference:							*
 * Schneider, Eifler, Krause 2010                               * 
 * ============================================================ */

#include <stdio.h>
#include <stdlib.h>
#include "errorlist.h"
#include "io.h"
#include "maths.h"
#include "cosmo.h"
#include "lensing.h"
#include "nofz.h"
#include "decomp_eb.h"

#include <getopt.h>
#include <gsl/gsl_sf_gamma.h>


void write_tpmlog(int n, double thmin, double thmax, const double *c, error **err)
{
   FILE *F;
   double z, theta, tp, tm, dth;
   char name[128];

   sprintf(name, "T_%d_%.3f_%.3f.dat", n, thmin/arcmin, thmax/arcmin);

   F = fopen_err(name, "w", err); forwardError(*err, __LINE__,);

   fprintf(F, "# z=log(thmax/thmin) theta[arcmin] T_+^log T_-^log, (n, thmin, thmax = %d, %g, %g)\n",
	   n, thmin/arcmin, thmax/arcmin);

   dth = 1.01;
   for (theta=thmin; theta<=thmax; theta*=dth) {

      z = log(theta / thmin);

      tp = Tplog_c(z, c, n, err);
      forwardError(*err, __LINE__,);

      tm = Tmlog(z, c, n, err);
      forwardError(*err, __LINE__,);

      fprintf(F, "%f %g % g % g\n", z, theta/arcmin, tp, tm);
   }

   fclose(F);
}


void write_xi_pm_EB(int N, double thmin, double thmax, const double *c_cosebi,
		    const double *E_cos, const double *B_cos, error **err)
{
   double xi_pm_EB[4], theta, fth;
   const char name[128] = "xi_pm_EB.dat";
   int n;
   FILE *F;

   F = fopen_err(name, "w", err);       forwardError(*err, __LINE__,);
   fprintf(F, "# theta[arcmin] E+ E- B+ B-\n");

   fth = 1.1;
   for (theta=thmin; theta<=thmax; theta *= fth) {
      xipmEB(theta, thmin, thmax, c_cosebi, E_cos, B_cos, N, xi_pm_EB, err);
      forwardError(*err, __LINE__,);

      fprintf(F, "%12.3g  ", theta/arcmin);
      for (n=0; n<4; n++) fprintf(F, " % g", xi_pm_EB[n]);
      fprintf(F, "\n");
   }

   fclose(F);
}

double E_cosebi_data(const double *xip, const double *xim, const double *theta, int Nxi,
      int n, double Psimin, double Psimax, const double *c_cosebi,
      double *B_cosebi, error **err)
{
   double z, delta_theta, rp, rm, summand;
   int i;


   testErrorRetVA(Psimin < theta[0], mr_range, "COSEBI theta_min = %g' is smaller than xi theta_min = %g'",
         *err, __LINE__, -1.0, Psimin / arcmin, theta[0] / arcmin);
   testErrorRetVA(Psimax > theta[Nxi-1], mr_range, "COSEBI theta_max = %g' is larger than xi theta_max = %g'",
         *err, __LINE__, -1.0, Psimax / arcmin, theta[Nxi-1] / arcmin);

   rp = rm = 0.0;
   for (i=0; i<Nxi; i++) {

      if (theta[i] < Psimin) continue;
      if (theta[i] > Psimax) break;

      if (i > 0) {
         delta_theta = theta[i] - theta[i-1];
      } else {
         delta_theta = theta[1] - theta[0];
      }

      z = log(theta[i] / Psimin);

      summand  = xip[i];
      summand *= theta[i] * delta_theta;
      summand *= Tplog_c(z, c_cosebi, n, err);
      forwardError(*err, __LINE__, -1.0);
      rp      += summand;

      summand  = xim[i];
      summand *= theta[i] * delta_theta;
      summand *= Tmlog(z, c_cosebi, n, err);
      forwardError(*err, __LINE__, -1.0);
      rm      += summand;

      //printf("%g %g   %g %g\n", theta[i]/arcmin, z, xip[i], xim[i]);

   }

   if (B_cosebi != NULL) {
      *B_cosebi = 0.5 * (rp - rm);
   }

   return 0.5 * (rp + rm);
}

void usage(int ex)
{
   fprintf(stderr, "Calculates the COSEBIs from Schneider, Eifler, Krause 2010\n");
   fprintf(stderr, "Usage: coseb_demo [OPTIONS]\n");
   fprintf(stderr, "OPTIONS\n");
   fprintf(stderr, "  -o OUT     Output name OUT (default: stdout)\n");
   fprintf(stderr, "  -N N       Maximum COSEBIs mode N (default 10)\n");
   fprintf(stderr, "  -m THMIN   Minimum COSEBIs scale THMIN (in arcmin; default 1)\n");
   fprintf(stderr, "  -M THMAX   Maximum COSEBIs scale THMAX (in arcmin; default 250)\n");
   fprintf(stderr, "  -d XI      Data file XI (cosmo_pmc input ['xipm'] format; default none)\n");
   fprintf(stderr, "  -T         Write filter functions to files 'T_n_thmin_thmax.dat'\n");
   fprintf(stderr, "  -X         Write EB correlation functions to files 'xEB_n_thmin_thmax.dat'\n");
   fprintf(stderr, "  -P PATH    Path for coefficient files ('cosebi_tplog_rN_N_thmin_thmax';\n");
   fprintf(stderr, "              default $COSMOPMC/par_files/COSEBIs)\n");
   fprintf(stderr, "  -q         Quiet\n");
   fprintf(stderr, "  -h         This message\n");

   if (ex >= 0) exit(ex);
}


int main(int argc, char* argv[])
{
   error *myerr = NULL, **err;
   double *E_cos, *B_cos;
   double Psimin, Psimax;
   int c, N, n, Nxi, write_T, write_X, i_bin, j_bin, Nzbin, quiet;
   cosmo_lens *model;
   char *xi_name, *path, *strp, *output_name;
   FILE *F, *Fout;
   datcov *dc;
   double *xip, *xim, *theta, *theta2, *c_cosebi;


   err = &myerr;


   /* Command line arguments */
   Psimin  = 1.0;
   Psimax  = 250.0;
   N       = 10;
   xi_name = NULL;
   write_T = write_X = 0;
   path    = NULL;
   output_name = NULL;
   quiet   = 0;
   while (1) {

      static struct option long_options[] = {
         {"", required_argument, 0, 'N'},
         {"", required_argument, 0, 'd'},
         {"", required_argument, 0, 'P'},
         {"", required_argument, 0, 'o'},
         {0, 0, 0, 0}
      };

      int option_index = 0;

      c = getopt_long(argc, argv, "o:N:d:m:M:TXP:qh", long_options, &option_index);
      switch (c) {
         case 'N' :
            N = atoi(optarg);
            break;
         case 'd' :
            xi_name = optarg;
            break;
         case 'P' :
            path = optarg;
            break;
         case 'o' :
            output_name = optarg;
            break;
         case 'm' :
            Psimin = atof(optarg);
            break;
         case 'M' :
            Psimax = atof(optarg);
            break;
         case 'T' :
            write_T = 1;
            break;
         case 'X' :
            write_X = 1;
            break;
         case 'q' :
            quiet = 1;
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


   if (quiet == 0) {
      fprintf(stderr, "# cosebi_demo (M. Kilbinger, L. Fu, 2012)\n");
   }


   if (N<0 || N>NMAX_COSEBI) {
      fprintf(stderr, "N=%d out of range [0; %d]", N, NMAX_COSEBI);
      exit(4);
   }

   if (path == NULL) {
      strp = getenv("COSMOPMC");
      if (strp == NULL) {
         strp = ".";
      }
      path = malloc_err(sizeof(char) * 1024, err);
      quitOnError(*err, __LINE__, stderr);
      sprintf(path, "%s/par_files/COSEBIs", strp);
   }

   if (output_name == NULL) {
      Fout = stdout;
   } else {
      Fout = fopen_err(output_name, "w", err);
      quitOnError(*err, __LINE__, stderr);
   }

   E_cos = malloc_err(sizeof(double) * (N+1), err);      quitOnError(*err, __LINE__, stderr);
   B_cos = malloc_err(sizeof(double) * (N+1), err);      quitOnError(*err, __LINE__, stderr);
   
   if (xi_name == NULL) {

      /* COSEBIs from theoretical model */

      F = fopen_err("cosmo_lens.par", "r", err);
      quitOnError(*err, __LINE__, stderr);
      read_cosmological_parameters_lens(&model, F, err);
      quitOnError(*err, __LINE__, stderr);
      fclose(F);
      Nzbin = model->redshift->Nzbin;

      dump_param(model->cosmo, stdout);
      //fprintf(Fout, "# n E_n^{ij} B_n^{ij}\n# n");
      for (i_bin=0; i_bin<Nzbin; i_bin++) {
         for (j_bin=i_bin; j_bin<Nzbin; j_bin++) {
            if (model->tomo==tomo_auto_only && i_bin!=j_bin) continue;
            if (model->tomo==tomo_cross_only && i_bin==j_bin) continue;
            fprintf(Fout, "        E_n^{%d%d}        B_n^{%d%d}", i_bin, j_bin, i_bin, j_bin);
         }
      }
      fprintf(Fout, "\n");

      for (n=1; n<=N; n++) {
         fprintf(Fout, "%3d", n);

         for (i_bin=0; i_bin<Nzbin; i_bin++) {
            for (j_bin=i_bin; j_bin<Nzbin; j_bin++) {
               if (model->tomo==tomo_auto_only && i_bin!=j_bin) continue;
               if (model->tomo==tomo_cross_only && i_bin==j_bin) continue;
               E_cos[n] = E_cosebi(model, n, Psimin*arcmin, Psimax*arcmin, i_bin, j_bin, path, B_cos+n, err);
               quitOnError(*err, __LINE__, stderr);
               fprintf(Fout, " % 15.6g % 15.6g", E_cos[n], B_cos[n]);
            }
         }

         fprintf(Fout, "\n");
      }

      c_cosebi = model->c_cosebi;

   } else {

      /* COSEBIs from data */
      c_cosebi = read_zeros_norm_cosebi_auto_check(Psimin*arcmin, Psimax*arcmin, path, err);
      quitOnError(*err, __LINE__, stderr);

      dc = malloc_err(sizeof(datcov), err);
      quitOnError(*err, __LINE__, stderr);
      dc->format = angle_center;
      /* Reading xi file */
      read_data_tomo(dc, xi_name, 0, second_order, err);
      quitOnError(*err, __LINE__, stderr);

      Nzbin    = dc->Nzbin;
      dc->type = xipm;

      /* Header */
      //fprintf(Fout, "# n E_n^{ij} B_n^{ij}\n");
      fprintf(Fout, "# n   ");
      for (i_bin=0; i_bin<Nzbin; i_bin++) {
         for (j_bin=i_bin; j_bin<Nzbin; j_bin++) {
            fprintf(Fout, "E_n^{%d%d}        B_n^{%d%d}          ", i_bin, j_bin, i_bin, j_bin);
         }
      }
      fprintf(Fout, "\n");

      for (n=1; n<=N; n++) {
         fprintf(Fout, "%3d", n);

         for (i_bin=0; i_bin<Nzbin; i_bin++) {
            for (j_bin=i_bin; j_bin<Nzbin; j_bin++) {
               datcov2xipm(dc, i_bin, j_bin, &xip, &xim, &theta, &theta2, &Nxi, err);
               quitOnError(*err, __LINE__, stderr);

               E_cos[n] = E_cosebi_data(xip, xim, theta, Nxi, n, Psimin*arcmin, Psimax*arcmin, c_cosebi, B_cos+n, err);
               quitOnError(*err, __LINE__, stderr);
               fprintf(Fout, "   % .9g % .9g", E_cos[n], B_cos[n]);
            }
         }
         fprintf(Fout, "\n");

         free(xip);
         free(xim);
         free(theta);
      }

      del_data_cov(&dc);

   }

   if (output_name != NULL) {
      fclose(Fout);
   }


   if (write_T) {

      /* Write filter functions to file */

      for (n=1; n<=N; n++) {
         write_tpmlog(n, Psimin*arcmin, Psimax*arcmin, c_cosebi, err);
         quitOnError(*err, __LINE__, stderr);
      }

   }

   if (write_X) {

      /* Calculate xipmEB according to SEK 10 eq (40) */

      write_xi_pm_EB(N, Psimin*arcmin, Psimax*arcmin, c_cosebi, E_cos, B_cos, err);

   }

   free(E_cos);
   free(B_cos);

   return 0;
}
