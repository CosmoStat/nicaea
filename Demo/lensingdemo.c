/* ============================================================ *
 * lensingdemo.c						                               *
 * Martin Kilbinger 2008-2010					                      *
 * ============================================================ */

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>

#include "errorlist.h"
#include "io.h"
#include "maths.h"
#include "cosmo.h"
#include "coyote.h"
#include "cmb_bao.h"
#include "lensing.h"
#include "nofz.h"


#define TH_MIN  0.5
#define TH_MAX  400.0
#define NTHETA  71
#define ELL_MIN 0.01
#define ELL_MAX 10000.0
#define NELL    50

double mu_a_k(cosmo *cosmo, double a, double M_0, double sigma, double n, double k)
{
   double mu, M;

   // G.-B. Zhao et al. (2012)
   //M = M_0 * pow(a, -sigma);

   // Harnois-Deraps et al. (2015)
   M = M_0 * pow((4*cosmo->Omega_de + cosmo->Omega_m*pow(a, -3))/(4*cosmo->Omega_de + cosmo->Omega_m), (n+2.0)/2.0);

   mu = 1.0 / (3.0 + 3.0*dsqr(a*M/k));

   return mu;
} 

double phi(double z)
{
   double z0, a, b, A;

   z0 = 0.53;
   a  = 1.1;
   b  = 1.1;
   A  = 2.3;

   return pow(z/z0, a) * exp(-pow(z/z0, b)) * A;
}

double M_0_from_f_R0(cosmo *cosmo, double f_R0, double n)
{
   double M_0;

   M_0 = 1.0/R_HUBBLE * sqrt((cosmo->Omega_m + 4.0*cosmo->Omega_de)/((n+1.0)*f_R0));

   return M_0;
}

double mu_0(cosmo *cosmo, error **err)
{
   double kmax, a, M_0, f_R0, sigma, z, dz, mu0, n;

   kmax   = 0.2; // h/Mpc

   // Zhao
   M_0    = 0.02; //h/Mpc
   sigma  = 3;

   // Harnois-D
   n    = 1;
   f_R0 = 3.0e-7;
   M_0    = M_0_from_f_R0(cosmo, f_R0, n);
   printf("M_0(f_R0=%g) = %g\n", f_R0, M_0);

   dz = 0.05;


   for (mu0=0,z=0; z<2.5; z+=dz) {
      a    = 1.0/(z + 1.0);
      mu0 += phi(z) * mu_a_k(cosmo, a, M_0, sigma, n, kmax) * cosmo->Omega_de / Omega_de_a(cosmo, a, -1, err) * dz;
   }

   return mu0;
}

void nofz_bins_equal_number(redshift_t *nofz, error **err)
{
   double zmin, zmax, val, z, norm, dz;
   int nbin, i;
   redshiftANDint rANDi;

   zmin = 0;
   zmax = 2.5;
   nbin = 10;

   rANDi.self = nofz;
   rANDi.i    = 2;    // redshift bin #2

   norm = sm2_qromberg(prob_unnorm, (void*)&rANDi, zmin, zmax, 1.0e-6, err);
   printf("Redshift bins for equal number densities:\n");
   i = 0;
   dz = 0.00001;
   for (z=zmin; z<=zmax+0.2; z+=dz) {
      val = sm2_qromberg(prob_unnorm, (void*)&rANDi, zmin, z, 1.0e-6, err);
      if (val > norm/(double)nbin) {
         printf("%.2f ", zmin);
         zmin = z;
         i++;
      }
   }
   val = sm2_qromberg(prob_unnorm, (void*)&rANDi, zmin, zmax, 1.0e-6, err);
   printf("%.2f ", zmin);
   printf("%.2f\n", zmax);
}

void write_header_p_kappa(FILE *F, int header, cosmo_lens *model, int dimensionless, int Nzbin, error **err)
{
   int i_bin, j_bin;

   if (header >= 2) {
      dump_param_lens(model, F, 1, err);
      forwardError(*err, __LINE__,);

      fprintf(F, "# Shear cross-power spectra ");
      if (dimensionless == 0) {
         fprintf(F, "P(l)");
      } else {
         fprintf(F, "l(l+1)/(2pi)P(l)");
      }
      fprintf(F, " for %d z-bin%s\n", Nzbin, Nzbin==1?"":"s");
   }
   if (header >= 1) {
      fprintf(F, "# l            ");
      for (i_bin=0; i_bin<Nzbin; i_bin++) {
         for (j_bin=i_bin; j_bin<Nzbin; j_bin++) {
            if (model->tomo==tomo_auto_only && i_bin!=j_bin) continue;
            if (model->tomo==tomo_cross_only && i_bin==j_bin) continue;
            fprintf(F, "P_%s^%d%d(l)    ", model->reduced==reduced_K10?"t":"k", i_bin, j_bin);
         }
      }
      fprintf(F, "   ");
      for (i_bin=0; i_bin<Nzbin; i_bin++) {
         for (j_bin=i_bin; j_bin<Nzbin; j_bin++) {
            if (model->tomo==tomo_auto_only && i_bin!=j_bin) continue;
            if (model->tomo==tomo_cross_only && i_bin==j_bin) continue;
            if (model->reduced!=reduced_none) {
               fprintf(F, "P_g1^%d%d(l)   ", i_bin, j_bin);
            }
         }
      }
      fprintf(F, "\n");
   }
}

void usage(int ex)
{
   fprintf(stderr, "Usage: lensingdemo [OPTIONS]\n");
   fprintf(stderr, "OPTIONS\n");
   fprintf(stderr, "  -t FNAME               File containing angular scales in arcmin\n");
   fprintf(stderr, "  -T 'TH_MIN TH_MAX NTH' NTH angular bins between TH_MIN and TH_MAX arcmin\n");
   fprintf(stderr, "                          (default: '%g %g %d')\n",TH_MIN, TH_MAX, NTHETA);
   fprintf(stderr, "  -L 'L_MIN L_MAX NL'    NL Fourier modes between L_MIN and L_MAX (default: '%g %g %d')\n", ELL_MIN, ELL_MAX, NELL);
   fprintf(stderr, "  --linlog MODE          Linear (MODE='LIN') or logarithmic ('LOG'; default) bins for xi+, xi-\n");
   fprintf(stderr, "  -D DIM                 Dimensionless (DIM=1; default) or dimensional (DIM=0) 3D power spectrum\n");
   fprintf(stderr, "  -H LEVEL               Output file header level: 0 (no header), 1 (one line of column description), 2 (full; default)\n");
   fprintf(stderr, "  -q                     Quiet mode\n");
   fprintf(stderr, "  -h                     This message\n");

   if (ex >= 0) exit(ex);
}


int main(int argc, char** argv)
{
   cosmo_lens* model;
   double k, logk, ell, theta, a, pk, pg, f, z, cumG, Ga, ww, incr;
   int i_bin, j_bin, Nzbin, c, dimensionless, verbose, lin, header, iell;
   error *myerr = NULL, **err;
   char *theta_fname, *theta_info, *ell_info, *substr, *slinlog;
   FILE *F;
   size_t Nrec, Ntheta, Nell;
   double *theta_list, Theta_min, Theta_max, ell_min, ell_max, ell_fac;
   double z_mean, z_median;


   err = &myerr;


   theta_fname = theta_info = ell_info = NULL;
   slinlog = "LOG";
   dimensionless = 1;
   verbose = 1;
   lin = 0;
   header = 2;
   while (1) {

      static struct option long_option[] = {
         {"", required_argument, 0, 't'},
         {"", required_argument, 0, 'T'},
         {"", required_argument, 0, 'L'},
         {"", required_argument, 0, 'D'},
         {"", required_argument, 0, 'H'},
         {"linlog", required_argument, 0, 'l'},
         {0, 0, 0, 0}
      };

      int option_index = 0;

      c = getopt_long(argc, argv, "t:T:L:l:qH:D:h", long_option, &option_index);
      switch (c) {
         case 't' :
            theta_fname = optarg;
            break;
         case 'T' :
            theta_info = optarg;
            break;
         case 'L':
            ell_info = optarg;
            break;
         case 'l':
            slinlog = optarg;
            break;
         case 'D':
            dimensionless = atoi(optarg);
            break;
         case 'H':
            header = atoi(optarg);
            break;
         case 'q' :
            verbose = 0;
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

   if (verbose) {
      printf("# lensingdemo (M. Kilbinger, 2006-2012)\n");
   }

   if (dimensionless > 0) {
      dimensionless = 1;
   }

   if (strcmp(slinlog, "lin")==0 || strcmp(slinlog, "LIN")==0) {
      lin = 1;
   } else if (strcmp(slinlog, "log")==0 || strcmp(slinlog, "LOG")==0) {
      lin = 0;
   } else {
      fprintf(stderr, "Wrong linlog mode '%s'\n", slinlog);
      exit(1);
   }



   /* Setting up cosmology */

   F = fopen_err("cosmo_lens.par", "r", err);
   quitOnError(*err, __LINE__, stderr);
   read_cosmological_parameters_lens(&model, F, err);
   quitOnError(*err, __LINE__, stderr);
   fclose(F);

   Nzbin = model->redshift->Nzbin;


   if (verbose) {
      printf("# Cosmological parameters:\n");
      dump_param_lens(model, stdout, 1, err);
      quitOnError(*err, __LINE__, stderr);

      for (i_bin=0; i_bin<Nzbin; i_bin++) {
         z_mean = zmean(model->redshift, i_bin, err);
         quitOnError(*err, __LINE__, stderr);
         z_median = zmedian(model->redshift, i_bin, err);
         quitOnError(*err, __LINE__, stderr);
         printf("Mean/median redshift (bin %d) = %.5f/%.5f\n", i_bin, z_mean, z_median);
      }
   }

   /* w(z) */
   F = fopen_err("w_eos", "w", err);
   quitOnError(*err, __LINE__, stderr);
   dump_param(model->cosmo, F);
   for (a=0.01; a<1.0; a+=0.01) {
      fprintf(F, "%f %f\n", a, w_de(model->cosmo, a, err));
      quitOnError(*err, __LINE__, stderr);
   }
   fclose(F);

   /* Growth factor */
   F = fopen("D_plus", "w");
   for (a=0.1; a<1.0; a+=0.01) {
      fprintf(F, "%f %.6f\n", a, D_plus(model->cosmo, a, 1, err));
      quitOnError(*err, __LINE__, stderr);
   }
   fclose(F);

   /* IST WL code comparison */
   char name[1024];
   int cc_case = 2;
   char author_id[] = "MKNic";

   // WL_case_C0
   /*
   F = fopen("growthfac-MK.dat", "w");
   for (z=0.0; z<2.5; z+=0.01) {
      a = 1.0 / (1.0 + z);
      fprintf(F, "%f %.6f\n", z, D_plus(model->cosmo, a, 1, err));
      quitOnError(*err, __LINE__, stderr);
   }
   fclose(F);
   */

   // lensing efficiency
   sprintf(name, "lensw-C%d-%s.dat", cc_case, author_id);
   F = fopen(name, "w");
   for (z=0.0; z<1.0; z+=0.01) {
      fprintf(F, "%g", z);
      a  = 1.0 / (1.0 + z);
      for (i_bin=0; i_bin<Nzbin; i_bin++) {
         ww = G_bis(model, a, i_bin, err);
         quitOnError(*err, __LINE__, stderr);
         fprintf(F, " %g", ww);
      }
      fprintf(F, "\n");
   }
   fclose(F);

   /* Lensing efficiency */
   F = fopen("G", "w");
   for (a=1.0,cumG=0.0; a>0.1; a-=0.01) {
      ww = w(model->cosmo, a, 0, err);
      quitOnError(*err, __LINE__, stderr);
      Ga = G(model, a, 0, err);
      quitOnError(*err, __LINE__, stderr);
      cumG += Ga * ww;
      fprintf(F, "%f %g %g %g\n", a, Ga * ww, cumG, ww);
   }
   fclose(F);

   /* a(w) */
   F = fopen("a_of_w", "w");
   fprintf(F, "# a w(a)\n");
   double wmax = w(model->cosmo, model->cosmo->a_min, 0, err);
   quitOnError(*err, __LINE__, stderr);
   for (ww=0; ww<wmax; ww+=10) {
      a = a_of_w(model->cosmo, ww, err);
      quitOnError(*err, __LINE__, stderr);
      fprintf(F, "%g %g\n", ww, a);
   }

   /* w_limber2 */
   double w0, w0i, w1, w2, w3;
   F = fopen("w_limber2", "w");
   fprintf(F, "# chi wl2(chi) wl2(chi, interp) wl2' wl2'' wl2'''\n");
   for (ww=0.01; ww<wmax/1.1; ww*=1.02) {
      w0 = w_limber2(model, ww, -1, 0, 0, err);         quitOnError(*err, __LINE__, stderr);
      w0i = w_limber2(model, ww, -1, 0, 1, err);        quitOnError(*err, __LINE__, stderr);
      w1 = deriv_w_23(model, ww, 0, &w2, &w3, err);     quitOnError(*err, __LINE__, stderr);
      fprintf(F, "%g  %g %g  %g %g %g\n", ww, w0, w0i, w1, w2, w3);
   }


   /* Density power spectrum */

   if (verbose) {
      printf("Calculating density power spectrum, writing file P_delta\n");
   }

   if (model->cosmo->nonlinear == coyote10) {
      set_H0_Coyote(model->cosmo, err);
      quitOnError(*err, __LINE__, stderr);
      if (verbose) {
         printf("Coyote10: Hubble parameter set to h = %g\n", model->cosmo->h_100);
      }
   }

   F = fopen("P_delta", "w");
   if (header >= 2) {
      dump_param(model->cosmo, F);
      quitOnError(*err, __LINE__, stderr);

      if (dimensionless == 0) { 
         fprintf(F, "# k[h/Mpc] P_NL(k,a)[(Mpc/h)^3]\n");
      } else {
         fprintf(F, "# k[h/Mpc] Delta_NL(k,a)\n");
      }
   }
   if (header >= 1) {
      fprintf(F, "# k           ");
      for (a=1.0; a>=0.5; a-=0.5) {
         fprintf(F, "        a=%g", a);
      }
      fprintf(F, "\n");
   }

   for (k=0.0001; k<=100.0; k*=1.05) {
      if (dimensionless == 0) {
         f = 1.0;
      } else {
         f = k*k*k/(2.0*pi_sqr);
      }
      fprintf(F, "%e ", k);
      for (a=1.0; a>=0.5; a-=0.5) {
         if (a > model->cosmo->a_min) {
            fprintf(F, "%e ", P_NL(model->cosmo, a, k, err)*f);
            quitOnError(*err, __LINE__, stderr);
         } else {
	    fprintf(F, "undef ");
	 }
      }
      fprintf(F, "\n");
      quitOnError(*err, __LINE__, stderr);
   }
   fclose(F);

   // WL_case_C0
   F = fopen("pkzero-lin-MK.dat", "w");
   quitOnError(*err, __LINE__, stderr);

   for (logk=-3; logk<=1.0; logk+=0.01) {
      k  = pow(10, logk) / model->cosmo->h_100; // k [Mpc^[-1] = k/h [h/Mpc] 
      pk = P_L(model->cosmo, 1.0, k, err) / pow(model->cosmo->h_100, 3);  // Pk [Mpc^3/h^3] = Pk/h^3
      quitOnError(*err, __LINE__, stderr);
      fprintf(F, "%8g %8g\n", logk, log10(pk));
   }
   fclose(F);

   // XC_case_C1
   F = fopen("P_delta_nicaea_XC_1.dat", "w");
   quitOnError(*err, __LINE__, stderr);

   for (logk=-5; logk<=3.0; logk+=0.01) {
      k  = pow(10, logk); // / model->cosmo->h_100; // k [Mpc^[-1] = k/h [h/Mpc] 
      pk = P_L(model->cosmo, 0.5, k, err); // / pow(model->cosmo->h_100, 3);  // Pk [Mpc^3/h^3] = Pk/h^3
      quitOnError(*err, __LINE__, stderr);
      fprintf(F, "%8g %8g\n", k, pk);
   }
   fclose(F);



   /* Redshift distribution */
   if (verbose) {
      printf("Printing redshift distribution to file nz\n");
   }
   F = fopen("nz", "w");
   dump_param_lens(model, F, 1, err);
   quitOnError(*err, __LINE__, stderr);

   for (z=0; z<=10; z+=0.01) {
      fprintf(F, "%.3f", z);
      for (i_bin=0; i_bin<Nzbin; i_bin++) {
         fprintf(F, " %g", prob(model->redshift, z, i_bin, err));
         quitOnError(*err, __LINE__, stderr);
      }
      fprintf(F, "\n");
   }
   fclose(F);

   //nofz_bins_equal_number(model->redshift, err); quitOnError(*err, __LINE__, stderr);

   /* Convergence power spectrum */

   if (verbose) {
      printf("Calculating convergence power spectrum for discrete ell, writing file P_kappa_d\n");
   }

   F = fopen("P_kappa_d", "w");
   write_header_p_kappa(F, header, model, dimensionless, Nzbin, err);
   quitOnError(*err, __LINE__, stderr);
   iell = 2;
   //while (iell <= 1e5) {
   while (iell < 400) {
      if (dimensionless == 0) {
         f = 1.0;
      } else {
         f = iell*(iell+1.0)/twopi;
      }
      fprintf(F, "%d  ", iell);
      for (i_bin=0; i_bin<Nzbin; i_bin++) {
         for (j_bin=i_bin; j_bin<Nzbin; j_bin++) {
            if (model->tomo==tomo_auto_only && i_bin!=j_bin) continue;
            if (model->tomo==tomo_cross_only && i_bin==j_bin) continue;
            pk = Pshear_spherical(model, iell, i_bin, j_bin, err);
            quitOnError(*err, __LINE__, stderr);
            fprintf(F, " %e", pk*f);
         }
      }

      fprintf(F, "\n");
      fflush(F);

      //iell++;
      if (iell<=10) iell++;
      else if (iell<=20) iell += 2;
      else if (iell<=50) iell += 5;
      else iell += 10;
   }

   fclose(F);

   if (model->projection != full) {

      if (verbose) {
         printf("Calculating convergence power spectrum for continuous ell, writing file P_kappa\n");
      }

      F = fopen("P_kappa", "w");
      write_header_p_kappa(F, header, model, dimensionless, Nzbin, err);
      quitOnError(*err, __LINE__, stderr);

      if (ell_info != NULL) {
         substr    = strsep(&ell_info, " "); ell_min = atof(substr);
         substr    = strsep(&ell_info, " "); ell_max = atof(substr);
         substr    = strsep(&ell_info, " "); Nell    = atoi(substr);
      } else {
         ell_min   = ELL_MIN;
         ell_max   = ELL_MAX;
         Nell      = NELL;   
      }
      ell_fac = pow( ell_max / ell_min, 1.0 / ((double)Nell - 1));


      for (ell=ell_min; ell<=ell_max*1.001; ell*=ell_fac) {
         if (dimensionless == 0) {
            f = 1.0;
         } else {
            f = ell*(ell+1.0)/twopi;
         }
         fprintf(F, "%e  ", ell);
         for (i_bin=0; i_bin<Nzbin; i_bin++) {
            for (j_bin=i_bin; j_bin<Nzbin; j_bin++) {
               if (model->tomo==tomo_auto_only && i_bin!=j_bin) continue;
               if (model->tomo==tomo_cross_only && i_bin==j_bin) continue;
               pk = Pshear(model, ell, i_bin, j_bin, err);
               if (getErrorValue(*err)==reduced_fourier_limit) {
                  purgeError(err);
               }
               quitOnError(*err, __LINE__, stderr);
               fprintf(F, " %e", pk*f);
            }
         }
         fprintf(F, "   ");

         if (model->reduced!=reduced_none) {
            for (i_bin=0; i_bin<Nzbin; i_bin++) {
               for (j_bin=i_bin; j_bin<Nzbin; j_bin++) {
                  if (model->tomo==tomo_auto_only && i_bin!=j_bin) continue;
                  if (model->tomo==tomo_cross_only && i_bin==j_bin) continue;
                  pg = Pg1(model, ell, i_bin, j_bin, err);
                  if (getErrorValue(*err)==reduced_fourier_limit) {
                     purgeError(err);
                  }
                  quitOnError(*err, __LINE__, stderr);
                  fprintf(F, " %e", pg*f);
               }
            }
         }

         fprintf(F, "\n");
      }

      fclose(F);

   } else {

      printf("For projection = full, real-space correlation functions are not yet implemented, exiting now...\n");
      return 0;
   }



   /* Shear correlation functions */

   if (theta_fname != NULL) {

      Nrec = 0;
      theta_list = (double*)read_any_list_count(theta_fname, &Nrec, "%lg", sizeof(double), &Ntheta, err);
      quitOnError(*err, __LINE__, stderr);

   } else {

      /* Create angular bins */
      if (theta_info != NULL) {
         substr    = strsep(&theta_info, " "); Theta_min = atof(substr);
         substr    = strsep(&theta_info, " "); Theta_max = atof(substr);
         substr    = strsep(&theta_info, " "); Ntheta    = atoi(substr);
      } else {
         Theta_min  = TH_MIN;
         Theta_max  = TH_MAX;
         Ntheta     = NTHETA;
      }

      theta_list = malloc_err(sizeof(double) * Ntheta, err);
      quitOnError(*err, __LINE__, stderr);
      for (c=0; c<Ntheta; c++) {
         incr = (double)c / ((double)Ntheta - 1.0);
         if (lin == 0) {
            theta_list[c] = Theta_min * pow(Theta_max / Theta_min, incr);
         } else {
            theta_list[c] = Theta_min + incr * (Theta_max - Theta_min);
         }
      }

   }


   if (verbose) {
      printf("Calculating shear correlation function, writing files xi_p, xi_m\n");
   }

   F = fopen("xi_p", "w");
   if (header >= 2) {
      dump_param_lens(model, F, 1, err);
      quitOnError(*err, __LINE__, stderr);
   }
   if (header >= 1) {
      fprintf(F, "# theta[arcmin]   xi+^{mn}(theta)\n# th[']  ");
      for (i_bin=0; i_bin<Nzbin; i_bin++) {
         for (j_bin=i_bin; j_bin<Nzbin; j_bin++) {
            if (model->tomo==tomo_auto_only && i_bin!=j_bin) continue;
            if (model->tomo==tomo_cross_only && i_bin==j_bin) continue;
            fprintf(F, "           %d%d", i_bin, j_bin);
         }
      }
      fprintf(F, "\n");
   }

   for (c=0; c<Ntheta; c++) {
      theta = theta_list[c] * arcmin;
      fprintf(F, "%20.8f", theta/arcmin);
      for (i_bin=0; i_bin<Nzbin; i_bin++) {
         for (j_bin=i_bin; j_bin<Nzbin; j_bin++) {
            if (model->tomo==tomo_auto_only && i_bin!=j_bin) continue;
            if (model->tomo==tomo_cross_only && i_bin==j_bin) continue;
            fprintf(F, " %e", xi(model, +1, theta, i_bin, j_bin, err));
            quitOnError(*err, __LINE__, stderr);
         }
      }
      fprintf(F, "\n");
   }
   fclose(F);

   F = fopen("xi_m", "w");
   if (header >= 2) {
      dump_param_lens(model, F, 1, err);
      quitOnError(*err, __LINE__, stderr);
   }
   if (header >= 1) {
      fprintf(F, "# theta[arcmin] xi-^{mn}(theta)\n# th[']  ");
      for (i_bin=0; i_bin<Nzbin; i_bin++) {
         for (j_bin=i_bin; j_bin<Nzbin; j_bin++) {
            if (model->tomo==tomo_auto_only && i_bin!=j_bin) continue;
            if (model->tomo==tomo_cross_only && i_bin==j_bin) continue;
            fprintf(F, "           %d%d", i_bin, j_bin);
         }
      }
      fprintf(F, "\n");
   }

   for (c=0; c<Ntheta; c++) {
      theta = theta_list[c] * arcmin;
      fprintf(F, "%20.8f", theta/arcmin);
      for (i_bin=0; i_bin<Nzbin; i_bin++) {
         for (j_bin=i_bin; j_bin<Nzbin; j_bin++) {
            if (model->tomo==tomo_auto_only && i_bin!=j_bin) continue;
            if (model->tomo==tomo_cross_only && i_bin==j_bin) continue;
            fprintf(F, " %e", xi(model, -1, theta, i_bin, j_bin, err));
            quitOnError(*err, __LINE__, stderr);
         }
      }
      fprintf(F, "\n");
   }
   fclose(F);


   /* Tophat shear variance */

   if (verbose) {
      printf("Calculating tophat shear variance, writing file gammasqr\n");
   }

   F = fopen("gammasqr", "w");
   dump_param_lens(model, F, 1, err);
   quitOnError(*err, __LINE__, stderr);

   fprintf(F, "# theta[arcmin] <|gamma|^2>^{mn}(theta)\n#        ");
   for (i_bin=0; i_bin<Nzbin; i_bin++) {
      for (j_bin=i_bin; j_bin<Nzbin; j_bin++) {
         if (model->tomo==tomo_auto_only && i_bin!=j_bin) continue;
         if (model->tomo==tomo_cross_only && i_bin==j_bin) continue;
         fprintf(F, "           %d%d", i_bin, j_bin);
      }
   }
   fprintf(F, "\n");

   for (c=0; c<Ntheta; c++) {
      theta = theta_list[c] * arcmin;
      fprintf(F, "%12.6f", theta/arcmin);
      for (i_bin=0; i_bin<Nzbin; i_bin++) {
         for (j_bin=i_bin; j_bin<Nzbin; j_bin++) {
            if (model->tomo==tomo_auto_only && i_bin!=j_bin) continue;
            if (model->tomo==tomo_cross_only && i_bin==j_bin) continue;
            fprintf(F, " %e", gamma2(model, theta, i_bin, j_bin, err));
            quitOnError(*err, __LINE__, stderr);
         }
      }
      fprintf(F, "\n");
   }
   fclose(F);


   /* Aperture mass variance */

   if (verbose) {
      printf("Calculating aperture mass variance (polynomial filter), writing files mapsqr, mapsqr_gauss\n");
   }

   F = fopen("mapsqr", "w");
   if (header >= 2) {
      dump_param_lens(model, F, 1, err);
      quitOnError(*err, __LINE__, stderr);

      fprintf(F, "# M = <M_ap^2>(theta){polynomial}\n");
   }

   if (header >= 1) {
      fprintf(F, "# theta[arcmin]");
      for (i_bin=0; i_bin<Nzbin; i_bin++) {
         for (j_bin=i_bin; j_bin<Nzbin; j_bin++) {
            if (model->tomo==tomo_auto_only && i_bin!=j_bin) continue;
            if (model->tomo==tomo_cross_only && i_bin==j_bin) continue;
            fprintf(F, "      M^{%d,%d}", i_bin, j_bin);
         }
      }
      fprintf(F, "\n");
   }

   for (c=0; c<Ntheta; c++) {
      theta = theta_list[c] * arcmin;
      fprintf(F, "%15.6f", theta/arcmin);
      for (i_bin=0; i_bin<Nzbin; i_bin++) {
         for (j_bin=i_bin; j_bin<Nzbin; j_bin++) {
            if (model->tomo==tomo_auto_only && i_bin!=j_bin) continue;
            if (model->tomo==tomo_cross_only && i_bin==j_bin) continue;
            fprintf(F, " %e", map2_poly(model, theta, i_bin, j_bin, err));
            quitOnError(*err, __LINE__, stderr);
         }
      }
      fprintf(F, "\n");
   }
   fclose(F);

   F = fopen("mapsqr_gauss", "w");
   dump_param_lens(model, F, 1, err);
   quitOnError(*err, __LINE__, stderr);

   fprintf(F, "# theta[arcmin] <M_ap^2>^{mn}(theta){Gaussian}\n#        ");
   for (i_bin=0; i_bin<Nzbin; i_bin++) {
      for (j_bin=i_bin; j_bin<Nzbin; j_bin++) {
	 if (model->tomo==tomo_auto_only && i_bin!=j_bin) continue;
	 if (model->tomo==tomo_cross_only && i_bin==j_bin) continue;
	 fprintf(F, "           %d%d", i_bin, j_bin);
      }
   }
   fprintf(F, "\n");

   for (c=0; c<Ntheta; c++) {
      theta = theta_list[c] * arcmin;
      fprintf(F, "%12.6f", theta/arcmin);
      for (i_bin=0; i_bin<Nzbin; i_bin++) {
         for (j_bin=i_bin; j_bin<Nzbin; j_bin++) {
            if (model->tomo==tomo_auto_only && i_bin!=j_bin) continue;
            if (model->tomo==tomo_cross_only && i_bin==j_bin) continue;
            fprintf(F, " %e", map2_gauss(model, theta, i_bin, j_bin, err));
            quitOnError(*err, __LINE__, stderr);
         }
      }
      fprintf(F, "\n");
   }
   fclose(F);


   if (verbose) {
      printf("Done\n");
   }

   free_parameters_lens(&model);

   return 0;
}
#undef TH_MIN
#undef TH_MAX
#undef NTHETA
#undef ELL_MIN
#undef ELL_MAX
#undef NELL
