/* ============================================================ *
 * halomodeldemo.c						*
 * Martin Kilbinger 2008-2010					*
 * ============================================================ */

#include <stdio.h>
#include <stdlib.h>

#include "errorlist.h"
#include "io.h"
#include "maths.h"
#include "cosmo.h"
#include "nofz.h"
#include "halomodel.h"
#include "hod.h"


int main()
{
   cosmo_hm *hm;
   error *myerr = NULL, **err;
   FILE *F;
   double a;
   int i_bin;
   double k, fac, P1h, P2h, Pl, Pnl, P1hg, P2hg, zm;


   err = &myerr;

   printf("# halomodeldemo (M. Kilbinger, 2009)\n");

   F = fopen_err("halomodel.par", "r", err);
   quitOnError(*err, __LINE__, stderr);
   read_cosmological_parameters_hm(&hm, F, err);
   quitOnError(*err, __LINE__, stderr);
   fclose(F);

   dump_param_hm(hm, stdout, err);
   quitOnError(*err, __LINE__, stderr);

   double R, sigma_R, dsigma2_R_dR;
   for (R=0.1; R<700; R*=1.2) {
      sigma_R = sqrt(sigma_R_sqr(hm, R, err));
      quitOnError(*err, __LINE__, stderr);
      dsigma2_R_dR = dsigma_R_sqr_dR(hm, R, err);
      quitOnError(*err, __LINE__, stderr);
      printf("%g %g %g\n", R, sigma_R, dsigma2_R_dR);
   }


   for (i_bin=0; i_bin<hm->redshift->Nzbin; i_bin++) {

      zm = zmean(hm->redshift, i_bin, err);
      quitOnError(*err, __LINE__, stderr);

      fprintf(stderr, "Mean redshift (bin %d) = %.3f\n", i_bin, zm);
      fprintf(stderr, "Comoving volume between z=0 and %f (bin %d) = %g (h^-3 Mpc^3 per solid angle\n",
              zm, i_bin, vc(hm, 0.0, zm, err));
      quitOnError(*err, __LINE__, stderr);
   }

   fprintf(stderr, "Calculating power spectra...\n");
   a = 0.99;
   F = fopen_err("Phalo", "w", err);
   quitOnError(*err, __LINE__, stderr);
   fprintf(F, "# k P1h_dm P2h_dm P1h_g P2h_g P_L P_NL\n");
   for (k=0.01; k<300; k*=1.25) {
      P1hg = P2hg = P1h = P2h = Pl = Pnl = 0.0;
      fac  = k*k*k/(2*pi*pi);
      P1h  = P1h_dm(hm, a, k, err);                      quitOnError(*err, __LINE__, stderr);
      P2h  = P2h_dm(hm, a, k, err);                      quitOnError(*err, __LINE__, stderr);
      //P1hg = P1h_g(hm, a, k, err);                       quitOnError(*err, __LINE__, stderr);
      //P2hg = FFTLog_P2h_g_hexcl(hm, a, k, 20, err);          quitOnError(*err, __LINE__, stderr);

      Pl   = P_L(hm->cosmo, a, k, err);                  quitOnError(*err, __LINE__, stderr); 
      Pnl  = P_NL(hm->cosmo, a, k, err);                 quitOnError(*err, __LINE__, stderr); 
      fprintf(F, "%f %e %e %e %e %e %e\n", k,
	      P1h*fac, P2h*fac, P1hg*fac, P2hg*fac, Pl*fac, Pnl*fac);
      //fprintf(stdout, "%f %e %e %e %e %e %e\n", k,
      //      P1h*fac, P2h*fac, P1hg*fac, P2hg*fac, Pl*fac, Pnl*fac); fflush(stdout);
   }
   fclose(F);


   /* Clean up */

   free_parameters_hm(&hm);

   return 0;
}
