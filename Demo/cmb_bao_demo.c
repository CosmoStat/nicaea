/* ============================================================ *
 * cmb_bao_demo.c						*
 * Martin Kilbinger 2010					*
 * ============================================================ */

#include <stdio.h>
#include <stdlib.h>

#include "errorlist.h"
#include "io.h"
#include "maths.h"
#include "cosmo.h"
#include "cmb_bao.h"


int main(int argc, char** argv)
{
   cosmo *model;
   FILE *F;
   error *myerr = NULL, **err;
   double z, z2, f, f2;

   err = &myerr;


   F = fopen_err("cosmoDP.par", "r", err);
   quitOnError(*err, __LINE__, stderr);
   read_cosmological_parameters(&model, F, err);
   quitOnError(*err, __LINE__, stderr);
   fclose(F);

   printf("# Cosmological parameters:\n");
   dump_param(model, stdout);
   printf("a_min = %g\n", model->a_min);

   z = z_drag(model);
   fprintf(stderr, "Drag epoch, baryon decoupling redshift z_d = %g\n", z);

   f = r_sound_integral(model, 1.0/(1.0+z), err);
   quitOnError(*err, __LINE__, stderr);
   fprintf(stderr, "Sound horizon num integ  (at z_d=%g) = %g Mpc/h (%g Mpc)\n", z, f, f/model->h_100);

   f = r_sound_drag_analytical(model, err);
   quitOnError(*err, __LINE__, stderr);
   fprintf(stderr, "Sound horizon analytical (at z_d=%g) = %g Mpc/h (%g Mpc)\n", z, f, f/model->h_100);

   f = r_sound_drag_fit(model, err);
   quitOnError(*err, __LINE__, stderr);
   fprintf(stderr, "Sound horizon EH fit     (at z_d=%g) = %g Mpc/h (%g Mpc)\n", z, f, f/model->h_100);

   z = 0.35;
   f = D_V(model, 1.0/(1.0+z), err);
   quitOnError(*err, __LINE__, stderr);
   printf("Spherically-averaged distance D_V(z=%g) = %g Mpc/h (%g Mpc)\n", z, f, f/model->h_100);

   z2 = 0.2;
   f2 = D_V(model, 1.0/(1.0+z2), err);
   quitOnError(*err, __LINE__, stderr);
   printf("Spherically-averaged distance D_V(z=%g) = %g Mpc/h  (%g Mpc)\n", z2, f2, f2/model->h_100);
   printf("Ratio D_V(%g)/D_V(%g) = %g\n", z, z2, f/f2);

   z = 0.275;
   f = D_V(model, 1.0/(1.0+z), err);
   quitOnError(*err, __LINE__, stderr);
   printf("Spherically-averaged distance D_V(z=%g) = %g Mpc/h (%g Mpc)\n", z, f, f/model->h_100);


   return 0;   
}
