
#include "sn1a.h"
#include "cosmo.h"
#include "io.h"
#include "maths.h"
#include "errorlist.h"


#define CHI2_BIG 10000.0


void print_dl(error **err)
{
   cosmo_SN *self;
   double z, dl;
   FILE *F;

   F = fopen_err("cosmo_SN.par", "r", err);
   forwardError(*err, __LINE__,);
   read_cosmological_parameters_SN(&self, F, err);
   forwardError(*err, __LINE__,);
   fclose(F);

   F = fopen_err("D_Lum", "w", err);
   forwardError(*err, __LINE__,);
   dump_param_SN(self, F);
   fprintf(F, "# z D_{lum}\n");

   for (z=0.01; z<1.1; z+=0.01) {
      dl = D_lum(self->cosmo, 1/(1.0+z), err);
      forwardError(*err, __LINE__,);
      //mu = distance_module(self->cosmo, dl, err);
      //forwardError(*err, __LINE__,);
      fprintf(F, "%g %g\n", z, dl);
   }
   fclose(F);
}

int main(int argc, char *argv[])
{
   error *myerr = NULL, **err;

   err = &myerr;

   print_dl(err);
   quitOnError(*err, __LINE__, stderr);

   return 0;
}
