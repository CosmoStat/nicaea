/* ============================================================ *
 * nofz.c						        *
 * Martin Kilbinger 2009					*
 * ============================================================ */

#include "nofz.h"

/* ============================================================ *
 * Initialised redshift without filling in data.		*
 * ============================================================ */
redshift_t *init_redshift_empty(int Nzbin, int Nnz_max, error **err)
{
   redshift_t *res;

   res = malloc_err(sizeof(redshift_t), err);   forwardError(*err, __LINE__, NULL);
 
   res->Nzbin     = Nzbin;
   res->Nnz       = calloc_err(Nzbin, sizeof(int), err);                forwardError(*err, __LINE__, NULL);
   res->Nnz_max   = Nnz_max;
   res->nofz      = calloc_err(Nzbin, sizeof(nofz_t), err);             forwardError(*err, __LINE__, NULL);
   res->photz     = calloc_err(Nzbin, sizeof(photz_t), err);            forwardError(*err, __LINE__, NULL);
   res->par_nz    = calloc_err(Nzbin*Nnz_max, sizeof(double), err);     forwardError(*err, __LINE__, NULL);
   res->prob_norm = calloc_err(Nzbin, sizeof(double), err);             forwardError(*err, __LINE__, NULL);
   res->z_rescale = calloc_err(Nzbin, sizeof(double), err);             forwardError(*err, __LINE__, NULL);

   return res;
}

/* ============================================================= *
 * Initialises redshift, copies content of the Nzbin-dimensional *
 * arrays Nnz, nofz, par_nz and z_rescale.			 *
 * If z_rescale==NULL, initialises with 1.			 *
 * ============================================================= */
redshift_t *init_redshift(int Nzbin, const int *Nnz, const nofz_t *nofz, const photz_t *photz,
            const double *par_nz, const double *z_rescale, error **err)
{
   redshift_t *res;
   int n;

   /* The upper limit 1000 is arbitrary, mainly a test whether Nzbin is initialised properly. *
    * However, we don't claim that the code works with 999 redshift bins.		      */
   testErrorRetVA(Nzbin<=0 || Nzbin>1000, redshift_Nzbin, "Number of redshift bins (Nzbin=%d) not valid",
		  *err, __LINE__, NULL, Nzbin);

   for (n=0; n<Nzbin; n++) {
      testErrorRetVA(Nnz[n]<2, redshift_nnz, "Number of redshift parameters (%d, bin #%d) must be >= 2",
		     *err, __LINE__, 0, Nnz[n], n);
   }

   res = malloc_err(sizeof(redshift_t), err);   forwardError(*err, __LINE__, NULL);

   res->Nzbin = Nzbin;

   res->Nnz = malloc_err(Nzbin*sizeof(int), err); forwardError(*err, __LINE__, NULL);
   memcpy(res->Nnz, Nnz, Nzbin*sizeof(int));

   for (n=0,res->Nnz_max=-1; n<Nzbin; n++) {
      if (res->Nnz[n]>res->Nnz_max) res->Nnz_max = res->Nnz[n];
   }
   
   res->nofz = malloc_err(Nzbin*sizeof(nofz_t), err); forwardError(*err, __LINE__, NULL);
   memcpy(res->nofz, nofz, Nzbin*sizeof(nofz_t));

   res->photz = malloc_err(Nzbin*sizeof(photz_t), err); forwardError(*err, __LINE__, NULL);
   memcpy(res->photz, photz, Nzbin*sizeof(photz_t));

   res->par_nz = malloc_err(Nzbin*res->Nnz_max*sizeof(double), err);
   forwardError(*err, __LINE__, NULL);
   memcpy(res->par_nz, par_nz, Nzbin*res->Nnz_max*sizeof(double));

   res->prob_norm = calloc_err(Nzbin, sizeof(double), err);
   forwardError(*err, __LINE__, NULL);

   res->z_rescale = malloc_err(Nzbin*sizeof(double), err);
   forwardError(*err, __LINE__, NULL);
   for (n=0; n<Nzbin; n++) {
      if (z_rescale!=NULL) {
         res->z_rescale[n] = z_rescale[n];
      } else {
         res->z_rescale[n] = 1.0;
      }
   }

   return res;
}

redshift_t *copy_redshift(const redshift_t *REDSHIFT, error **err)
{
   redshift_t *res;

   res = init_redshift(REDSHIFT->Nzbin, REDSHIFT->Nnz, REDSHIFT->nofz, REDSHIFT->photz, REDSHIFT->par_nz,
                       REDSHIFT->z_rescale, err);
   forwardError(*err, __LINE__, NULL);

   testErrorRetVA(REDSHIFT->Nnz_max!=res->Nnz_max, redshift_nnz, "Inconsistent Nz_max (%d!=%d)",
		  *err, __LINE__, NULL, REDSHIFT->Nnz_max, res->Nnz_max);

   return res;
}


void updateFrom_redshift(redshift_t *avant, redshift_t *apres)
{
   if (change_redshift(avant, apres)) {
      free(apres->prob_norm);
      apres->prob_norm = NULL;
   }
}

void free_redshift(redshift_t **self)
{
   redshift_t *s;

   s = *self;
   free(s->Nnz);          s->Nnz = NULL;
   free(s->nofz);         s->nofz = NULL;
   free(s->photz);        s->photz = NULL;
   free(s->par_nz);       s->par_nz = NULL;
   free(s->prob_norm);    s->prob_norm = NULL;
   free(s->z_rescale);    s->z_rescale = NULL;
   free(s);
   s = NULL;
}

/* Useful if redshift changes */
void free_and_reinit_redshift(redshift_t **self, int Nzbin, int Nnz_max, error **err)
{
   if (*self!=NULL) free_redshift(self);
   *self = init_redshift_empty(Nzbin, Nnz_max, err);
   forwardError(*err, __LINE__,);
}

/* ============================================================ *
 * The mother of all functions to read redshift stuff from      *
 * files. Also used by wrappers/nz.c:nz_init.			*
 * Reads from a config-type file:				*
 *    Nzbin	[integer]					*
 *    snzmode	[string]					*
 * Next, depending on nzmode, does				*
 * - nothing (nz_fit_hist)					*
 * - reads							*
 *      nzfile     [Nzbin strings]				*
 *   and reads redshift from files (nz_read_from_files)         *
 * ============================================================ */
void read_redshift_info(redshift_t **self, FILE *F, error **err)
{
   config_element c = {0, 0.0, ""};
   struct {
     int Nzbin;
     char snzmode[128], **nzfile;
     nzmode_t nzmode;
   } tmp;
   int j;
   char s[128];

   CONFIG_READ(&tmp, Nzbin, i, F, c, err);
   CONFIG_READ_S(&tmp, snzmode, s, F, c, err);
   STRING2ENUM(tmp.nzmode, tmp.snzmode, nzmode_t, snzmode_t, j, Nnzmode_t, err);

   switch (tmp.nzmode) {

      case nz_fit_hist :
         *self = NULL;
         return;           /* Used only in wrappers:nz.c */

      case nz_read_from_files :
         tmp.nzfile = malloc_err(sizeof(char*)*tmp.Nzbin, err);
         forwardError(*err, __LINE__,);
         for (j=0; j<tmp.Nzbin; j++) {
            tmp.nzfile[j] = malloc_err(sizeof(char)*128*tmp.Nzbin, err);
         }
         CONFIG_READ_S_ARR(&tmp, nzfile, s, j, tmp.Nzbin, F, c, err);

         *self = init_redshift_from_files((const char**)(tmp.nzfile), tmp.Nzbin, err);
         forwardError(*err, __LINE__,);
         break;

      default :
         *err = addErrorVA(redshift_unknown, "Unknown nzmode type %d",
               *err, __LINE__, tmp.nzmode);
         return;
   }
}

/* Reads parameters of slice n_bin from file name. Used by init_redshift_from_files. */
void read_redshift_slice(redshift_t *self, int n_bin, const char *name, error **err)
{
   nofz_t nofz;
   photz_t photz;
   int i;
   size_t Nrec, Nlines;
   double *ptr;

   get_nofz_t_file(name, &nofz, &photz, err);
   forwardError(*err, __LINE__,);

   self->nofz[n_bin]  = nofz;
   self->photz[n_bin] = photz;

   if (nofz == single) {
      testErrorRet(photz != photz_no, redshift_type,
         "Single-redshift n(z) cannot have photometric error type other than 'sphotz=photz_no'",
         *err, __LINE__,);
   }

   if (nofz==hist) {

      testErrorRet(photz != photz_no, redshift_type,
         "Histogram n(z) cannot have photometric error type other than 'sphotz=photz_no'",
         *err, __LINE__,);
      ptr = read_par_nz_hist(name, self->Nnz+n_bin, err);
      forwardError(*err, __LINE__,);
      memcpy(self->par_nz + n_bin*self->Nnz_max, ptr, self->Nnz[n_bin]*sizeof(double));

    } else {

      Nrec = 0;
      ptr = (double*)read_any_list_count(name, &Nrec, "%lg", sizeof(double), &Nlines, err);
      testErrorRetVA(Nrec!=Nlines, redshift_nnz,
            "Number of records (%d) not equal to number of (uncommented) lines (%d) in n(z) file '%s'",
            *err, __LINE__,, Nrec, Nlines, name);
      self->Nnz[n_bin] = Nlines;
      for (i=0; i<Nrec; i++) {
         self->par_nz[n_bin*self->Nnz_max+i] = ptr[i];
      }

   }

   free(ptr);
}

void dump_redshift(redshift_t *self, FILE *F, error **err)
{
   int n_bin;

   dump_redshift_nostruct(self->Nzbin, self->Nnz, self->nofz, self->photz, self->par_nz, self->Nnz_max, F);
   fprintf(F, "# Mean redshift =");
   for (n_bin=0; n_bin<self->Nzbin; n_bin++) {
      fprintf(F, " %.3f", zmean(self, n_bin, err));
      forwardError(*err, __LINE__,);
   }
   fprintf(F, "\n");
}

void dump_redshift_nostruct(int Nzbin, int *Nnz, nofz_t *nofz, photz_t *photz, double *par_nz, int Nnz_max, FILE *F)
{
   int n, m;

   /* Number of redshift bins */
   fprintf(F, "# Nzbin = %d\n#", Nzbin);

   /* Redshift distribution type */
   for (n=0; n<Nzbin; n++) fprintf(F, " nofz%d", n);
   fprintf(F, "\n#");
   for (n=0; n<Nzbin; n++) fprintf(F, "     %d", nofz[n]);
   fprintf(F, " (");
   for (n=0; n<Nzbin-1; n++) fprintf(F, "%s ", snofz_t(nofz[n]));
   fprintf(F, "%s)", snofz_t(nofz[n]));
   fprintf(F, "\n# ");

   /* Photometric error type */
   for (n=0; n<Nzbin; n++) fprintf(F, " photz%d", n);
   fprintf(F, "\n#");
   for (n=0; n<Nzbin; n++) fprintf(F, "       %d", photz[n]);
   fprintf(F, " (");
   for (n=0; n<Nzbin-1; n++) fprintf(F, "%s ", sphotz_t(photz[n]));
   fprintf(F, "%s)", sphotz_t(photz[n]));
   fprintf(F, "\n# ");

   /* Number of redshift parameters */
   for (n=0; n<Nzbin; n++) fprintf(F, " Nnz%d ", n);
   fprintf(F, "\n#");
   for (n=0; n<Nzbin; n++) fprintf(F, "     %d", Nnz[n]);
   fprintf(F, "\n#  ");

   /* Redshift parameters */
   for (n=0; n<Nzbin; n++) {
      for (m=0; m<Nnz[n]; m++) {
         fprintf(F, "pnz%d   ", m);
      }
      fprintf(F, "   ");
   }
   fprintf(F, "\n# ");
   for (n=0; n<Nzbin; n++) {
      for (m=0; m<Nnz[n]; m++) {
         fprintf(F, "%.3f  ", par_nz[n*Nnz_max+m]);
      }
      fprintf(F, "   ");
   }
   fprintf(F, "\n");
}


/* Used in likeli_Lensing and likeli_Nz, to copy parameter values into redshift structure */
void fill_redshift_slice(redshift_t *self, int n_bin, nofz_t nofz, photz_t photz, error **err, ...)
{
   va_list args;
   int Nnz_max;

   va_start(args, err);

   self->nofz[n_bin]  = nofz;
   self->photz[n_bin] = photz;
   Nnz_max            = self->Nnz_max;

   switch (nofz) {

      case ymmk :
         self->Nnz[n_bin] = 5;
         self->par_nz[n_bin*Nnz_max+0] = va_arg(args, double);   /* z_min */
         self->par_nz[n_bin*Nnz_max+1] = va_arg(args, double);   /* z_max */
         self->par_nz[n_bin*Nnz_max+2] = va_arg(args, double);   /* a_ymmk */
         self->par_nz[n_bin*Nnz_max+3] = va_arg(args, double);   /* b_ymmk */
         self->par_nz[n_bin*Nnz_max+4] = va_arg(args, double);   /* c_ymmk */
         testErrorRet(photz != photz_no, redshift_type,
            "Histogram n(z) cannot have photometric error type other than 'sphotz=photz_no'",
            *err, __LINE__,);
         break;

      default :
         *err = addErrorVA(redshift_unknown, "Unknown redshift type %d\n", *err, __LINE__, nofz);
         return;
   }
}

/* ============================================================ *
 * Reads a redshift histogram from the file "name" with columns	*
 *    z_i(left bin corner)   n_i				*
 * The last row contains					*
 *    z_max x							*
 * where x is not used and can be any number.			*
 * Returns vector with redshift parameters and sets Nnz         *
 * accordingly.							*
 * ============================================================ */
double *read_par_nz_hist(const char *name, int *Nnz, error **err)
{
   FILE *F;
   char c, cpre;
   int i, n, nread;
   double *par_nz, dummy;
   char str[1024];
   char *dummy2;
   size_t size;

   /* Number of lines */
   F = fopen(name, "r");
   if (F==NULL) {
      sprintf(str, "Could not open file %s", name);
      *err = addError(io_file, str, *err, __LINE__);
      return NULL;
   }
   cpre = '\n';
   n = 0;
   do {
      c = fgetc(F);
      if (c=='#' && cpre=='\n') n--;
      if (c=='\n') n++;
      cpre = c;
   } while (c!=EOF);
   testErrorRet(n==0, io_file, "File corrupt: no valid line found in n(z) file", *err, __LINE__, NULL);
   n--;   /* Last line contains z_n 0 */

   *Nnz = 2*n+1;
   par_nz = malloc_err(sizeof(double)*(*Nnz), err);  forwardError(*err, __LINE__, NULL);

   /* Read parameters */
   /* TODO: Replace following by read_any_list */
   rewind(F);

   /* Header line */
   // fgetln(F, &size);

   nread = fscanf(F, "%lg %lg\n", par_nz, par_nz+n+1);           /* z_0  n_0 */
   testErrorRetVA(nread!=2, io_file, "Error while reading n(z) file ('%s'): two doubles expected (first line)",
		  *err, __LINE__, NULL, name);

   for (i=1; i<n; i++) {
      nread = fscanf(F, "%lg %lg\n", par_nz+i+1, par_nz+n+i+1);  /* z_i  n_i */
   //printf("%g %g  %d\n", par_nz[i+1], par_nz[n+i+1], nread);
      testErrorRetVA(nread!=2, io_file, "Error while reading n(z) file ('%s'): two doubles expected",
		     *err, __LINE__, NULL, name);
   }

   nread= fscanf(F, "%lg %lg\n", par_nz+1, &dummy);            /* z_{n-1}=z_max  0 */
   testErrorRetVA(nread!=2, io_file, "Error while reading n(z) file ('%s'): two doubles expected (last line)",
		  *err, __LINE__, NULL, name);

   testErrorRetVA(par_nz[1]-par_nz[0]<MINIMUM_ZBIN_WIDTH, redshift_narrow_hist,
		  "Redshift histogram too narrow, required is z_max-z_min>=%g. Use nofz type 'single' instead",
		  *err, __LINE__, NULL, MINIMUM_ZBIN_WIDTH);

   return par_nz;
}

/* Determines the redshift type of file 'name' from its header  *
 * # snofz [sphotz] 						                            */
void get_nofz_t_file(const char *name, nofz_t *nofz, photz_t *photz, error **err)
{
   FILE *F;
   int n, j, s_length, only_ws;
   struct { char snofz[128], sphotz[128]; nofz_t nofz; photz_t photz; } tmp;
   char dummy[128];

   /* Read first line with nofz type */
   F = fopen_err(name, "r", err);
   forwardError(*err, __LINE__,);

   n = fscanf(F, "%s %s", dummy, tmp.snofz);
   testErrorRetVA(strcmp(dummy, "#")!=0, io_file,
		  "Wrong format of redshift info file '%s': First line has to be '# snofz_t [sphotz_t]'",
		  *err, __LINE__,, name);
   testError(n!=2, io_file, "Wrong format of redshift info file '%s': First line has to be '# snofz_t [sphotz_t]'", *err, __LINE__);

   n = fscanf(F, "%[^\n]", dummy);
   testErrorVA(n!=0 && n!=1, io_file,
      "Wrong format of redshift info file '%s': First line has to be '# snofz_t [sphotz_t]'", *err, __LINE__, name);

   fclose(F);

   /* Set photz error string to default if not indicated */
   s_length = strlen(dummy);
   only_ws = 1;
   while (j < s_length)  {
      if (dummy[j] != ' ' && dummy[j] != '#') {
         only_ws = 0; // False, found non-white-space character
      }
      j++;
   }
   if (only_ws == 1) {
      strcpy(tmp.sphotz, "photz_no\0");
   } else {
      /* Copy while removing white-space characters */
      for (j=0,n=0; j<s_length; j++) {
         if (dummy[j] != ' ') tmp.sphotz[n++] = dummy[j];
      }
      tmp.sphotz[n] = '\0';
   }


   STRING2ENUM(tmp.nofz, tmp.snofz, nofz_t, snofz_t, j, Nnofz_t, err);
   *nofz = tmp.nofz;

   STRING2ENUM(tmp.photz, tmp.sphotz, photz_t, sphotz_t, j, Nphotz_t, err);
   *photz = tmp.photz;
}

/* Sets on output number of parameters Nnz corresponding to file 'name' */
void Nnz_from_file(const char *name, int *Nnz, error **err)
{
   size_t Nrec, Nlines;
   nofz_t nofz;
   photz_t photz;
   double *res;

   get_nofz_t_file(name, &nofz, &photz, err);
   forwardError(*err, __LINE__,);

   /* Read complete file for record count */
   Nrec = 0;
   res = (double*)read_any_list_count(name, &Nrec, "%lg", sizeof(double), &Nlines, err);
   forwardError(*err, __LINE__,);
   free(res);

   switch (nofz) {
      case ludo : case ymmk : case ymmk0const : case cfhtlens : case single :
         /* All params, include zmin, zmax */
         *Nnz = Nrec;
         break;
      case hist :
         *Nnz = Nrec - 1;                   /* All entries - last (=0) */
         break;
      default :
         *err = addErrorVA(redshift_unknown, "Unknown redshift type %d", *err, __LINE__, nofz);
   }
}

/* Reads n(z) file(s) of any nofz_t type and returns the corresponding redshift structure, one or more bins */
redshift_t *init_redshift_from_files(const char **name, int Nzbin, error **err)
{
   redshift_t *res;
   int n, this_Nnz, Nnz_max;

   /* Get maximum number of redshift parameters */
   for (n=0,Nnz_max=-1; n<Nzbin; n++) {
      Nnz_from_file(name[n], &this_Nnz, err);
      forwardError(*err, __LINE__, NULL);
      if (Nnz_max<this_Nnz) Nnz_max = this_Nnz;
   }

   res = init_redshift_empty(Nzbin, Nnz_max, err);
   forwardError(*err, __LINE__, NULL);
   for (n=0; n<Nzbin; n++) {
      read_redshift_slice(res, n, name[n], err);
      forwardError(*err, __LINE__, NULL);
   }

   return res;
}

redshift_t *init_redshift_from_histogram_file(const char *name, error **err)
{
   redshift_t *self;
   self = init_redshift_from_files(&name, 1, err);
   forwardError(*err, __LINE__, NULL);
   return self;
}

/* Returns (true) z_min for redshift bin #n_bin */
double get_zmin(const redshift_t *redshift, int n_bin)
{
   double zmin;

   if (redshift->photz[n_bin] != photz_no) {
      zmin = 0.0;
   } else {
      zmin = redshift->par_nz[n_bin*redshift->Nnz_max + 0];
   }

   return zmin;
}

/* Returns (true) z_max for redshift bin #n_bin */
double get_zmax(const redshift_t *redshift, int n_bin)
{
   double zmax;
   if (redshift->photz[n_bin] != photz_no) {
      zmax = ZMAX;
   } else {
      zmax = redshift->par_nz[n_bin*redshift->Nnz_max + 1];
   }

   return zmax;
}

/* Returns minimum scale factor corresponding to (true) *
 * overall maximum redshift                             */
double get_amin(const redshift_t *redshift, error **err)
{
   int n;
   double amin, amintmp;

   for (n=0; n<redshift->Nzbin; n++) {
      if (redshift->photz[n] != photz_no) {
         amin = 1.0 / ZMAX - 1.0;
         return amin;
      }
   }

   for (n=0,amin=1.0; n<redshift->Nzbin; n++) {
      amintmp = 1.0/(1.0 + get_zmax(redshift, n));
      if (amintmp<amin) amin = amintmp;
   }

   testErrorRet(amin>=1.0, redshift_zmax, "zmax not valid (has to be different from zero)",
		*err, __LINE__, 0.0);

   /* To be on the safe side, return a slightly smaller scale factor. *
    * This was a problem for the highest bin (if narrow) of lensing tomography */
   if (amin>0.01) amin = amin-0.01;

   return amin;
}

/* Returns (photometric) z_min for redshift bin #n_bin (used to be get_zmin) */
double get_zmin_ph(const redshift_t *redshift, int n_bin)
{
   double zmin;

   zmin = redshift->par_nz[n_bin*redshift->Nnz_max + 0];

   return zmin;
}

/* Returns (photometric) z_max for redshift bin #n_bin (used to be get_zmax) */
double get_zmax_ph(const redshift_t *redshift, int n_bin)
{
   double zmax;

   zmax = redshift->par_nz[n_bin*redshift->Nnz_max + 1];

   return zmax;
}


/* Return the unnormalised redshift probability distribution */
double prob_unnorm(double z, void *intpar, sm2_error **err)
{
   double res;
   int n_par[Nnofz_t], n_par_photz[Nphotz_t];
   redshiftANDint* rANDi;
   redshift_t *self;
   int n_bin;
   size_t neval;
   int n, status;
   double epsabs, epsrel;
   double abserr;
   gsl_function F;
   params_photz_t params_photz;
   gsl_integration_workspace *w;

   rANDi = (redshiftANDint*)intpar;
   self  = rANDi->self;
   n_bin = rANDi->i;

   /* Number of parameters for all nofz and photz types */
   n_par[ludo] = n_par[jonben] = n_par[ymmk] = 5;
   n_par[ymmk0const] = n_par[cfhtlens] = 7;
   n_par[single] = 2;
   n_par[hist] = self->Nnz_max;

   n_par_photz[photz_no]    = 0;
   n_par_photz[photz_gauss] = 7;

   if (self->nofz[n_bin] != hist) {
      testErrorRetVA(self->Nnz[n_bin] != n_par[self->nofz[n_bin]] + n_par_photz[self->photz[n_bin]], redshift_nnz,
                  "Wrong number of redshift parameters %d, should be %d for nofz=%s(%d) and photz=%s(%d)",
                  *err, __LINE__, 0,
                  self->Nnz[n_bin], n_par[self->nofz[n_bin]], n_par_photz[self->photz[n_bin]],
                  snofz_t(self->nofz[n_bin]), self->nofz[n_bin], sphotz_t(self->photz[n_bin]), self->photz[n_bin]);
   }

   if (self->photz[n_bin] == photz_no) {

      if (z < get_zmin(self, n_bin) || z > get_zmax(self, n_bin)) return 0.0;

      res = prob_unnorm_nofz(z, n_bin, self, err);
      forwardError(*err, __LINE__, 0.0);

   } else {

      /* Convolve with photometric error distribution */

      params_photz.redshift = self;
      params_photz.n_bin    = n_bin;
      params_photz.z        = z;
      params_photz.Npz      = n_par_photz[self->photz[n_bin]];
      params_photz.err      = err;

      n          = 1000;
      epsabs     = 1.0e-6;
      epsrel     = 1.0e-3;
      w          = gsl_integration_workspace_alloc(n);
      F.function = &dpn_dz; 
      F.params   = (void*)&params_photz;

      status = gsl_integration_qng(&F, get_zmin_ph(self, n_bin), get_zmax_ph(self, n_bin), epsabs, epsrel, &res, &abserr, &neval);
      forwardError(*err, __LINE__, 0.0);

      testErrorVA(status != 0, redshift_integration, "GSL integration returned status %d, should be zero, something went wrong",
         *err, __LINE__, 0.0, status);
      //res = sm2_qromberg(dpn_dz, (void*)args, get_zmin_ph(self, n_bin), get_zmax_ph(self, n_bin), 1.0e-6, err);

      /* TODO Could speed up by re-using w */
      gsl_integration_workspace_free(w);
   }

   return res;
}

double dpn_dz(double zp, void *intpar)
{
   params_photz_t *params_photz;
   redshift_t *redshift;
   error **err;
   double res;
   int n_bin, Npz;

   params_photz = (params_photz_t *)intpar;
   redshift     = params_photz->redshift;
   n_bin        = params_photz->n_bin;
   Npz          = params_photz->Npz;
   err          = params_photz->err;

   res  = prob_photz(params_photz->z, zp, redshift->photz[n_bin],
            redshift->par_nz + n_bin*redshift->Nnz_max + redshift->Nnz[n_bin] - Npz, err);
   forwardError(*err, __LINE__, 0.0);
   res *= prob_unnorm_nofz(zp, params_photz->n_bin, redshift, err);
   forwardError(*err, __LINE__, 0.0);

   return res;
}

double gaussian(double z, double mu, double sigma)
{
   double res;

   res = exp(-0.5 * dsqr((z - mu)/sigma)) / sqrt(twopi) / sigma;

   return res; 
}

double prob_photz(double z, double zp, photz_t photz, const double *par_pz, error **err)
{
   double res;
   params_photz_gauss par;

   switch (photz) {

      case photz_gauss:
         par.sigma_z     = par_pz[0];
         par.z_bias      = par_pz[1];
         par.c_cal       = par_pz[2];
         par.f_out       = par_pz[3];
         par.sigma_z_out = par_pz[4];
         par.z_bias_out  = par_pz[5];
         par.c_cal_out   = par_pz[6];

         //printf("MKDEBUG prob_photz %g %g %g   %g   %g %g %g\n", par.sigma_z, par.z_bias, par.c_cal, par.f_out,
            //par.sigma_z_out, par.z_bias_out, par.c_cal_out);

         res = (1.0 - par.f_out) * gaussian(z, par.c_cal*zp - par.z_bias, par.sigma_z * (1 + zp))
               + par.f_out * gaussian(z, par.c_cal_out*zp - par.z_bias_out, par.sigma_z_out * (1 + zp));
         break;

      case photz_no:
         *err = addError(redshift_unknown, "Invalid photometric redshift error distribution 'photoz_no', would be delta function",
               *err, __LINE__);
         res = 0.0;
         break;

      default:
         *err = addErrorVA(redshift_unknown, "Unknown photometric redshift error distribution '%d'",
               *err, __LINE__, photz);
         res = 0.0;

   }

   return res;
}

/* Zero-Truncation for p <= exp(-TRUNC) */
#define TRUNC 100.0

/* Returns the 'true' unnormalized redshift probability distribution n(z) */
double prob_unnorm_nofz(double z, int n_bin, const redshift_t *self, error **err)
{
   double x, p, res, a, b, c, z0, zmin, zmax, ac, d, z1, z2;
   int n, m;

   switch (self->nofz[n_bin]) {
      case ludo :
         a  = self->par_nz[n_bin*self->Nnz_max+2];
         b  = self->par_nz[n_bin*self->Nnz_max+3];
         z0 = self->par_nz[n_bin*self->Nnz_max+4];

         testErrorRet(b<EPSILON, redshift_betap, "beta_p<0", *err, __LINE__, 0);

         x = z/z0;
         p = pow(x, b);

         if (p>TRUNC) {
            res = 0.0;
         } else {
            res = pow(x, a)*exp(-p);
         }
         break;

      case jonben :
         a = self->par_nz[n_bin*self->Nnz_max+2];
         b = self->par_nz[n_bin*self->Nnz_max+3];
         c = self->par_nz[n_bin*self->Nnz_max+4];
         res = pow(z, a)/(pow(z, b) + c);
         break;

      case ymmk :
         a = self->par_nz[n_bin*self->Nnz_max+2];
         b = self->par_nz[n_bin*self->Nnz_max+3];
         c = self->par_nz[n_bin*self->Nnz_max+4];
         res = (pow(z, a) + pow(z, a*b))/(pow(z, b) + c);
         break;

      case ymmk0const :
         a = self->par_nz[n_bin*self->Nnz_max+2];
         b = self->par_nz[n_bin*self->Nnz_max+3];
         c = self->par_nz[n_bin*self->Nnz_max+4];
         res = (pow(z, a) + pow(z, a*b))/(pow(z, b) + c);

         testErrorRet(self->par_nz[n_bin*self->Nnz_max+5]<0, redshift_parnzerror, "z0<0!",
               *err, __LINE__, 0);

         /* For z>z0 multiply n(z) with constant 0<c<1 */
         if (z>self->par_nz[n_bin*self->Nnz_max+5]) {
            if (self->par_nz[n_bin*self->Nnz_max+6]<=0) res = 0;
            else res *= self->par_nz[n_bin*self->Nnz_max+6];
         }
         break;

      case cfhtlens :
         z1 = self->par_nz[n_bin*self->Nnz_max+2];
         z2 = self->par_nz[n_bin*self->Nnz_max+3];
         ac = self->par_nz[n_bin*self->Nnz_max+4];
         b  = self->par_nz[n_bin*self->Nnz_max+5];
         d  = self->par_nz[n_bin*self->Nnz_max+6];
         res = ac * exp(-dsqr((z-z1) / b)) + exp(-dsqr((z-z2) / d));
         break;

      case hist :
         zmin = get_zmin_ph(self, n_bin);
         zmax = get_zmax_ph(self, n_bin);
         testErrorRetVA(z<zmin, redshift_zmax, "z=%g smaller than zmin=%g", *err, __LINE__, 0.0, z, zmin);
         testErrorRetVA(z>zmax, redshift_zmax, "z=%g larger than zmax=%g", *err, __LINE__, 0.0, z, zmax);

         if (self->Nnz[n_bin]==3) {

            /* Only one redshift bin */
            res = self->par_nz[n_bin*self->Nnz_max+2];

         } else {

            n = (self->Nnz[n_bin]-1)/2;

            if (z>=zmin && z<self->par_nz[n_bin*self->Nnz_max+2]) {
               /* z_0=z_min <= z < z_1: ---> n_0 */
               res = self->par_nz[n_bin*self->Nnz_max+n+1];
            } else {
               for (m=2; m<n; m++) {
                  if (z>=self->par_nz[n_bin*self->Nnz_max+m] && z<self->par_nz[n_bin*self->Nnz_max+m+1]) break;
               }
               /* z_m <= z < z_{m+1} ---> n_m */
               res = self->par_nz[n_bin*self->Nnz_max+n+m];
            }

         }
         break;

      case single:

         *err = addError(redshift_unknown, "Invalid redshift distribution 'single' in prob_unnorm",
               *err, __LINE__);
         res = 0.0;

      default :
         *err = addErrorVA(redshift_unknown, "Unknown redshift distribution %d (type nofz_t). "
               "See lensing.h for possible values", *err, __LINE__, self->nofz[n_bin]);
         res = 0.0;
   }

   return res;
}
#undef TRUNC

int change_prob(redshift_t *avant, redshift_t *apres)
{
   return change_redshift(avant, apres);
}

int change_redshift(redshift_t *avant, redshift_t *apres)
{
   int n, m;

   if (NCOEQ(avant, apres, Nzbin)) return 1;
   for (n=0; n<avant->Nzbin; n++) {
      if (NCOEQ(avant, apres, Nnz[n])) return 1;
      if (NCOEQ(avant, apres, nofz[n])) return 1;
      for (m=0; m<avant->Nnz[n]; m++) {
	      if (NCOCLOSE(avant, apres, par_nz[n*avant->Nnz_max+m])) return 1;
      }
   }
   return 0;
}

/* ============================================================ *
 * Redshift probability distribution of source galaxies.	       *
 * Note that a single redshift sheet (beta_p=0) is		          *
 * implemented in g_source.					                      *
 * ============================================================ */
double prob(redshift_t *self, double z, int n_bin, error **err)
{
   double res, zz;
   redshiftANDint rANDi;
   double zmin, zmax;

   zz = z*self->z_rescale[n_bin];

   testErrorRet(zz<0, ce_negative, "Negative redshift z", *err, __LINE__ , 0);
		
   rANDi.self = self;
   rANDi.i    = n_bin;

   if (self->nofz[n_bin] == single) {
         /* Note: This case should only appear if prob is called directly. The lensing projection
          * (g_source) with nofz=single is done without calling prob.
          */
      return 0.0;
   }

   if (self->prob_norm==NULL) {
      self->prob_norm = calloc_err(self->Nzbin, sizeof(double), err);
      forwardError(*err, __LINE__, 0.0);
   }

   if (self->prob_norm[n_bin] <= 0.0) {

      zmin = get_zmin(self, n_bin);
      zmax = get_zmax(self, n_bin);

      res = sm2_qromberg(prob_unnorm, (void*)&rANDi, zmin, zmax, 1.0e-6, err);
      //res = sm2_qromberg(prob_unnorm, (void*)&rANDi, 0.0, ZMAX, 1.0e-6, err);
      forwardError(*err,__LINE__, 0.0);

      testErrorRetVA(fabs(res)<1.0e-10, ce_infnan, "Int[ n(z) dz] in bin %d is zero, division by 0",
         *err, __LINE__, 0.0, n_bin);

      self->prob_norm[n_bin] = 1.0/res;
   }

   res = self->prob_norm[n_bin]*prob_unnorm(zz, (void*)&rANDi, err);
   forwardError(*err, __LINE__, 0);

   return res;
}

double int_for_zmean(double z, void *intpar, error **err)
{
   double res, p;
   redshiftANDint *rANDi;
   int n_bin;

   rANDi = (redshiftANDint*)intpar;
   n_bin = rANDi->i;

   p   = prob(rANDi->self, z, n_bin, err);
   forwardError(*err, __LINE__, 0);

   res = z*p;

   return res;	
}

#define EPS 1.0e-5
int change_zmean(redshift_t *avant, redshift_t *apres, error **err)
{
   double a_avant, a_apres;

   a_avant = 1.0 / (1.0 + zmean(avant, 0, err)); forwardError(*err, __LINE__, -1);
   a_apres = 1.0 / (1.0 + zmean(apres, 0, err)); forwardError(*err, __LINE__, -1);

   if (fabs(a_avant - a_apres) > EPS) return 1;
   return 0;
}

double zmean(redshift_t *self, int n_bin, error **err)
{
   double res, zmin, zmax;
   redshiftANDint rANDi;

   testErrorRetVA(n_bin<0, redshift_Nzbin, "Redshift bin (%d) smaller than zero",
		  *err, __LINE__, -1.0, n_bin);
   testErrorRetVA(n_bin>=self->Nzbin, redshift_Nzbin, "Redshift bin (%d) larger than maximum (Nzbin=%d)",
		  *err, __LINE__, -1.0, n_bin, self->Nzbin);

   /* Single lens plane */
   if (self->nofz[n_bin]==single)
     return get_zmin(self, n_bin);   /* zmin = zmax = z0 */

   rANDi.self = self;
   rANDi.i    = n_bin;
   zmin       = get_zmin(self, n_bin);
   zmax       = get_zmax(self, n_bin);
   res = sm2_qromberg(int_for_zmean, (void*)&rANDi, zmin, zmax, 1.0e-6, err);
   forwardError(*err, __LINE__, 0);

   return res;	
}

double zmedian(redshift_t *self, int n_bin, error **err)
{
   const double dz = 0.00001;
   double res, zmed;

   testErrorRetVA(n_bin<0, redshift_Nzbin, "Redshift bin (%d) smaller than zero",
		  *err, __LINE__, -1.0, n_bin);
   testErrorRetVA(n_bin>=self->Nzbin, redshift_Nzbin, "Redshift bin (%d) larger than maximum (Nzbin=%d)",
		  *err, __LINE__, -1.0, n_bin, self->Nzbin);

   if (self->nofz[n_bin] == single) {
      return get_zmin(self, n_bin);
   }

   res = zmed = 0;
   while (res<=0.5) {
      res += prob(self, zmed, n_bin, err) * dz;
      forwardError(*err, __LINE__, 0);
      zmed += dz;
   }

   return zmed;
}


/* Checks consistency of Nzcorr with the number of redshift bins n, and *
 * checks consistency of this number n and Nzbin.			*/
int get_and_check_Nzbin(int Nzcorr, int Nzbin, error **err)
{
   double dNzbin;
   int n;

   /* Solution to Nzcorr = Nzbin (Nzbin + 1) / 2 */
   dNzbin = 0.5*(sqrt(8.0*Nzcorr+1.0) - 1.0);
   testErrorRetVA(fabs(dNzbin-floor(dNzbin+0.5))>EPSILON, redshift_Nzbin,
		  "Number of z-correlations inconsistent, %d != Nz*(Nz+1)/2", *err, __LINE__, 0, Nzcorr);
   n = (int)floor(dNzbin+0.5);

   if (Nzbin>0) {
      testErrorRetVA(Nzbin!=n, redshift_Nzbin,
		     "Number of z bins in file (%d) inconsistent with input value (%d)",
		     *err, __LINE__, 0, n, Nzbin);
   }

   return n;
}

// fgetln
