/* ------------------------------------------------------------ *
 * getHODModel.c
 * Jean Coupon 2015
 *
 * Versions:
 * v1.0: first version (March. 2015)
 * v0.1: prototype version (Feb. 2013)
 * ------------------------------------------------------------ */

#include "getHODModel.h"

int main(int argc, char** argv)
{
  /* initialize */
  error *myerr = NULL, **err;
  err = &myerr;
  
  ConfigGetModel para;
  
  initPara(argc, argv, &para, err);
  quitOnError(*err, __LINE__, stderr);

  /* Output type */
  switch (para.type) {
  case COSMO :
    writeCosmo(para, err);
	quitOnError(*err, __LINE__, stderr);
    break;
  case HMF :
    writeHMF(para, err);
    quitOnError(*err, __LINE__, stderr);
    break;
  case DM :
    writeTPCF_DM(para, err);
    quitOnError(*err, __LINE__, stderr);
    break;
  case SMF :
    writeSMF(para, err);
    quitOnError(*err, __LINE__, stderr);
    break;
  case WTHETA : case WP : case XI :
    writeTPCF_GAL(para, err);
    quitOnError(*err, __LINE__, stderr);
    break;
  case GGL:
    writeGGL(para, err);
    quitOnError(*err, __LINE__, stderr);
    break;
  case MSMH :
    writeMSMH(para, err);
    quitOnError(*err, __LINE__, stderr);
    break;
  case SHMR :
    writeSHMR(para, err);
    quitOnError(*err, __LINE__, stderr);
    break;
  case TESTS:
    writeTESTS(para, err);
    quitOnError(*err, __LINE__, stderr);
    break;
  }
  

  return EXIT_SUCCESS;
}

/*----------------------------------------------------------------*
 * Main fonctions
 *----------------------------------------------------------------*/

void writeCosmo(ConfigGetModel para, error **err)
/* Write cosmological quantities */
{  
  /* open fileOut */
  FILE *fileOut =  fopen_err(para.fileOutName, "w", err);
  forwardError(*err, __LINE__,);
  
  dump_param(para.model->cosmo, fileOut);
  fprintf(fileOut, "#  Mean redshift       = %f\n",  zmean(para.model->redshift, 0, err));
  forwardError(*err, __LINE__,);
  fprintf(fileOut, "#  comoving volume V_c = %g (h^-3 Mpc^3 per solid angle)\n", vc(para.model, para.model->zmin, para.model->zmax, err));    
  forwardError(*err, __LINE__,);
  fprintf(fileOut, "#  D+                  = %f\n", D_plus(para.model->cosmo, para.a, 1, err));
  forwardError(*err, __LINE__,);
  
  fclose(fileOut);

}


void writeHMF(ConfigGetModel para, error **err){
  /* Computes the halo mass fonction */
  size_t i;
  double M;
  
  double *Nh     = malloc(para.nbins*sizeof(double)); /* HMF */
  double *log10M = malloc(para.nbins*sizeof(double)); /* log10M */
  
  for(i=0;i<para.nbins;i++){
    M          = exp(para.min + (double)i*para.delta + para.delta/2.0);
    log10M[i]  = log10(M);
    Nh[i]      = dn_dlnM_uf(M, para.model, para.a, err)*log(10.0);
    forwardError(*err, __LINE__,);
  }
  
  /* open fileOut */
  FILE *fileOut =  fopen_err(para.fileOutName, "w", err);
  forwardError(*err, __LINE__,); 
  
  /* write fileOut */
  fprintf(fileOut, "# Cosmological parameters:\n");
  dump_param(para.model->cosmo, fileOut);	
  /* Multi redshift bin is not supported in halomodel 
      dump_redshift(model->redshift, F, err);
  */
  fprintf(fileOut, "# Mean redshift = %f\n",  zmean(para.model->redshift, 0, err));
  fprintf(fileOut, "#\n# Dark matter halo mass function)\n");
  fprintf(fileOut, "# Mstar[h^-1 M_sun] HMF[h^3 Mpc^-3 dex^-1]\n"); 
  for(i=0;i<para.nbins;i++) {
    fprintf(fileOut,"%e %e\n", pow(10.0, log10M[i]), Nh[i]);
  }
  
  free(log10M);
  free(Nh);
  fclose(fileOut);
}





void writeTPCF_DM(ConfigGetModel para, error **err){
  /* Computes two-point correlation for dark matter and haloes */
  size_t i;
  
  /* initilization */
  double *r      = malloc(para.nbins*sizeof(double));
  double *theta  = malloc(para.nbins*sizeof(double));
  double *k      = malloc(para.nbins*sizeof(double));

  double *xi_dm  = malloc(para.nbins*sizeof(double));
  double *wt_dm  = malloc(para.nbins*sizeof(double));
  double *wp_dm  = malloc(para.nbins*sizeof(double));
  double *pk_dm  = malloc(para.nbins*sizeof(double));
  
  for(i=0;i<para.nbins;i++){
    r[i] = exp(para.min + (double)i*para.delta + para.delta/2.0);
	 theta[i] = exp(para.min + (double)i*para.delta + para.delta/2.0) * arcmin;
  }
  for(i=0;i<para.nbins;i++){
    k[i] = 2.0*pi/r[para.nbins-i-1];
  }
  
  xi_dm   = xi_P_NL(para.model, para.a, r, para.nbins, err);
  forwardError(*err, __LINE__,);
  wt_dm   = woftheta(para.model, pnl, theta, para.nbins, 0, 0, err);
  forwardError(*err, __LINE__,);
  wp_dm   = wp(para.model, pnl, r, para.nbins, para.model->pi_max, GG, err);
  forwardError(*err, __LINE__,);
  for(i=0;i<para.nbins;i++){
    pk_dm[i] = P_NL(para.model->cosmo, para.a, k[i], err);
    forwardError(*err, __LINE__,);
  }
  
  /* open fileOut */
  FILE *fileOut =  fopen_err(para.fileOutName, "w", err);
  forwardError(*err, __LINE__,);
  
  /* write fileOut */
  fprintf(fileOut, "# Cosmological parameters:\n");
  dump_param(para.model->cosmo, fileOut);	
  /* Multi redshift bin is not supported in halomodel 
      dump_redshift(model->redshift, F, err);
  */
  fprintf(fileOut, "# Mean redshift = %f\n",  zmean(para.model->redshift, 0, err));
  fprintf(fileOut, "# M_halo [h^-1 M_sun] = %5.2e (log10 = %5.2f)\n", pow(10.0,para.log10Mhalo), para.log10Mhalo); 
  fprintf(fileOut, "# bh(M_halo)          = %5.2f\n", halo_bias(para.model, pow(10.0,para.log10Mhalo), para.a, 1, err)); 
  fprintf(fileOut, "# pi_max [h^-1 Mpc]   = %5.2f\n", para.model->pi_max);
  fprintf(fileOut, "#\n# Dark matter two-point correlation function, xi(r), wp_dm(rp) and P_K(k)\n");
  fprintf(fileOut, "#r[h^-1 Mpc]  xi_dm(r)  theta[arcmin] w_dm(theta) wp_dm(r)[h^-1 Mpc]  k[h Mpc^-1]  P(k)[(h Mpc^-1)^3]\n"); 
  for(i=0;i<para.nbins;i++) {
    fprintf(fileOut,"%f %e %f %e %f %f %e\n", r[i], xi_dm[i], theta[i]/arcmin, wt_dm[i], wp_dm[i], k[i], pk_dm[i]);
  }

  fclose(fileOut);
}


void writeSMF(ConfigGetModel para, error **err){
  /* Computes the stellar mass fonction */
  size_t i;
  
  double *log10Mstellar = malloc(para.nbins*sizeof(double)); /* Mstar */
  
  for(i=0;i<para.nbins;i++)  log10Mstellar[i]  = log10(exp(para.min + (double)i*para.delta + para.delta/2.0));

  double *Ngal   = dNdlog10Mstellar(para.model, log10Mstellar, para.nbins, err);
  double *Ngal_c = dNdlog10Mstellar_c(para.model, log10Mstellar, para.nbins, err);
  double *Ngal_s = dNdlog10Mstellar_s(para.model, log10Mstellar, para.nbins, err);
  
  /* open fileOut */
  FILE *fileOut =  fopen_err(para.fileOutName, "w", err);
  forwardError(*err, __LINE__,);
  
  /* write fileOut */
  dump_param_hm(para.model, fileOut, err); quitOnError(*err, __LINE__, stderr);
  fprintf(fileOut, "#\n# Stellar mass function\n");
  fprintf(fileOut, "# Mstar in [h^-1 M_sun], SMF in [h^3 Mpc^-3 dex^-1]\n");
  fprintf(fileOut, "# Mstar      SMF          SMF_centrals SMF_satellites\n");
  for(i=0;i<para.nbins;i++) {
    fprintf(fileOut,"%e %e %e %e\n", pow(10.0, log10Mstellar[i]), Ngal[i], Ngal_c[i], Ngal_s[i]);
  }

  free(Ngal);
  free(Ngal_c);
  free(Ngal_s);
  fclose(fileOut);
}

void  writeTPCF_GAL(ConfigGetModel para, error **err)
/* Computes two-point correlation for galaxies */
{
  size_t i;
      
  /* initilization */
  double *y1hgcs = NULL, *y1hgss = NULL, *y2hg = NULL;
  double *x      = malloc(para.nbins*sizeof(double)); /* theta, r or rp */
  double *y      = malloc(para.nbins*sizeof(double)); /* w(theta), wp(rp) or xi(r) */
  
  for (i=0;i<para.nbins;i++) {
    x[i] = exp(para.min + (double)i*para.delta + para.delta/2.0);
  }

  switch (para.type) {
	  case WTHETA :
		  y1hgcs = woftheta(para.model, p1hgcs, x, para.nbins, 0, 0, err);  forwardError(*err, __LINE__,);
		  y1hgss = woftheta(para.model, p1hgss, x, para.nbins, 0, 0, err);  forwardError(*err, __LINE__,);
		  y2hg   = woftheta(para.model, p2hg,   x, para.nbins, 0, 0, err);  forwardError(*err, __LINE__,);
		  break;
	  case WP:
		  y1hgcs = wp(para.model, p1hgcs, x, para.nbins, para.model->pi_max, GG, err);  forwardError(*err, __LINE__,);
		  y1hgss = wp(para.model, p1hgss, x, para.nbins, para.model->pi_max, GG, err);  forwardError(*err, __LINE__,);
		  y2hg   = wp(para.model, p2hg,   x, para.nbins, para.model->pi_max, GG, err);  forwardError(*err, __LINE__,);
		  break;
	  case XI:
		  y1hgcs = xiofr(para.model, p1hgcs, x, para.nbins, GG, err);  forwardError(*err, __LINE__,);
		  y1hgss = xiofr(para.model, p1hgss, x, para.nbins, GG, err);  forwardError(*err, __LINE__,);
		  y2hg   = xiofr(para.model, p2hg,   x, para.nbins, GG, err);  forwardError(*err, __LINE__,);
		  break;
  }


  for(i=0;i<para.nbins;i++) y[i] = y1hgcs[i] + y1hgss[i] + y2hg[i]; 
  
  /* for the occupation function N(M) */
  double M, min = 10.0, max = 16.0, delta;
  delta = (max - min)/(double)para.nbins;


  /* open fileOut */
  FILE *fileOut =  fopen_err(para.fileOutName, "w", err);
  forwardError(*err, __LINE__,);
  double Mc, Ms;


  /* Print out cosmological parameters + derived quantities in fileOut */ 
  dump_param_hm(para.model, fileOut, err); quitOnError(*err, __LINE__, stderr);
  fprintf(fileOut, "#\n# Deduced parameters (all masses in [h^-1 M_sun])\n");
  fprintf(fileOut, "# galaxy density n_g                = %e (in h^-3 Mpc^-3)\n", ngal_den(para.model, para.a, logMmax, para.model->log10Mstar_min, para.model->log10Mstar_max, err));
  forwardError(*err, __LINE__,);
  fprintf(fileOut, "# average galaxy bias <b_g> (z=0)   = %f (%f)\n", av_gal_bias(para.model, para.a, err), av_gal_bias(para.model, 1.0, err));
  forwardError(*err, __LINE__,);
  if (para.model->hod==berwein02_hexcl || para.model->hod==berwein02 || para.model->hod==leauthaud11 || para.model->hod==coupon15){
  	fprintf(fileOut, "# satellite fraction frsat          = %f\n", av_frsat(para.model, para.a,err));
    forwardError(*err, __LINE__,);
  }
  if (para.model->hod==leauthaud11 || para.model->hod==coupon15){
  	fprintf(fileOut, "# mean stellar mass log10(<Mstar>)  = %f\n", log10(av_stellar_mass(para.model, para.a, err)));
    forwardError(*err, __LINE__,);
  }
  fprintf(fileOut, "#\n# Mean halo properties for central galaxies (masses in [h^-1 M_sun]) \n");
  double Mav = av_halo_mass_cen(para.model, para.a, err); forwardError(*err, __LINE__,);
  fprintf(fileOut, "# average halo mass log10(<Mh>)     = %f\n", log10(Mav));
  fprintf(fileOut, "# concentration(<Mh>)               = %f\n", concentration(para.model, Mav, para.a,err));
  forwardError(*err, __LINE__,);
  fprintf(fileOut, "# bh(<Mh>)                          = %f\n", halo_bias(para.model, Mav, para.a, 1, err)); 
  forwardError(*err, __LINE__,);
  fprintf(fileOut, "# Delta_h (wrt mean density)        = %f\n", Delta_vir(para.model, para.a));

  /* Write fileout */
  switch (para.type) {
  case WTHETA :
    fprintf(fileOut, "#\n# Two-point correlation function (w(theta))\n");
    fprintf(fileOut,"#theta(deg)  w        w_1hgcs      w_1hgss      w_2hg        M            n(M)         n(M)_c       n(M)_s   total_Mstar_cen      total_Mstar_sat\n"); 
    break;
  case WP :
    fprintf(fileOut, "#\n# Two-point correlation function wp(rp)\n");
    fprintf(fileOut,"#r(h^-1 Mpc) w(h^1 Mpc)  w_1hgcs   w_1hgss      w_2hg        M            n(M)         n(M)_c       n(M)_s   total_Mstar_cen      total_Mstar_sat\n"); 
    break;
  case XI :
    fprintf(fileOut, "#\n# Two-point correlation function xi(r)\n");
    fprintf(fileOut,"#r(h^-1 Mpc) xi          xi_1hgcs  xi_1hgss     xi_2hg       M            n(M)         n(M)_c       n(M)_s   total_Mstar_cen      total_Mstar_sat\n"); 
    break;
  }
  
  for(i=0;i<para.nbins;i++) {
    M  = pow(10.0, min + (double)i*delta);     
    Mc = Mstar_tot_c(para.model, M, para.model->log10Mstar_min, para.model->log10Mstar_max, err);
    forwardError(*err, __LINE__,);
    Ms = Mstar_tot_s(para.model, M, para.model->log10Mstar_min, para.model->log10Mstar_max, err);
    forwardError(*err, __LINE__,);

    fprintf(fileOut,"%f %e %e %e %e %e %e %e %e %e %e\n", 
	    x[i], y[i], y1hgcs[i], y1hgss[i], y2hg[i],
	    M, Ngal(para.model, M, para.model->log10Mstar_min, para.model->log10Mstar_max, err), 
	    Ngal_c(para.model, M, para.model->log10Mstar_min, para.model->log10Mstar_max, err), 
	    Ngal_s(para.model, M, para.model->log10Mstar_min, para.model->log10Mstar_max, err),
	    Mc, Ms);
	forwardError(*err, __LINE__,);

  }
  
  fclose(fileOut);

  free(y); free(x);
  free(y1hgcs); free(y1hgss); free(y2hg);
  
  



  return;
}



void  writeGGL(ConfigGetModel para, error **err)
/* Computes galaxy-galaxy lensing. Delta Sigma [h M_sun pc^-2]*/
{
  size_t i;
  double Mav;
  
  /* initilization */
  double *ystellar, *y1hgcs, *y1hgss, *y2hg;
  double *x      = malloc(para.nbins*sizeof(double)); /* rp */
  double *y      = malloc(para.nbins*sizeof(double)); /* DeltaSigma */
  
  for(i=0;i<para.nbins;i++) x[i] = exp(para.min + (double)i*para.delta + para.delta/2.0);

  ystellar = DeltaSigma(para.model, pstellar, x, para.nbins, err);  forwardError(*err, __LINE__,);
  y1hgcs   = DeltaSigma(para.model, p1hgcs,   x, para.nbins, err);  forwardError(*err, __LINE__,);
  y1hgss   = DeltaSigma(para.model, p1hgss,   x, para.nbins, err);  forwardError(*err, __LINE__,);
  y2hg     = DeltaSigma(para.model, p2hg,     x, para.nbins, err);  forwardError(*err, __LINE__,);
  
  for(i=0;i<para.nbins;i++) y[i] = ystellar[i] + y1hgcs[i] + y1hgss[i] + y2hg[i]; 
  
  /* open fileOut */
  FILE *fileOut =  fopen_err(para.fileOutName, "w", err);
  forwardError(*err, __LINE__,);
  
 /* Print out cosmological parameters + derived quantities in fileOut */ 
  dump_param_hm(para.model, fileOut, err); quitOnError(*err, __LINE__, stderr);
  fprintf(fileOut, "# Delta_h (wrt mean density)        = %f\n", Delta_vir(para.model, para.a));
  fprintf(fileOut, "#\n# Deduced parameters for the full galaxy sample\n");
  fprintf(fileOut, "# galaxy density n_g                = %e (in h^-3 Mpc^-3)\n", ngal_den(para.model, para.a, logMmax, para.model->log10Mstar_min, para.model->log10Mstar_max, err));
  forwardError(*err, __LINE__,);
  fprintf(fileOut, "# average galaxy bias <b_g> (z=0)   = %f (%f)\n", av_gal_bias(para.model, para.a, err), av_gal_bias(para.model, 1.0, err));
  forwardError(*err, __LINE__,);
  if (para.model->hod==berwein02_hexcl || para.model->hod==berwein02 || para.model->hod==leauthaud11 || para.model->hod==coupon15){
  	fprintf(fileOut, "# satellite fraction frsat          = %f\n", av_frsat(para.model, para.a,err));
    forwardError(*err, __LINE__,);
  }
  if (para.model->hod==leauthaud11 || para.model->hod==coupon15){
  	fprintf(fileOut, "# mean stellar mass log10(<Mstar>)  = %f\n", log10(av_stellar_mass(para.model, para.a, err)));
    forwardError(*err, __LINE__,);
  }
  
  fprintf(fileOut, "#\n# Mean halo properties for central galaxies (masses in [h^-1 M_sun]) \n");

  Mav = av_halo_mass_cen(para.model, para.a, err); forwardError(*err, __LINE__,);
  fprintf(fileOut, "# average halo mass log10(<Mh_cen>) = %f\n", log10(Mav));
  fprintf(fileOut, "# concentration(<Mh_cen>)           = %f\n", concentration(para.model, Mav, para.a,err));
  forwardError(*err, __LINE__,);
  fprintf(fileOut, "# bh(<Mh_cen>)                      = %f\n", halo_bias(para.model, Mav, para.a, 1, err));
  forwardError(*err, __LINE__,);
  
  fprintf(fileOut, "#\n# Mean halo properties for satellite galaxies (masses in [h^-1 M_sun]) \n");

  Mav = av_halo_mass_sat(para.model, para.a, err); forwardError(*err, __LINE__,);
  fprintf(fileOut, "# average halo mass log10(<Mh_sat>) = %f\n", log10(Mav));
  fprintf(fileOut, "# concentration(<Mh_sat>)           = %f\n", concentration(para.model, Mav, para.a,err));
  forwardError(*err, __LINE__,);
  fprintf(fileOut, "# bh(<Mh_sat>)                      = %f\n", halo_bias(para.model, Mav, para.a, 1, err));
  forwardError(*err, __LINE__,);
  
  
  /* Write fileout */
  fprintf(fileOut, "#\n# Delta Sigma (r)\n");
  fprintf(fileOut, "#r(h^-1 Mpc) Dsig[h M_sun pc^-2] stellar  1hgcs   1hgss          2hg\n"); 
  for(i=0;i<para.nbins;i++) {
    fprintf(fileOut,"%f %e %e %e %e %e\n", x[i], y[i], ystellar[i], y1hgcs[i], y1hgss[i], y2hg[i]);
  }
  
  fclose(fileOut);
  
  free(y); free(x);
  free(ystellar); free(y1hgcs); free(y1hgss); free(y2hg);
  
  return;
}

void writeMSMH(ConfigGetModel para, error **err){
  /* Computes the Mstar-Mh relation (MSMH)
   	 Mstar_given_Mh = SH
     Mh_given_Mstar = HS
  */
  
  size_t i;
  
  if(strcmp(para.configpmcFileInName,"")){
  // TO DO
  }else if(strcmp(para.configmcmcFileInName,"")){
  // TO DO   
  }else{

    double *SH_log10Mstar = malloc(para.nbins*sizeof(double));
    double *SH_log10Mh    = malloc(para.nbins*sizeof(double));               
    double *HS_log10Mstar = malloc(para.nbins*sizeof(double));               
    double *HS_log10Mh    = malloc(para.nbins*sizeof(double));
    

    /* The range (para.min,para.max) is given for the halo masses. Below is the corresponding 
       range for the stellar masses. */
    double logMstar_min = MAX(9.0*log(10), log10_fSHMR(para.min/log(10), para.model)*log(10));
    double logMstar_max = log10_fSHMR(para.max/log(10), para.model)*log(10); 
  	double deltaMstar   = (logMstar_max - logMstar_min)/(double)para.nbins;

    /* loop over bins */
    for(i=0;i<para.nbins;i++){
		SH_log10Mh[i]        = log10(exp(para.min + (double)i*para.delta + para.delta/2.0));
		SH_log10Mstar[i]     = log10_fSHMR(SH_log10Mh[i], para.model);
		HS_log10Mstar[i]     = log10(exp(logMstar_min + (double)i*deltaMstar + deltaMstar/2.0));
		HS_log10Mh[i]        = log10(av_Mh_given_Mstar(para.model, HS_log10Mstar[i], para.a, err));
    }

    /* open fileOut */
    FILE *fileOut =  fopen_err(para.fileOutName, "w", err);
    forwardError(*err, __LINE__,); 
    
    /* write fileOut */
	dump_param_hm(para.model, fileOut, err); quitOnError(*err, __LINE__, stderr);
  	fprintf(fileOut, "#\n# Mstar-Mh relation (MSMH)\n");
    fprintf(fileOut, "# Mh      Mstar|Mh           Mstar      Mh|Mstar\n");

    for(i=0;i<para.nbins;i++) {
      fprintf(fileOut,"%e %e %e %e\n", pow(10.0, SH_log10Mh[i]), pow(10.0, SH_log10Mstar[i]), pow(10.0, HS_log10Mstar[i]), pow(10.0, HS_log10Mh[i]));
    }

    fclose(fileOut);
        
  	free(SH_log10Mstar);
  	free(SH_log10Mh);
  	free(HS_log10Mstar);
  	free(HS_log10Mh);
  }
}

void writeSHMR(ConfigGetModel para, error **err){
  /* Computes the stellar-to-halo mass ratio (SHMR)
   */
  
  size_t i;
  
  if(strcmp(para.configpmcFileInName,"")){
  // TO DO
  }else if(strcmp(para.configmcmcFileInName,"")){
  // TO DO   
  }else{
  
    double *Mh       = malloc(para.nbins*sizeof(double)); /* Mhalo */
    double *Mstar_c  = malloc(para.nbins*sizeof(double)); /* Total Mstar central gals */
    double *Mstar_s  = malloc(para.nbins*sizeof(double)); /* Total Mstar satellite gals */
    double *Mstar_t  = malloc(para.nbins*sizeof(double)); /* Total Mstar satellite gals */

    /* loop over bins */
    for(i=0;i<para.nbins;i++){
		Mh[i]      = exp(para.min + (double)i*para.delta + para.delta/2.0);
		Mstar_c[i] = Mstar_tot_c(para.model, Mh[i], para.model->log10Mstar_min, para.model->log10Mstar_max, err)/Mh[i];
		Mstar_s[i] = Mstar_tot_s(para.model, Mh[i], para.model->log10Mstar_min, para.model->log10Mstar_max, err)/Mh[i];
		Mstar_t[i] = Mstar_c[i] + Mstar_s[i];
    }

    /* open fileOut */
    FILE *fileOut =  fopen_err(para.fileOutName, "w", err);
    forwardError(*err, __LINE__,); 
    
    /* write fileOut */
	dump_param_hm(para.model, fileOut, err); quitOnError(*err, __LINE__, stderr);
  	fprintf(fileOut, "#\n# Stellar-to-halo mass ratio (SHMR)\n");
    fprintf(fileOut, "# Mh      Mstar_t           Mstar_c      Mstar_s\n");

    for(i=0;i<para.nbins;i++) {
      fprintf(fileOut,"%e %e %e %e\n", Mh[i], Mstar_t[i], Mstar_c[i], Mstar_s[i]);
    }

    fclose(fileOut);
        
  	free(Mh);
  	free(Mstar_c);
  	free(Mstar_s);
  	free(Mstar_t);

  }
}

void writeTESTS(ConfigGetModel para, error **err){
  /* This function is used for tests only */
  size_t i;
  
  para.nbins = 100;
  para.min   = log(1.0e-3);
  para.max   = log(0.8e2);
  para.delta = (para.max - para.min)/(double)para.nbins;
  
  double *x     = malloc(para.nbins*sizeof(double)); /* r */
  double *y     = malloc(para.nbins*sizeof(double)); 
  
  double z = para.z;
  double como_to_phys = 1.0+z;
  double a = 1.0/(1.0+z);

  for(i=0;i<para.nbins;i++){
    x[i] = exp(para.min + (double)i*para.delta + para.delta/2.0);
    //y[i] = DeltaSigma_WB2000(para.model, x[i], a, pow(10.0, 12.0)*0.72,  7.093888, err);
    y[i] = DeltaSigma_WB2000(para.model, x[i], a, pow(10.0, 14.17)*0.72, 3.59, err);
    forwardError(*err, __LINE__,);
  }
  
  double *y_DS = DeltaSigma(para.model, p1hgcs,   x, para.nbins, err);  forwardError(*err, __LINE__,);
  

  /* open fileOut */
  FILE *fileOut =  fopen_err(para.fileOutName, "w", err);
  forwardError(*err, __LINE__,); 

  for(i=0;i<para.nbins;i++) {
    fprintf(fileOut,"%e %e %e %e\n", x[i]/como_to_phys, y[i]*dsqr(como_to_phys),x[i], y_DS[i]);
    //fprintf(fileOut,"%e %e\n", x[i], y[i]);
  }
  
  free(x);
  free(y);
  return;
  
}


/*----------------------------------------------------------------*
 * Utils
 *----------------------------------------------------------------*/

void initPara(int argc, char **argv, ConfigGetModel *para, error **err){
  /* Function to initialize parameters for halolpot */
  size_t i, Ncol, rangeIsDefined = 0;
  char list[NFIELD*NCHAR];
  
  /* initialisation */
  strcpy(MYNAME, "getHODModel");
  strcpy(para->fileOutName,          "model.out");
  strcpy(para->configpmcFileInName,  "");
  strcpy(para->configmcmcFileInName, "");
  strcpy(para->sampleFileInName,     "");

  para->type        = COSMO;
  para->nbins       = 20;
  para->log10Mhalo  = 13.0;

  if(argc < 2 || !strcmp(argv[1],"-h") || !strcmp(argv[1],"--help")) {
    
    fprintf(stderr,"\n\n\
                          G E T H O D M O D E L\n\n\
      Program to compute cosmological and HOD quantities. Version 1.0\n\n\
To compute quantities:\n\
Usage:  %s halomodel.par -t TYPE [OPTIONS] [-para PARAMETERS]\n\n\
To compute quantities with error bars from pmc or mcmc run (only for -t [msmh,shmr]) :\n\
Usage:  %s halomodel.par -t TYPE [OPTIONS] [-para PARAMETERS] -c_[p,mc]mc CONFIG_FILE -p [chain,PMCSIM]_FILE \n", MYNAME, MYNAME);
    fprintf(stderr, "\nOPTIONS:\n");
    fprintf(stderr, "  -o OUT                 Output file name\n");
    fprintf(stderr, "  -t TYPE                Output type, TYPE in [cosmo,dm,wtheta,wp,xi,smf,shmr,ggl,hmf], default: cosmo\n");
    fprintf(stderr, "  -nbins                 Number of bins (logarithmic intervals)\n");
    fprintf(stderr, "  -range                 Range: min,max\n");
    fprintf(stderr, "  -para                  Parameter set, see below (supersedes halomodel.par)\n");
    fprintf(stderr, "  -c_pmc  CONFIG_FILE    pmc config file\n");
    fprintf(stderr, "  -c_mcmc CONFIG_FILE    mcmc config file\n");
    fprintf(stderr, "  -p [chain,PMCSIM]_FILE mcmc chain or pmc output file to compute model confidence limits\n");
    fprintf(stderr, "  -h                     This message\n");
    fprintf(stderr,"\nTYPE : PARAMETERS:\n");
    fprintf(stderr,"  cosmo  : no para\n");
    fprintf(stderr,"  hmf    : no para\n");
    fprintf(stderr,"  dm     : log10Mhalo,pi_max\n");
    fprintf(stderr,"\nIf berwein02 or berwein02_hexcl (hamana04 model is no longer maintained):\n");
    fprintf(stderr,"  xi,wtheta : log10M_min,log10M1,log10M0,sigma_log_M,alpha_halo,eta\n");
    fprintf(stderr,"  wp        : log10M_min,log10M1,log10M0,sigma_log_M,alpha_halo,eta,pi_max\n");
    fprintf(stderr,"  ggl       : log10M_min,log10M1,log10M0,sigma_log_M,alpha_halo,eta,pi_max,log10Mstar_min,log10Mstar_max\n");
    fprintf(stderr,"\nIf leauthaud11:\n");
    fprintf(stderr,"  msmh      : log10M1,log10Mstar0,beta,delta,gamma\n");
    fprintf(stderr,"  xi,wtheta : log10M1,log10Mstar0,beta,delta,gamma,log10Mstar_min,log10Mstar_max,sigma_log_M,B_cut,B_sat,beta_cut,beta_sat,alpha,fcen1,fcen2\n");
    fprintf(stderr,"  smf       : log10M1,log10Mstar0,beta,delta,gamma,log10Mstar_min,log10Mstar_max,sigma_log_M,B_cut,B_sat,beta_cut,beta_sat,alpha,fcen1,fcen2\n");
    fprintf(stderr,"  shmr      : log10M1,log10Mstar0,beta,delta,gamma,log10Mstar_min,log10Mstar_max,sigma_log_M,B_cut,B_sat,beta_cut,beta_sat,alpha,fcen1,fcen2\n");
    fprintf(stderr,"  ggl,wp    : log10M1,log10Mstar0,beta,delta,gamma,log10Mstar_min,log10Mstar_max,sigma_log_M,B_cut,B_sat,beta_cut,beta_sat,alpha,fcen1,fcen2,pi_max\n");
    fprintf(stderr,"\nIf coupon15:\n");
    fprintf(stderr,"  msmh      : log10M1,log10Mstar0,beta,delta,gamma\n");
    fprintf(stderr,"  xi,wtheta : log10M1,log10Mstar0,beta,delta,gamma,log10Mstar_min,log10Mstar_max,sigma_log_M,B_cut,B_sat,beta_cut,beta_sat,alpha,fcen1,fcen2\n");
    fprintf(stderr,"  smf       : log10M1,log10Mstar0,beta,delta,gamma,log10Mstar_min,log10Mstar_max,sigma_log_M,B_cut,B_sat,beta_cut,beta_sat,alpha,fcen1,fcen2\n");
    fprintf(stderr,"  shmr      : log10M1,log10Mstar0,beta,delta,gamma,log10Mstar_min,log10Mstar_max,sigma_log_M,B_cut,B_sat,beta_cut,beta_sat,alpha,fcen1,fcen2\n");
    fprintf(stderr,"  ggl,wp    : log10M1,log10Mstar0,beta,delta,gamma,log10Mstar_min,log10Mstar_max,sigma_log_M,B_cut,B_sat,beta_cut,beta_sat,alpha,fcen1,fcen2,pi_max\n");
    fprintf(stderr,"\nInput masses must be in log10 with unit of h^-1 M_{sol} and pi_max in Mpc.\n");
    fprintf(stderr,"All quantities are computed at mean redshift or single value from nzfile.\n");
    fprintf(stderr,"nzfile may of any type as describe in the doc, except for\n\
-t wtheta for which nzfile must be type \"hist\" (see cosmo_pmc documentation).\n");
    exit(EXIT_FAILURE);
  }
  
  /* read over the command line options */
  for(i=0;i<argc;i++){                
    if(!strcmp(argv[i],"-o"))       strcpy(para->fileOutName, argv[i+1]);
    if(!strcmp(argv[i],"-c_pmc"))   strcpy(para->configpmcFileInName, argv[i+1]);
    if(!strcmp(argv[i],"-c_mcmc"))  strcpy(para->configmcmcFileInName, argv[i+1]);
    if(!strcmp(argv[i],"-p"))       strcpy(para->sampleFileInName, argv[i+1]);

    if(!strcmp(argv[i],"-t")) {
      if(!strcmp(argv[i+1],"cosmo"))  para->type = COSMO;
      if(!strcmp(argv[i+1],"wtheta")) para->type = WTHETA;
      if(!strcmp(argv[i+1],"wp"))     para->type = WP;
      if(!strcmp(argv[i+1],"xi"))     para->type = XI;
      if(!strcmp(argv[i+1],"dm"))     para->type = DM;
      if(!strcmp(argv[i+1],"msmh"))   para->type = MSMH;
      if(!strcmp(argv[i+1],"smf"))    para->type = SMF;
      if(!strcmp(argv[i+1],"ggl"))    para->type = GGL;
      if(!strcmp(argv[i+1],"shmr"))   para->type = SHMR;
      if(!strcmp(argv[i+1],"hmf"))    para->type = HMF;
      if(!strcmp(argv[i+1],"tests"))  para->type = TESTS;
    }
    if(!strcmp(argv[i], "-range")){
      getStrings(argv[i+1], list, ",", &Ncol, err);
      forwardError(*err, __LINE__,);
      para->min = getDoubleValue(list, 1);
      para->max = getDoubleValue(list, 2);
      rangeIsDefined = 1;
    }
    if(!strcmp(argv[i], "-nbins")) para->nbins = atoi(argv[i+1]);
  }
  


  /* if "range" is not set in the command line, set default values */
  if(!rangeIsDefined){
    switch (para->type) {
    case COSMO  :   para->min = 0.0;    para->max = 10.0;   break;
    case WTHETA :   para->min = 0.001;  para->max = 10.0;   break;
    case WP     :   para->min = 0.01;   para->max = 60.0;   break;
    case XI     :   para->min = 0.01;   para->max = 60.0;   break;
    case DM     :   para->min = 0.01;   para->max = 60.0;   break;
    case MSMH   :   para->min = 10.0;   para->max = 16.0;   break;
    case SHMR   :   para->min = 10.0;   para->max = 16.0;   break;
    case SMF    :   para->min = 8.0;    para->max = 12.0;   break;
    case GGL    :   para->min = 0.001;  para->max = 60.0;   break;
    case HMF    :   para->min = 8.0;    para->max = 16.0;   break;
    }
  }
  
  /* cosmological model read from halomodel.par */
  FILE *cosmoFileIn = fopen_err(argv[1], "r", err);                   forwardError(*err, __LINE__,);
  read_cosmological_parameters_hm(&(para->model), cosmoFileIn, err);  forwardError(*err, __LINE__,);
  fclose(cosmoFileIn);
  
  /* mean redshift */
  para->z = zmean(para->model->redshift, 0, err);
  forwardError(*err, __LINE__,);
  
  /* set new parameters given from command line */
  for(i=0;i<argc;i++){
    if(!strcmp(argv[i], "-para")){
      	getStrings(argv[i+1], list, ",", &Ncol, err);
      	forwardError(*err, __LINE__,);
    	switch (para->type) {
      		case COSMO :  case HMF :
      			/* no parameters needed */
	  			break;

	  		case DM :
				para->log10Mhalo         = getDoubleValue(list, 1);
	   			para->model->pi_max      = getDoubleValue(list, 2);
				/* Set minimum value for a. P_NL crashes otherwise */
				para->model->cosmo->a_min = 0.001;
				break;
      
		    case WTHETA : case XI : case SHMR : case SMF : case MSMH : case GGL : case WP :
				switch (para->model->hod) {
					case berwein02 : case berwein02_hexcl: case hamana04 :
			  			para->model->log10M_min     = getDoubleValue(list, 1);
			  			para->model->log10M1        = getDoubleValue(list, 2);
			  			para->model->log10M0        = getDoubleValue(list, 3);
			  			para->model->sigma_log_M    = getDoubleValue(list, 4);
			  			para->model->alpha          = getDoubleValue(list, 5);
			  			para->model->eta            = getDoubleValue(list, 6);
				  		if(para->type == GGL || para->type == WP){
				  			para->model->pi_max         = getDoubleValue(list, 7);
					  		if(para->type == GGL){
				  				para->model->log10Mstar_min = getDoubleValue(list, 8);	
				  				para->model->log10Mstar_max = getDoubleValue(list, 9);
							}
			  			}
			  			break;
					case leauthaud11:
		      			para->model->log10M1        = getDoubleValue(list, 1);
			  			para->model->log10Mstar0    = getDoubleValue(list, 2);
			  			para->model->beta           = getDoubleValue(list, 3);
			  			para->model->delta          = getDoubleValue(list, 4);
			  			para->model->gamma          = getDoubleValue(list, 5);
			  			if(para->type != MSMH){
				  			para->model->log10Mstar_min = getDoubleValue(list, 6);	
				  			para->model->log10Mstar_max = getDoubleValue(list, 7);
				  			para->model->sigma_log_M    = getDoubleValue(list, 8);
				  			para->model->B_cut          = getDoubleValue(list, 9);
				  			para->model->B_sat          = getDoubleValue(list, 10);
				  			para->model->beta_cut       = getDoubleValue(list, 11);
				  			para->model->beta_sat       = getDoubleValue(list, 12);
				  			para->model->alpha          = getDoubleValue(list, 13);
				  			para->model->fcen1          = getDoubleValue(list, 14);
				  			para->model->fcen2          = getDoubleValue(list, 15);
				  			if(para->type == GGL || para->type == WP){
						  		para->model->pi_max         = getDoubleValue(list, 16);
							}
				  		}
			  			break;
		    	case coupon15:
						para->model->log10M1        = getDoubleValue(list, 1);
			  			para->model->log10Mstar0    = getDoubleValue(list, 2);
			  			para->model->beta           = getDoubleValue(list, 3);
			  			para->model->delta          = getDoubleValue(list, 4);
			  			para->model->gamma          = getDoubleValue(list, 5);
			  			if(para->type != MSMH){
				  			para->model->log10Mstar_min = getDoubleValue(list, 6);	
				  			para->model->log10Mstar_max = getDoubleValue(list, 7);
				  			para->model->sigma_log_M    = getDoubleValue(list, 8);
				  			para->model->B_cut          = getDoubleValue(list, 9);
				  			para->model->B_sat          = getDoubleValue(list, 10);
				  			para->model->beta_cut       = getDoubleValue(list, 11);
				  			para->model->beta_sat       = getDoubleValue(list, 12);
				  			para->model->alpha          = getDoubleValue(list, 13);
				  			para->model->fcen1          = getDoubleValue(list, 14);
				  			para->model->fcen2          = getDoubleValue(list, 15);
				  			if(para->type == GGL || para->type == WP){
						  		para->model->pi_max         = getDoubleValue(list, 16);
							}
				  		}
				  		break;
				case hod_none:
						*err = addError(hm_hodtype, "Invalid HOD type 'hod_none'", *err, __LINE__);
						return;
				}
				break;
    		}
	    }
  	}

  
  /* wtheta must have a broad redshift distribution */
  testErrorExit(para->type == WTHETA && para->model->redshift->nofz[0] != hist , GM_nz_type,  
	       "For WTHETA, nz type must be \"hist\"", *err, __LINE__);
    
  /* rmax should be smaller than pi_max */
  testErrorExit((para->type == WP || para->type == DM) && para->max > para->model->pi_max, GM_pi_max,  
	       "rmax should be smaller than pi_max", *err, __LINE__);
  
  /* shmr only for leauthaud11 and coupon15*/
  testErrorExit(para->type == MSMH && (para->model->hod != leauthaud11 && para->model->hod != coupon15), GM_SHMR,  
	       "msmh only available for leauthaud11 and coupon15 models", *err, __LINE__);

  /* tshmr only for leauthaud11 */
  testErrorExit(para->type == SHMR && (para->model->hod != leauthaud11 && para->model->hod != coupon15), GM_SHMR,  
	       "shmr only available for leauthaud11 and coupon15 models", *err, __LINE__);
  
  /* smf only for leauthaud11 */
  testErrorExit(para->type == SMF && para->model->hod != leauthaud11 && para->model->hod != coupon15, GM_SMF,  
	       "smf only available for leauthaud11 and coupon15 models", *err, __LINE__);
  
  testErrorRet((para->model->hod == leauthaud11 || para->model->hod == coupon15) && para->model->log10Mstar_min < 0, GM_MSTAR_MAX,
                "Minimum stellar mass cannot be set to a negative value (undef) for leauthaud11 and coupon15 models",
                *err, __LINE__,);
  testErrorRet((para->model->hod == leauthaud11 || para->model->hod == coupon15) && para->model->log10Mstar_max < 0, GM_MSTAR_MAX,
                "Maximum stellar mass cannot be set to a negative value (undef) for leauthaud11 and coupon15 models. If threshold sample, simply put log10Mstar_max=12.5",
                *err, __LINE__,);

  testErrorRet(para->model->log10Mstar_min > para->model->log10Mstar_max && para->model->log10Mstar_max > 0, GM_MSTAR_MAX,
                "Maximum has to be larger than minimum stellar mass",
                *err, __LINE__,);

  
  /* logscale */
  if(para->type == SMF || para->type == MSMH || para->type == SHMR || para->type == HMF){
    para->min *= log(10.0); /* base ln */
    para->max *= log(10.0); /* base ln */
  }else{
    para->min = log(para->min);
    para->max = log(para->max);
  }

  /* scale factor at redshift z */
  para->a = 1.0/(1.0 + para->z);
  
  /* binsize */
  para->delta = (para->max - para->min)/(double)para->nbins;
  
}

#define NFIELD_LARGE 1000
#define NCHAR_LARGE  100

long fileSize(FILE *fileIn){
  char line[NFIELD_LARGE*NCHAR_LARGE];
  
  /* get size of file and allocate data */
  long N = 0;
  
  rewind(fileIn);
  while(fgets(line, NFIELD_LARGE*NCHAR_LARGE, fileIn) != NULL)
    if(*(line) != '#' && *(line) != '\0' && *(line) != '\n') N++;
  rewind(fileIn);
 
  return N;
}

#undef NFIELD_LARGE
#undef NCHAR_LARGE
