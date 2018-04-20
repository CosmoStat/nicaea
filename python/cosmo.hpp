/*
 * <one line to give the program's name and a brief idea of what it does.>
 * Copyright (C) 2015  Francois <email>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef COSMO_HPP
#define COSMO_HPP

#include "cosmo.h"

using namespace nicaea;

class cosmology: public cosmo
{
    void init( double OMEGAM, double OMEGADE, double W0_DE, double W1_DE,
               double *W_POLY_DE, int N_POLY_DE,
               double H100, double OMEGAB, double OMEGANUMASS,
               double NEFFNUMASS, double NORM, double NSPEC,
               nonlinear_t NONLINEAR, transfer_t TRANSFER, growth_t GROWTH,
               de_param_t DEPARAM, norm_t normmode, double AMIN);
public:

    /* Parameters are:
    Om Od w0 w1 h Ob Onu Neffnu s8 ns
    nonlin transfer growth deparam norm amin */
    cosmology():myerr(NULL) {
        init(0.25, 0.75, -1.0, 0.0, NULL, 0, 0.70, 0.044, 0.0, 0.0, 0.80, 0.96,
             smith03, eisenhu, growth_de, linder, norm_s8, 0.0);
    }

    ~cosmology() {
        if (de_param == poly_DE) {
            free(w_poly_de);
        }
        del_interTable(&linearGrowth);
        del_interTable(&transferFct);
        del_interTable(&transferBE);
        del_interTable2D(&P_NL);
        del_interTable(&slope);
        del_interTable(&w);
        delete *err;
    }

    /* ============================================================ *
     * Cosmology.                           *
     * ============================================================ */

    /* Homogeneous Universe */
    double da_dtau(double a) {
        return ::da_dtau(this,a,err);
    }
    double w_de(double a) {
        return ::w_de(this,a,err);
    }
    double f_de(double a) {
        return ::f_de(this,a,err);
    }
    double Esqr(double a, int wOmegar) {
        return ::Esqr(this,a, wOmegar,err);
    }
    double Omega_m_a(double a, double Esqrpre) {
        return ::Omega_m_a(this,a,Esqrpre,err);
    }
    double Omega_de_a(double a, double Esqrpre) {
        return ::Omega_de_a(this,a,Esqrpre,err);
    }
    double w_nu_mass(double a) {
        return ::w_nu_mass(this,a);
    }

    /* Geometry, distances */
    double w_a(double a, int wOmegar) {
        return ::w(this,a,wOmegar,err);
    }
    double dwoverda(double a) {
        return ::dwoverda(this,a,err);
    }
    double drdz(double a) {
        return ::drdz(this,a,err);
    }
    double dvdz(double a) {
        return ::dvdz(this,a,err);
    }
    double f_K(double wr) {
        return ::f_K(this,wr,err);
    }
    double D_a(double a) {
        return ::D_a(this,a,err);
    }
    double D_a12(double a1, double a2) {
        return ::D_a12(this,a1,a2,err);
    }
    double D_lum(double a) {
        return ::D_lum(this,a,err);
    }


//
//     /* Growth factor */
//     double D_plus(cosmo*, double a, int normalised, error **err);
//     double g(cosmo*, double a);
//
//     /* Transfer function */
//     double r_sound(cosmo *model);
//     double r_sound_integral(cosmo *self, double a, error **err);
//     double r_sound_drag_fit(cosmo *model, error **err);
//     double r_sound_drag_analytical(cosmo *self, error **err);
//     double k_silk(const cosmo *model);
//     double ratio_b_gamma(cosmo *self, double a);
//     double z_drag(cosmo *self);
//     double T_tilde(const cosmo *self, double k, double alpha_c, double beta_c);
//     double Tsqr_one(cosmo*,double k,double Gamma_eff,error **err);
//     double Tsqr(cosmo*,double k,error **err);*/


    /* Linear power spectrum */
    double P_L(double a, double k){
        return ::P_L(this,a,k,err);        
    }

    /* Smith et al. non-linear power spectrum (halofit) */
    double P(double a, double k){
        return ::P_NL(this,a,k,err);        
    }
    
private:
    error** err;
    error *myerr;
};

#endif // COSMO_HPP
