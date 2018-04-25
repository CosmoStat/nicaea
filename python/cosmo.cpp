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

#include "cosmo.hpp"
#include <iostream>


void cosmology::init(double OMEGAM, double OMEGADE, double W0_DE, double W1_DE, double* W_POLY_DE, int N_POLY_DE, double H100, double OMEGAB, double OMEGANUMASS, double NEFFNUMASS, double NORM, double NSPEC, nonlinear_t NONLINEAR, transfer_t TRANSFER, growth_t GROWTH, de_param_t DEPARAM, norm_t normmode, double AMIN)
{
    err = &myerr;
    
    int N_a_min;
    Omega_m       = OMEGAM;
    Omega_de      = OMEGADE;
    
    if (DEPARAM != poly_DE) {
        w0_de         = W0_DE;
        w1_de         = W1_DE;
        w_poly_de     = NULL;
        N_poly_de     = 0;
    } else {
        set_w_poly_de(&w_poly_de, &N_poly_de, W_POLY_DE, N_POLY_DE, 1, err);
        
        /* MDKDEBUG: know what we are doing here! */
        w0_de = w_poly_de[0];
        w1_de = -1.0e30;
    }
    
    h_100         = H100;
    Omega_b       = OMEGAB;
    Omega_nu_mass = OMEGANUMASS;
    Neff_nu_mass  = NEFFNUMASS;
    normalization = NORM;
    n_spec        = NSPEC;
    
    nonlinear     = NONLINEAR;
    transfer      = TRANSFER;
    growth        = GROWTH;
    de_param      = DEPARAM;
    normmode      = normmode;
    
    if (AMIN>0.0) a_min = AMIN;
    else a_min = a_minmin;
    
    /* Reset pre-computed numbers and tables */
    
    /* New 06/2014: Set N_a such that da >= 0.001 */
    N_a_min = (int)((1.0 - a_min) / 0.001) + 1;
    if ( N_a_min < _N_a) {
        N_a = _N_a;
    } else {
        N_a = N_a_min;
    }
    
    //printf("N_a, da = %d %g\n", N_a, (1.0 - a_min) / (N_a - 1)); //, N_a_min);
    
    linearGrowth         = NULL;
    growth_delta0        = -1;
    transferFct          = NULL;
    transfer_alpha_Gamma = -1;
    transfer_s           = -1;
    transferBE           = NULL;
    cmp_sigma8           = -1;
    P_NL                 = NULL;
    slope                = NULL;
    w                    = NULL;
    //k_NL                 = NULL;
    ystar_allz           = NULL;
    
    //dump_param(res, stdout);
    
    set_norm(this, err);   
    // MKDEBUG NEW
    consistency_parameters(this, err);    
}

