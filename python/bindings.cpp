// This file has been generated by Py++.

#include "boost/python.hpp"

#include "cosmo.hpp"

#include "../Cosmo/include/cosmo.h"

namespace bp = boost::python;

struct cosmo_wrapper : nicaea::cosmo, bp::wrapper< nicaea::cosmo > {

    cosmo_wrapper(nicaea::cosmo const & arg )
    : nicaea::cosmo( arg )
      , bp::wrapper< nicaea::cosmo >(){
        // copy constructor
        
    }

    cosmo_wrapper()
    : nicaea::cosmo()
      , bp::wrapper< nicaea::cosmo >(){
        // null constructor
        
    }

    static ::nicaea::interTable2D * get_P_NL(nicaea::cosmo const & inst ){
        return inst.P_NL;
    }
    
    static void set_P_NL( nicaea::cosmo & inst, ::nicaea::interTable2D * new_value ){ 
        inst.P_NL = new_value;
    }

    static ::nicaea::interTable * get_linearGrowth(nicaea::cosmo const & inst ){
        return inst.linearGrowth;
    }
    
    static void set_linearGrowth( nicaea::cosmo & inst, ::nicaea::interTable * new_value ){ 
        inst.linearGrowth = new_value;
    }

    static ::nicaea::interTable * get_slope(nicaea::cosmo const & inst ){
        return inst.slope;
    }
    
    static void set_slope( nicaea::cosmo & inst, ::nicaea::interTable * new_value ){ 
        inst.slope = new_value;
    }

    static ::nicaea::interTable * get_transferBE(nicaea::cosmo const & inst ){
        return inst.transferBE;
    }
    
    static void set_transferBE( nicaea::cosmo & inst, ::nicaea::interTable * new_value ){ 
        inst.transferBE = new_value;
    }

    static ::nicaea::interTable * get_transferFct(nicaea::cosmo const & inst ){
        return inst.transferFct;
    }
    
    static void set_transferFct( nicaea::cosmo & inst, ::nicaea::interTable * new_value ){ 
        inst.transferFct = new_value;
    }

    static ::nicaea::interTable * get_w(nicaea::cosmo const & inst ){
        return inst.w;
    }
    
    static void set_w( nicaea::cosmo & inst, ::nicaea::interTable * new_value ){ 
        inst.w = new_value;
    }

};

BOOST_PYTHON_MODULE(nicaea){
    bp::class_< cosmo_wrapper >( "cosmo" )    
        .def_readwrite( "As", &nicaea::cosmo::As )    
        .def_readwrite( "N_a", &nicaea::cosmo::N_a )    
        .def_readwrite( "N_poly_de", &nicaea::cosmo::N_poly_de )    
        .def_readwrite( "Neff_nu_mass", &nicaea::cosmo::Neff_nu_mass )    
        .def_readwrite( "Omega_b", &nicaea::cosmo::Omega_b )    
        .def_readwrite( "Omega_de", &nicaea::cosmo::Omega_de )    
        .def_readwrite( "Omega_m", &nicaea::cosmo::Omega_m )    
        .def_readwrite( "Omega_nu_mass", &nicaea::cosmo::Omega_nu_mass )    
        .add_property( "P_NL"
                    , bp::make_function( (::nicaea::interTable2D * (*)( ::nicaea::cosmo const & ))(&cosmo_wrapper::get_P_NL), bp::return_internal_reference< >() )
                    , bp::make_function( (void (*)( ::nicaea::cosmo &,::nicaea::interTable2D * ))(&cosmo_wrapper::set_P_NL), bp::with_custodian_and_ward_postcall< 1, 2 >() ) )    
        .def_readwrite( "a_min", &nicaea::cosmo::a_min )    
        .def_readwrite( "cmp_sigma8", &nicaea::cosmo::cmp_sigma8 )    
        .def_readwrite( "de_param", &nicaea::cosmo::de_param )    
        .def_readwrite( "growth", &nicaea::cosmo::growth )    
        .def_readwrite( "growth_delta0", &nicaea::cosmo::growth_delta0 )    
        .def_readwrite( "h_100", &nicaea::cosmo::h_100 )    
        .add_property( "linearGrowth"
                    , bp::make_function( (::nicaea::interTable * (*)( ::nicaea::cosmo const & ))(&cosmo_wrapper::get_linearGrowth), bp::return_internal_reference< >() )
                    , bp::make_function( (void (*)( ::nicaea::cosmo &,::nicaea::interTable * ))(&cosmo_wrapper::set_linearGrowth), bp::with_custodian_and_ward_postcall< 1, 2 >() ) )    
        .def_readwrite( "n_spec", &nicaea::cosmo::n_spec )    
        .def_readwrite( "nonlinear", &nicaea::cosmo::nonlinear )    
        .def_readwrite( "normalization", &nicaea::cosmo::normalization )    
        .def_readwrite( "normmode", &nicaea::cosmo::normmode )    
        .def_readwrite( "sigma_8", &nicaea::cosmo::sigma_8 )    
        .add_property( "slope"
                    , bp::make_function( (::nicaea::interTable * (*)( ::nicaea::cosmo const & ))(&cosmo_wrapper::get_slope), bp::return_internal_reference< >() )
                    , bp::make_function( (void (*)( ::nicaea::cosmo &,::nicaea::interTable * ))(&cosmo_wrapper::set_slope), bp::with_custodian_and_ward_postcall< 1, 2 >() ) )    
        .def_readwrite( "transfer", &nicaea::cosmo::transfer )    
        .add_property( "transferBE"
                    , bp::make_function( (::nicaea::interTable * (*)( ::nicaea::cosmo const & ))(&cosmo_wrapper::get_transferBE), bp::return_internal_reference< >() )
                    , bp::make_function( (void (*)( ::nicaea::cosmo &,::nicaea::interTable * ))(&cosmo_wrapper::set_transferBE), bp::with_custodian_and_ward_postcall< 1, 2 >() ) )    
        .add_property( "transferFct"
                    , bp::make_function( (::nicaea::interTable * (*)( ::nicaea::cosmo const & ))(&cosmo_wrapper::get_transferFct), bp::return_internal_reference< >() )
                    , bp::make_function( (void (*)( ::nicaea::cosmo &,::nicaea::interTable * ))(&cosmo_wrapper::set_transferFct), bp::with_custodian_and_ward_postcall< 1, 2 >() ) )    
        .def_readwrite( "transfer_alpha_Gamma", &nicaea::cosmo::transfer_alpha_Gamma )    
        .def_readwrite( "transfer_s", &nicaea::cosmo::transfer_s )    
        .add_property( "w"
                    , bp::make_function( (::nicaea::interTable * (*)( ::nicaea::cosmo const & ))(&cosmo_wrapper::get_w), bp::return_internal_reference< >() )
                    , bp::make_function( (void (*)( ::nicaea::cosmo &,::nicaea::interTable * ))(&cosmo_wrapper::set_w), bp::with_custodian_and_ward_postcall< 1, 2 >() ) )    
        .def_readwrite( "w0_de", &nicaea::cosmo::w0_de )    
        .def_readwrite( "w1_de", &nicaea::cosmo::w1_de );

    bp::class_< cosmology, bp::bases< nicaea::cosmo > >( "cosmology", bp::init< >() )    
        .def( 
            "D_a"
            , (double ( ::cosmology::* )( double ))( &::cosmology::D_a )
            , ( bp::arg("a") ) )    
        .def( 
            "D_a12"
            , (double ( ::cosmology::* )( double,double ))( &::cosmology::D_a12 )
            , ( bp::arg("a1"), bp::arg("a2") ) )    
        .def( 
            "D_lum"
            , (double ( ::cosmology::* )( double ))( &::cosmology::D_lum )
            , ( bp::arg("a") ) )    
        .def( 
            "Esqr"
            , (double ( ::cosmology::* )( double,int ))( &::cosmology::Esqr )
            , ( bp::arg("a"), bp::arg("wOmegar") ) )    
        .def( 
            "Omega_de_a"
            , (double ( ::cosmology::* )( double,double ))( &::cosmology::Omega_de_a )
            , ( bp::arg("a"), bp::arg("Esqrpre") ) )    
        .def( 
            "Omega_m_a"
            , (double ( ::cosmology::* )( double,double ))( &::cosmology::Omega_m_a )
            , ( bp::arg("a"), bp::arg("Esqrpre") ) )    
        .def( 
            "P"
            , (double ( ::cosmology::* )( double,double ))( &::cosmology::P )
            , ( bp::arg("a"), bp::arg("k") ) )    
        .def( 
            "P_L"
            , (double ( ::cosmology::* )( double,double ))( &::cosmology::P_L )
            , ( bp::arg("a"), bp::arg("k") ) )    
        .def( 
            "da_dtau"
            , (double ( ::cosmology::* )( double ))( &::cosmology::da_dtau )
            , ( bp::arg("a") ) )    
        .def( 
            "drdz"
            , (double ( ::cosmology::* )( double ))( &::cosmology::drdz )
            , ( bp::arg("a") ) )    
        .def( 
            "dvdz"
            , (double ( ::cosmology::* )( double ))( &::cosmology::dvdz )
            , ( bp::arg("a") ) )    
        .def( 
            "dwoverda"
            , (double ( ::cosmology::* )( double ))( &::cosmology::dwoverda )
            , ( bp::arg("a") ) )    
        .def( 
            "f_K"
            , (double ( ::cosmology::* )( double ))( &::cosmology::f_K )
            , ( bp::arg("wr") ) )    
        .def( 
            "f_de"
            , (double ( ::cosmology::* )( double ))( &::cosmology::f_de )
            , ( bp::arg("a") ) )    
        .def( 
            "w_a"
            , (double ( ::cosmology::* )( double,int ))( &::cosmology::w_a )
            , ( bp::arg("a"), bp::arg("wOmegar") ) )    
        .def( 
            "w_de"
            , (double ( ::cosmology::* )( double ))( &::cosmology::w_de )
            , ( bp::arg("a") ) )    
        .def( 
            "w_nu_mass"
            , (double ( ::cosmology::* )( double ))( &::cosmology::w_nu_mass )
            , ( bp::arg("a") ) );
}
