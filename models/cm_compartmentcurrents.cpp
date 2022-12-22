/*
 *  cm_compartmentcurrents.cpp
 *
 *  This file is part of NEST.
 *
 *  Copyright (C) 2004 The NEST Initiative
 *
 *  NEST is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  NEST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#include "cm_compartmentcurrents.h"
bool debug = false;
bool Larkum_Ca = false;    // Larkum vs Hay/Branco Ca current dynamics
bool Hay_K_Ca = false;     // Hay vs Branco K_Ca current dynamics

nest::Na::Na()
  // initialization state state variables
  : m_Na_( 0.0 )
  , h_Na_( 0.0 )
  // initialization parameters
  , gbar_Na_( 0.0 )
  , e_Na_( 50. )
{
}
nest::Na::Na( const DictionaryDatum& channel_params )
  // initialization state state variables
  : m_Na_( 0.0 )
  , h_Na_( 0.0 )
  // initialization parameters
  , gbar_Na_( 0.0 )
  , e_Na_( 50. )
{
  // update sodium channel parameters
  if ( channel_params->known( "gbar_Na" ) )
  {
    gbar_Na_ = getValue< double >( channel_params, "gbar_Na" );
  }
  if ( channel_params->known( "e_Na" ) )
  {
    e_Na_ = getValue< double >( channel_params, "e_Na" );
  }
}

void
nest::Na::append_recordables( std::map< Name, double* >* recordables, const long compartment_idx )
{
  ( *recordables )[ Name( "m_Na_" + std::to_string( compartment_idx ) ) ] = &m_Na_;
  ( *recordables )[ Name( "h_Na_" + std::to_string( compartment_idx ) ) ] = &h_Na_;
}

std::pair< double, double >
nest::Na::f_numstep( const double v_comp )
{
  const double dt = Time::get_resolution().get_ms();
  double g_val = 0., i_val = 0.;

  if ( gbar_Na_ > 1e-9 )
  {
    /**
     * Channel rate equations from the following .mod file:
     * https://senselab.med.yale.edu/ModelDB/ShowModel?model=140828&file=/Branco_2010/mod.files/na.mod#tabs-2
     */
    // auxiliary variables
    double v_comp_plus_35 = v_comp + 35.013;

    // trap the case where alpha_m and beta_m are 0/0 by substituting explicitly
    // precomputed limiting values
    double alpha_m, frac_alpha_plus_beta_m;
    if ( std::abs( v_comp_plus_35 ) > 1e-5 )
    {
      double exp_vcp35_div_9 = std::exp( 0.111111111111111 * v_comp_plus_35 );
      double frac_evcp35d9 = 1. / ( exp_vcp35_div_9 - 1. );

      alpha_m = 0.182 * v_comp_plus_35 * exp_vcp35_div_9 * frac_evcp35d9;
      double beta_m = 0.124 * v_comp_plus_35 * frac_evcp35d9;
      frac_alpha_plus_beta_m = 1. / ( alpha_m + beta_m );
    }
    else
    {
      alpha_m = 1.638;
      frac_alpha_plus_beta_m = 1. / ( alpha_m + 1.116 );
    }

    double v_comp_plus_50 = v_comp + 50.013;
    double v_comp_plus_75 = v_comp + 75.013;

    // trap the case where alpha_h or beta_h are 0/0 by substituting
    // precomputed limiting values
    double alpha_h, beta_h;
    if ( std::abs( v_comp_plus_50 ) > 1e-5 )
    {
      alpha_h = 0.024 * v_comp_plus_50 / ( 1.0 - std::exp( -0.2 * v_comp_plus_50 ) );
    }
    else
    {
      alpha_h = 0.12;
    }
    if ( std::abs( v_comp_plus_75 ) > 1e-9 )
    {
      beta_h = -0.0091 * v_comp_plus_75 / ( 1.0 - std::exp( 0.2 * v_comp_plus_75 ) );
    }
    else
    {
      beta_h = 0.0455;
    }

    // activation and timescale for state variable 'm'
    double tau_m_Na = q10_ * frac_alpha_plus_beta_m;
    double m_inf_Na = alpha_m * frac_alpha_plus_beta_m;

    // activation and timescale for state variable 'h'
    double tau_h_Na = q10_ / ( alpha_h + beta_h );
    double h_inf_Na = 1. / ( 1. + std::exp( ( v_comp + 65. ) / 6.2 ) );

    // advance state variable 'm' one timestep
    double p_m_Na = std::exp( -dt / tau_m_Na );
    m_Na_ *= p_m_Na;
    m_Na_ += ( 1. - p_m_Na ) * m_inf_Na;

    // advance state variable 'h' one timestep
    double p_h_Na = std::exp( -dt / tau_h_Na );
    h_Na_ *= p_h_Na;
    h_Na_ += ( 1. - p_h_Na ) * h_inf_Na;

    // compute the conductance of the sodium channel
    double g_Na = gbar_Na_ * std::pow( m_Na_, 3 ) * h_Na_;

    // add to variables for numerical integration
    g_val += g_Na / 2.;
    i_val += g_Na * ( e_Na_ - v_comp / 2. );
  }
  if (debug)
    printf("i_Na     = %f\n",i_val);
  return std::make_pair( g_val, i_val );
}


nest::K::K()
  // initialization state variables
  : n_K_( 0.0 )
  // initialization parameters
  , gbar_K_( 0.0 )
  , e_K_( -85. )
{
}
nest::K::K( const DictionaryDatum& channel_params )
  // initialization state variables
  : n_K_( 0.0 )
  // initialization parameters
  , gbar_K_( 0.0 )
  , e_K_( -85. )
{
  // update potassium channel parameters
  if ( channel_params->known( "gbar_K" ) )
  {
    gbar_K_ = getValue< double >( channel_params, "gbar_K" );
  }
  if ( channel_params->known( "e_Na" ) )
  {
    e_K_ = getValue< double >( channel_params, "e_K" );
  }
}

void
nest::K::append_recordables( std::map< Name, double* >* recordables, const long compartment_idx )
{
  ( *recordables )[ Name( "n_K_" + std::to_string( compartment_idx ) ) ] = &n_K_;
}

std::pair< double, double >
nest::K::f_numstep( const double v_comp )
{
  const double dt = Time::get_resolution().get_ms();
  double g_val = 0., i_val = 0.;

  if ( gbar_K_ > 1e-9 )
  {
    /**
     * Channel rate equations from the following .mod file:
     * https://senselab.med.yale.edu/ModelDB/ShowModel?model=140828&file=/Branco_2010/mod.files/kv.mod#tabs-2
     */
    // auxiliary variables
    double v_comp_minus_25 = v_comp - 25.;

    // trap the case where alpha_n and beta_n are 0/0 by substituting explicitly
    // precomputed limiting values
    double alpha_n, frac_alpha_plus_beta_n;
    if ( std::abs( v_comp_minus_25 ) > 1e-5 )
    {
      double exp_vm25_div_9 = std::exp( 0.111111111111111 * v_comp_minus_25 );
      double frac_evm25d9 = 1. / ( exp_vm25_div_9 - 1. );

      alpha_n = 0.02 * v_comp_minus_25 * exp_vm25_div_9 * frac_evm25d9;
      double beta_n = 0.002 * v_comp_minus_25 * frac_evm25d9;
      frac_alpha_plus_beta_n = 1. / ( alpha_n + beta_n );
    }
    else
    {
      alpha_n = 0.18;
      double beta_n = 0.018;
      frac_alpha_plus_beta_n = 1. / ( alpha_n + beta_n );
    }

    // activation and timescale of state variable 'n'
    double tau_n_K = q10_ * frac_alpha_plus_beta_n;
    double n_inf_K = alpha_n * frac_alpha_plus_beta_n;

    //printf("tau_n_K = %f\n",tau_n_K);

    // advance state variable 'm' one timestep
    double p_n_K = std::exp( -dt / tau_n_K );
    n_K_ *= p_n_K;
    n_K_ += ( 1. - p_n_K ) * n_inf_K;

    // compute the conductance of the potassium channel
    double g_K = gbar_K_ * n_K_;
    //double g_K = gbar_K_ * std::pow(n_K_,4);

    // add to variables for numerical integration
    g_val += g_K / 2.;
    i_val += g_K * ( e_K_ - v_comp / 2. );
  }
  if (debug)
    printf("i_K      = %f\n",i_val);
  return std::make_pair( g_val, i_val );
}


nest::Ca::Ca()
  // initialization state variables
  : m_Ca_( 0.0 )
  , h_Ca_( 0.0 )
  // initialization parameters
  , gbar_Ca_( 0.0 )
  , m_half_(0.0)
  , h_half_(0.0)
  , m_slope_(0.0)
  , h_slope_(0.0)
  , tau_m_(0.0)
  , tau_h_(0.0)
  , tau_decay_Ca_( 120 )
  , scale_( 1. )
  , v_rest_( 0.0 )
    //, g_AHP_( 0.0 )
    //, e_K_AHP_( -90 )
    //, Ca_conc_(0.0)
    //, tau_K_Ca_(80.0)
    //, scale_( 1. )
{
}
nest::Ca::Ca( const DictionaryDatum& channel_params )
  // initialization state state variables
  : m_Ca_( 0.0 )
  , h_Ca_( 0.0 )
  // initialization parameters
  , gbar_Ca_( 0.0 )
  , m_half_(0.0)
  , h_half_(0.0)
  , m_slope_(0.0)
  , h_slope_(0.0)
  , tau_m_(0.0)
  , tau_h_(0.0)
  , tau_decay_Ca_( 120 )
  , scale_( 1. )
  , v_rest_( 0.0 )
    //, g_AHP_(0.0 )
    //, e_K_AHP_( -90 )
    //, Ca_conc_(0.0)
    //, tau_K_Ca_(80.0)
    //, scale_( 1. )
{
  // update AHP channel parameters
  if ( channel_params->known( "gbar_Ca" ) )
  {
    gbar_Ca_ = getValue< double >( channel_params, "gbar_Ca" );
  }
  if ( channel_params->known( "m_half" ) )
  {
    m_half_ = getValue< double >( channel_params, "m_half" );
  }
  if ( channel_params->known( "h_half" ) )
  {
    h_half_ = getValue< double >( channel_params, "h_half" );
  }
  if ( channel_params->known( "m_slope" ) )
  {
    m_slope_ = getValue< double >( channel_params, "m_slope" );
  }
  if ( channel_params->known( "h_slope" ) )
  {
    h_slope_ = getValue< double >( channel_params, "h_slope" );
  }
  if ( channel_params->known( "tau_m" ) )
  {
    tau_m_ = getValue< double >( channel_params, "tau_m" );
  }
  if ( channel_params->known( "tau_h" ) )
  {
    tau_h_ = getValue< double >( channel_params, "tau_h" );
  }
  if ( channel_params->known( "tau_decay_Ca" ) )
  {
    tau_decay_Ca_ = getValue< double >( channel_params, "tau_decay_Ca" );
  }
  if ( channel_params->known( "scale" ) )
  {
    scale_ = getValue< double >( channel_params, "scale" );
  }
  if ( channel_params->known( "e_L" ) )
  {
    v_rest_ = getValue< double >( channel_params, "e_L" );
  }

  if(Larkum_Ca)
    {
      // Dynamics from Larkum et al. 2004 and Chua et al.2015
      // State variable 'm' initialization
      m_Ca_ = 1. / ( 1. + std::exp( m_slope_ * ( v_rest_ - m_half_ ) ) );
      // State variable 'h' initialization
      h_Ca_ = 1. / ( 1. + std::exp( h_slope_ * ( v_rest_ - h_half_ ) ) );
    }
  else
    {
      // Dynamics from https://senselab.med.yale.edu/ModelDB/ShowModel?model=140828&file=/Branco_2010/mod.files/ca.mod#tabs-2 and from Hay et al. 2011 (I_Ca_HVA)
      // State variable 'm' initialization
      double alpha_m_ = -0.055 * (v_rest_ + 27) / (std::exp(-(v_rest_+27)/3.8) - 1);
      double beta_m_ = 0.94 * std::exp(-(v_rest_+75)/17);
      m_Ca_ = alpha_m_ / (alpha_m_ + beta_m_);
      
      // State variable 'h' initialization
      double alpha_h_ =  0.000457 * std::exp(-(v_rest_+13)/50);
      double beta_h_ = 0.0065 / (std::exp(-(v_rest_+15)/28) + 1);
      h_Ca_ = alpha_h_ / (alpha_h_ + beta_h_);
    }
  
}

void
nest::Ca::append_recordables( std::map< Name, double* >* recordables, const long compartment_idx )
{
  ( *recordables )[ Name( "m_Ca_" + std::to_string( compartment_idx ) ) ] = &m_Ca_;
  ( *recordables )[ Name( "h_Ca_" + std::to_string( compartment_idx ) ) ] = &h_Ca_;
  ( *recordables )[ Name( "e_Ca_" + std::to_string( compartment_idx ) ) ] = &e_Ca_;
}

std::pair< double, double >
nest::Ca::f_numstep( const double v_comp )
{
  const double dt = Time::get_resolution().get_ms();
  double g_val = 0., i_val = 0.;

  double k = 1000;
  double R = 8.31441;
  double T = 309.15;
  double F = 96489;
  double Ca_o = 2.;
  double Ca_0_ = 0.0001;    // mM (from Hay et al. 2011 PlosCompBio and Gerstner 2014 book)

  if ( gbar_Ca_ > 1e-9 )
  {
    double tau_m_Ca;
    double m_inf_Ca;
    double tau_h_Ca;
    double h_inf_Ca;
    double g_Ca;
    double I_Ca_;

    if(Larkum_Ca)
      {      
	// Dynamics from Larkum et al. 2004 and Chua et al.2015
	// activation and timescale for state variable 'm' 
        tau_m_Ca = tau_m_;
	m_inf_Ca = 1. / ( 1. + std::exp( m_slope_ * ( v_comp - m_half_ ) ) );
	
	// activation and timescale for state variable 'h'
	tau_h_Ca = tau_h_;
	h_inf_Ca = 1. / ( 1. + std::exp( h_slope_ * ( v_comp - h_half_ ) ) );
      }
    else
      {
	// Dynamics from https://senselab.med.yale.edu/ModelDB/ShowModel?model=140828&file=/Branco_2010/mod.files/ca.mod#tabs-2 and from Hay et al. 2011 (I_Ca_HVA)
	// activation and timescale for state variable 'm'
	double alpha_m_ = -0.055 * (v_comp + 27) / (std::exp(-(v_comp+27)/3.8) - 1);
	double beta_m_ = 0.94 * std::exp(-(v_comp+75)/17);
    
	tau_m_Ca = 1 / (alpha_m_ + beta_m_);
	m_inf_Ca = alpha_m_ / (alpha_m_ + beta_m_);

	// activation and timescale for state variable 'h'
	double alpha_h_ =  0.000457 * std::exp(-(v_comp+13)/50);
	double beta_h_ = 0.0065 / (std::exp(-(v_comp+15)/28) + 1);

	tau_h_Ca = 1 / (alpha_h_ + beta_h_);
	h_inf_Ca = alpha_h_ / (alpha_h_ + beta_h_);
      }
       
    /*
    // activation and timescale for state variable 'm' - from Destexhe 1994
    double tau_m_Ca = 0.44 + 0.15 / (std::exp( (v_comp + 27.) / 10. ) + std::exp( -(v_comp + 102.) / 15. ));
    double m_inf_Ca = 1. / ( 1. + std::exp( -(v_comp + 52.) / 7.4 ));
    
    // activation and timescale for state variable 'h' - from Destexhe 1994
    double tau_h_Ca = 22.7 + 0.27 / (std::exp( (v_comp + 48.) / 4. ) + std::exp( -(v_comp + 407.) / 50. ));
    double h_inf_Ca = 1. / ( 1. + std::exp( (v_comp + 80.) / 5. )); 
    */

    // advance state variable 'm' one timestep
    double p_m_Ca = std::exp( -dt / tau_m_Ca );
    m_Ca_ *= p_m_Ca;
    m_Ca_ += ( 1. - p_m_Ca ) * m_inf_Ca;

    // advance state variable 'h' one timestep
    double p_h_Ca = std::exp( -dt / tau_h_Ca );
    h_Ca_ *= p_h_Ca;
    h_Ca_ += ( 1. - p_h_Ca ) * h_inf_Ca;

    // compute the conductance of the calcium channel
    if(Larkum_Ca)
      g_Ca = gbar_Ca_ * m_Ca_ * h_Ca_;
    else
      g_Ca = gbar_Ca_ * std::pow( m_Ca_, 2 ) * h_Ca_;

    // Calcium reversal potential
    e_Ca_ = k * R * T * log(Ca_o/Ca_conc_)/(2*F);
    //e_Ca_ = 30.;

    // Calcium current
    if(Larkum_Ca)
      I_Ca_ = gbar_Ca_ * m_Ca_ * h_Ca_ * (e_Ca_ - v_comp);
    else
      I_Ca_ = gbar_Ca_ * std::pow( m_Ca_, 2 ) * h_Ca_ * (e_Ca_ - v_comp);

    // Calcium concentration
    // Dynamics from https://senselab.med.yale.edu/ModelDB/ShowModel?model=140828&file=/Branco_2010/mod.files/cad.mod#tabs-2 and from Hay et al. 2011
    double Ca_Pump_ = (Ca_conc_ - Ca_0_) * std::exp( -dt / tau_decay_Ca_) + Ca_0_;
    double Ca_Influx_ = scale_ * I_Ca_ * dt;
    if (Ca_Influx_ <= 0.0)
      Ca_Influx_ = 0.0;
    Ca_conc_ = Ca_Pump_ + Ca_Influx_;
    
    /*
    // Calcium concentration
    // Dynamics from https://senselab.med.yale.edu/ModelDB/ShowModel?model=140828&file=/Branco_2010/mod.files/cad.mod#tabs-2 and from Hay et al. 2011
    // Implementation with propagators as in other NEST code
    double p_Ca_conc_ = std::exp( -dt / tau_decay_Ca_ );
    Ca_conc_ *= p_Ca_conc_;
    Ca_conc_ += ( 1. - p_Ca_conc_ ) * Ca_0_;
    double Ca_Influx_ = scale_ * I_Ca_ * dt;
    if (Ca_Influx_ <= 0.0)
      Ca_Influx_ = 0.0;
    Ca_conc_ += Ca_Influx_;
    */

    // add to variables for numerical integration
    g_val += g_Ca / 2.;
    i_val += g_Ca * ( e_Ca_ - v_comp / 2. );
    if (debug)
      printf("i_Ca     = %f\n",i_val);

    if (debug){
      printf("i_Ca = %f     -----    Ca_conc = %.10e \n",i_val,Ca_conc_);
      printf("E_Ca = %f\n",e_Ca_);
      printf("m_Ca = %f     -----    h_Ca = %.10e \n",m_Ca_,h_Ca_);
    }

    /*
    // Compute AHP current
    if ( g_AHP_ > 1e-9 ){
      Ca_conc_ *= std::exp( -dt / tau_K_Ca_);
      Ca_conc_ += scale_ * m_Ca_ * h_Ca_ * dt;
      g_val += Ca_conc_ * g_AHP_ / 2.;
      i_val += Ca_conc_ * g_AHP_ * ( e_K_AHP_ - v_comp / 2. );
      if (debug){
	double temp = Ca_conc_ * g_AHP_ * ( e_K_AHP_ - v_comp / 2. );
	printf("i_Ca_AHP = %f\n",i_val);
      }
    }
    */

  }

  return std::make_pair( g_val, i_val );
}



nest::K_Ca::K_Ca()
  // initialization state variables
  : m_K_Ca_( 0.0 )
  , m_Ca_( 0.0 )
  , h_Ca_( 0.0 )
  , Ca_conc_( 0.0001 )
  // initialization parameters
  , gbar_K_Ca_( 0.0 )
  , gbar_Ca_( 0.0 )
  , e_K_( -90 )
  , m_half_(0.0)
  , h_half_(0.0)
  , m_slope_(0.0)
  , h_slope_(0.0)
  , tau_m_(0.0)
  , tau_h_(0.0)
  , tau_decay_Ca_( 120 )
  , scale_( 1. )
  , v_rest_( 0.0 )
{
}
nest::K_Ca::K_Ca( const DictionaryDatum& channel_params )
  // initialization state state variables
  : m_K_Ca_( 0.0 )
  , m_Ca_( 0.0 )
  , h_Ca_( 0.0 )
  , Ca_conc_( 0.0001 )
  // initialization parameters
  , gbar_K_Ca_(0.0 )
  , gbar_Ca_( 0.0 )
  , e_K_( -90 )
  , m_half_(0.0)
  , h_half_(0.0)
  , m_slope_(0.0)
  , h_slope_(0.0)
  , tau_m_(0.0)
  , tau_h_(0.0)
  , tau_decay_Ca_( 120 )
  , scale_( 1. )
  , v_rest_( 0.0 )
{
  // update K_Ca channel parameters
  if ( channel_params->known( "gbar_K_Ca" ) )
  {
    gbar_K_Ca_ = getValue< double >( channel_params, "gbar_K_Ca" );
  }
  if ( channel_params->known( "gbar_Ca" ) )
  {
    gbar_Ca_ = getValue< double >( channel_params, "gbar_Ca" );
  }
  if ( channel_params->known( "e_K" ) )
  {
    e_K_ = getValue< double >( channel_params, "e_K" );
  }
  if ( channel_params->known( "m_half" ) )
  {
    m_half_ = getValue< double >( channel_params, "m_half" );
  }
  if ( channel_params->known( "h_half" ) )
  {
    h_half_ = getValue< double >( channel_params, "h_half" );
  }
  if ( channel_params->known( "m_slope" ) )
  {
    m_slope_ = getValue< double >( channel_params, "m_slope" );
  }
  if ( channel_params->known( "h_slope" ) )
  {
    h_slope_ = getValue< double >( channel_params, "h_slope" );
  }
  if ( channel_params->known( "tau_m" ) )
  {
    tau_m_ = getValue< double >( channel_params, "tau_m" );
  }
  if ( channel_params->known( "tau_h" ) )
  {
    tau_h_ = getValue< double >( channel_params, "tau_h" );
  }
  if ( channel_params->known( "tau_decay_Ca" ) )
  {
    tau_decay_Ca_ = getValue< double >( channel_params, "tau_decay_Ca" );
  }
  if ( channel_params->known( "scale" ) )
  {
    scale_ = getValue< double >( channel_params, "scale" );
  }
  if ( channel_params->known( "e_L" ) )
  {
    v_rest_ = getValue< double >( channel_params, "e_L" );
  }

  if(Larkum_Ca)
    {
      // Dynamics from Larkum et al. 2004 and Chua et al.2015
      // State variable 'm' initialization
      m_Ca_ = 1. / ( 1. + std::exp( m_slope_ * ( v_rest_ - m_half_ ) ) );
      // State variable 'h' initialization
      h_Ca_ = 1. / ( 1. + std::exp( h_slope_ * ( v_rest_ - h_half_ ) ) );
    }
  else
    {
      // Dynamics from https://senselab.med.yale.edu/ModelDB/ShowModel?model=140828&file=/Branco_2010/mod.files/ca.mod#tabs-2 and from Hay et al. 2011 (I_Ca_HVA)
      // State variable 'm' initialization
      double alpha_m_ = -0.055 * (v_rest_ + 27) / (std::exp(-(v_rest_+27)/3.8) - 1);
      double beta_m_ = 0.94 * std::exp(-(v_rest_+75)/17);
      m_Ca_ = alpha_m_ / (alpha_m_ + beta_m_);
      
      // State variable 'h' initialization
      double alpha_h_ =  0.000457 * std::exp(-(v_rest_+13)/50);
      double beta_h_ = 0.0065 / (std::exp(-(v_rest_+15)/28) + 1);
      h_Ca_ = alpha_h_ / (alpha_h_ + beta_h_);
    }

  if(Hay_K_Ca)
    {
      // Dynamics from Hay et al. 2011, current I_SK
      m_K_Ca_ = 1. / ( 1 + std::pow(0.00043/Ca_conc_,4.8));
    }
  else
    {
      // Dynamics from https://senselab.med.yale.edu/ModelDB/ShowModel?model=140828&file=/Branco_2010/mod.files/kca.mod#tabs-2
      double alpha_K_Ca = 0.01 * Ca_conc_;
      double beta_K_Ca = 0.02;
      m_K_Ca_ = alpha_K_Ca / (alpha_K_Ca + beta_K_Ca);
      }

}

void
nest::K_Ca::append_recordables( std::map< Name, double* >* recordables, const long compartment_idx )
{
  ( *recordables )[ Name( "m_K_Ca_" + std::to_string( compartment_idx ) ) ] = &m_K_Ca_;
  ( *recordables )[ Name( "Ca_conc_" + std::to_string( compartment_idx ) ) ] = &Ca_conc_;
  ( *recordables )[ Name( "e_Ca_K_Ca_" + std::to_string( compartment_idx ) ) ] = &e_Ca_;
}

std::pair< double, double >
nest::K_Ca::f_numstep( const double v_comp )
{
  const double dt = Time::get_resolution().get_ms();
  double g_val = 0., i_val = 0.;
  double Ca_0_ = 0.0001;    // mM (from Hay et al. 2011 PlosCompBio and Gerstner 2014 book)
  double Ca_th_ = 0.00043;  // Threshold for Ca channel opening

  double k = 1000;
  double R = 8.31441;
  double T = 309.15;
  double F = 96489;
  double d = 1;
  double Ca_o = 2.;
  double k_phi = 0.1*1000000;
  //double phi_Ca_ = k_phi / (2 * F * d);

  if ( gbar_K_Ca_ > 1e-9 )
  {
    double tau_m_Ca;
    double m_inf_Ca;
    double tau_h_Ca;
    double h_inf_Ca;
    double g_Ca;
    double I_Ca_;

    if(Larkum_Ca)
      {
	// Dynamics from Larkum et al. 2004 and Chua et al.2015
	// activation and timescale for state variable 'm' 
	tau_m_Ca = tau_m_;
	m_inf_Ca = 1. / ( 1. + std::exp( m_slope_ * ( v_comp - m_half_ ) ) );
	
	// activation and timescale for state variable 'h'
	tau_h_Ca = tau_h_;
	h_inf_Ca = 1. / ( 1. + std::exp( h_slope_ * ( v_comp - h_half_ ) ) );
      }
    else
      {
	// Dynamics from https://senselab.med.yale.edu/ModelDB/ShowModel?model=140828&file=/Branco_2010/mod.files/ca.mod#tabs-2 and from Hay et al. 2011 (I_Ca_HVA)
	// activation and timescale for state variable 'm'
	double alpha_m_ = -0.055 * (v_comp + 27) / (std::exp(-(v_comp+27)/3.8) - 1);
	double beta_m_ = 0.94 * std::exp(-(v_comp+75)/17);
    
	tau_m_Ca = 1 / (alpha_m_ + beta_m_);
	m_inf_Ca = alpha_m_ / (alpha_m_ + beta_m_);

	// activation and timescale for state variable 'h'
	double alpha_h_ =  0.000457 * std::exp(-(v_comp+13)/50);
	double beta_h_ = 0.0065 / (std::exp(-(v_comp+15)/28) + 1);

	tau_h_Ca = 1 / (alpha_h_ + beta_h_);
	h_inf_Ca = alpha_h_ / (alpha_h_ + beta_h_);
      }
       
    /*
    // activation and timescale for state variable 'm' - from Destexhe 1994
    double tau_m_Ca = 0.44 + 0.15 / (std::exp( (v_comp + 27.) / 10. ) + std::exp( -(v_comp + 102.) / 15. ));
    double m_inf_Ca = 1. / ( 1. + std::exp( -(v_comp + 52.) / 7.4 ));
    
    // activation and timescale for state variable 'h' - from Destexhe 1994
    double tau_h_Ca = 22.7 + 0.27 / (std::exp( (v_comp + 48.) / 4. ) + std::exp( -(v_comp + 407.) / 50. ));
    double h_inf_Ca = 1. / ( 1. + std::exp( (v_comp + 80.) / 5. )); 
    */

    // advance state variable 'm' one timestep
    double p_m_Ca = std::exp( -dt / tau_m_Ca );
    m_Ca_ *= p_m_Ca;
    m_Ca_ += ( 1. - p_m_Ca ) * m_inf_Ca;

    // advance state variable 'h' one timestep
    double p_h_Ca = std::exp( -dt / tau_h_Ca );
    h_Ca_ *= p_h_Ca;
    h_Ca_ += ( 1. - p_h_Ca ) * h_inf_Ca;

    // Calcium reversal potential
    e_Ca_ = k * R * T * log(Ca_o/Ca_conc_)/(2*F);
    //e_Ca_ = 30.;

    // Calcium current
    if(Larkum_Ca)
      I_Ca_ = gbar_Ca_ * m_Ca_ * h_Ca_ * (e_Ca_ - v_comp);
    else
      I_Ca_ = gbar_Ca_ * std::pow( m_Ca_, 2 ) * h_Ca_ * (e_Ca_ - v_comp);
    
    // Calcium concentration
    // Dynamics from https://senselab.med.yale.edu/ModelDB/ShowModel?model=140828&file=/Branco_2010/mod.files/cad.mod#tabs-2 and from Hay et al. 2011
    double Ca_Pump_ = (Ca_conc_ - Ca_0_) * std::exp( -dt / tau_decay_Ca_) + Ca_0_;
    double Ca_Influx_ = scale_ * I_Ca_ * dt;
    if (Ca_Influx_ <= 0.0)
      Ca_Influx_ = 0.0;
    Ca_conc_ = Ca_Pump_ + Ca_Influx_;
    
    /*
    // Calcium concentration
    // Dynamics from https://senselab.med.yale.edu/ModelDB/ShowModel?model=140828&file=/Branco_2010/mod.files/cad.mod#tabs-2 and from Hay et al. 2011
    // Implementation with propagators as in other NEST code
    double p_Ca_conc_ = std::exp( -dt / tau_decay_Ca_ );
    Ca_conc_ *= p_Ca_conc_;
    Ca_conc_ += ( 1. - p_Ca_conc_ ) * Ca_0_;
    double Ca_Influx_ = scale_ * I_Ca_ * dt;
    if (Ca_Influx_ <= 0.0)
      Ca_Influx_ = 0.0;
    Ca_conc_ += Ca_Influx_;
    */

    // Calculate state variable 'm_K_Ca_'
    if(Hay_K_Ca)
      {
	// Dynamics from Hay et al. 2011, current I_SK
	double m_inf_K_Ca_ = 1. / ( 1 + std::pow(Ca_th_/Ca_conc_,4.8));
	double tau_m_K_Ca_ = 1.;
	double p_m_K_Ca = std::exp( -dt / tau_m_K_Ca_ );
	m_K_Ca_ *= p_m_K_Ca;
	m_K_Ca_ += ( 1. - p_m_K_Ca ) * m_inf_K_Ca_;
      }
    else
      {
	// Dynamics from https://senselab.med.yale.edu/ModelDB/ShowModel?model=140828&file=/Branco_2010/mod.files/kca.mod#tabs-2
	double alpha_K_Ca = 0.01 * Ca_conc_;
	double beta_K_Ca = 0.02;
	double m_inf_K_Ca_ = alpha_K_Ca / (alpha_K_Ca + beta_K_Ca);
	double tau_m_K_Ca_ = 1 / (alpha_K_Ca + beta_K_Ca);
	double p_m_K_Ca = std::exp( -dt / tau_m_K_Ca_ );
	m_K_Ca_ *= p_m_K_Ca;
	m_K_Ca_ += ( 1. - p_m_K_Ca ) * m_inf_K_Ca_;
      }

    // compute the conductance of the K_Ca channel
    double g_K_Ca_ = gbar_K_Ca_ * m_K_Ca_;

    // add to variables for numerical integration
    g_val += g_K_Ca_ / 2.;
    i_val += g_K_Ca_ * ( e_K_ - v_comp / 2. );

    if (debug){
      printf("i_K_Ca = %f     -----    Ca_conc = %.10e \n",i_val,Ca_conc_);
      printf("E_Ca = %f\n",e_Ca_);
      printf("m_Ca = %f     -----    h_Ca = %.10e \n",m_Ca_,h_Ca_);
    }
    
  }

  return std::make_pair( g_val, i_val );
}


nest::Na_Adex::Na_Adex()
  // initialization state state variables
  : m_Na_Adex_( 0.0 )
  , h_Na_Adex_( 0.0 )
  // initialization parameters
  , gbar_Na_Adex_( 0.0 )
  , e_Na_Adex_( 0.0 )
  , delta_T_( 0.0 )
{
}
nest::Na_Adex::Na_Adex( const DictionaryDatum& channel_params )
  // initialization state state variables
  : m_Na_Adex_( 0.0 )
  , h_Na_Adex_( 0.0 )
  // initialization parameters
  , gbar_Na_Adex_( 0.0 )
  , e_Na_Adex_( 0.0 )
  , delta_T_( 0.0 )
{
  // update sodium channel parameters
  if ( channel_params->known( "gbar_Na_Adex" ) )
  {
    gbar_Na_Adex_ = getValue< double >( channel_params, "gbar_Na_Adex" );
  }
  if ( channel_params->known( "e_Na_Adex" ) )
  {
    e_Na_Adex_ = getValue< double >( channel_params, "e_Na_Adex" );
  }
  if ( channel_params->known( "delta_T" ) )
  {
    delta_T_ = getValue< double >( channel_params, "delta_T" );
  }
}

void
nest::Na_Adex::append_recordables( std::map< Name, double* >* recordables, const long compartment_idx )
{
  ( *recordables )[ Name( "m_Na_Adex_" + std::to_string( compartment_idx ) ) ] = &m_Na_Adex_;
  ( *recordables )[ Name( "h_Na_Adex_" + std::to_string( compartment_idx ) ) ] = &h_Na_Adex_;
}

std::pair< double, double >
nest::Na_Adex::f_numstep( const double v_comp )
{
  double g_val = 0., i_val = 0.;

  double i_exp_;
  double E_Na_ = 55.;

  if ( gbar_Na_Adex_ > 1e-9 )
  {
    //g_val += gbar_Na_Adex_ / 2.;
    i_exp_ = gbar_Na_Adex_ * delta_T_ * std::exp((v_comp - e_Na_Adex_) / delta_T_);

    g_val += i_exp_ / (2 * fabs(E_Na_ - v_comp)) ;
    i_val += i_exp_ / fabs(E_Na_ - v_comp) * (E_Na_ - v_comp/2);

  }
  
  if (debug)
    printf("i_Na_Adex     = %f\n",i_val);
  return std::make_pair( g_val, i_val );
}




nest::AMPA::AMPA( const long syn_index )
  // initialization state variables
  : g_r_AMPA_( 0.0 )
  , g_d_AMPA_( 0.0 )
  // initialization parameters
  , e_rev_( 0.0 )
  , tau_r_( 0.2 )
  , tau_d_( 3.0 )
{
  syn_idx = syn_index;

  double tp = ( tau_r_ * tau_d_ ) / ( tau_d_ - tau_r_ ) * std::log( tau_d_ / tau_r_ );
  g_norm_ = 1. / ( -std::exp( -tp / tau_r_ ) + std::exp( -tp / tau_d_ ) );
}
nest::AMPA::AMPA( const long syn_index, const DictionaryDatum& receptor_params )
  // initialization state variables
  : g_r_AMPA_( 0.0 )
  , g_d_AMPA_( 0.0 )
  // initialization parameters
  , e_rev_( 0.0 )
  , tau_r_( 0.2 )
  , tau_d_( 3.0 )
{
  syn_idx = syn_index;

  // update AMPA receptor parameters
  if ( receptor_params->known( "e_AMPA" ) )
  {
    e_rev_ = getValue< double >( receptor_params, "e_AMPA" );
  }
  if ( receptor_params->known( "tau_r_AMPA" ) )
  {
    tau_r_ = getValue< double >( receptor_params, "tau_r_AMPA" );
  }
  if ( receptor_params->known( "tau_d_AMPA" ) )
  {
    tau_d_ = getValue< double >( receptor_params, "tau_d_AMPA" );
  }

  double tp = ( tau_r_ * tau_d_ ) / ( tau_d_ - tau_r_ ) * std::log( tau_d_ / tau_r_ );
  g_norm_ = 1. / ( -std::exp( -tp / tau_r_ ) + std::exp( -tp / tau_d_ ) );
}

void
nest::AMPA::append_recordables( std::map< Name, double* >* recordables )
{
  ( *recordables )[ Name( "g_r_AMPA_" + std::to_string( syn_idx ) ) ] = &g_r_AMPA_;
  ( *recordables )[ Name( "g_d_AMPA_" + std::to_string( syn_idx ) ) ] = &g_d_AMPA_;
  ( *recordables )[ Name( "i_AMPA_" + std::to_string( syn_idx ) ) ] = &i_AMPA_;
}

std::pair< double, double >
nest::AMPA::f_numstep( const double v_comp, const long lag )
{
  // update conductance
  g_r_AMPA_ *= prop_r_;
  g_d_AMPA_ *= prop_d_;

  // add spikes
  double s_val = b_spikes_->get_value( lag ) * g_norm_;
  g_r_AMPA_ -= s_val;
  g_d_AMPA_ += s_val;

  // compute synaptic conductance
  double g_AMPA = g_r_AMPA_ + g_d_AMPA_;

  // total current
  double i_tot = g_AMPA * ( e_rev_ - v_comp );
  // voltage derivative of total current
  double d_i_tot_dv = -g_AMPA;

  // for numerical integration
  double g_val = -d_i_tot_dv / 2.;
  double i_val = i_tot + g_val * v_comp;

  i_AMPA_ = i_val;

  if (debug){
    //printf("g_r_AMPA = %f     g_d_AMPA = %f     g_AMPA = %f\n",g_r_AMPA_,g_d_AMPA_,g_AMPA);
    //printf("i_tot = %f     e_rev = %f     v_comp = %f\n",i_tot,e_rev_,v_comp);
    //printf("g_val     = %f\n",g_val);
    if (tau_r_ == 5.)
      printf("i_Beta   = %f\n",i_val);
    else
      printf("i_BAP/AMPA    = %f\n",i_val);
    }
  return std::make_pair( g_val, i_val );
}


nest::GABA::GABA( const long syn_index )
  // initialization state variables
  : g_r_GABA_( 0.0 )
  , g_d_GABA_( 0.0 )
  // initialization parameters
  , e_rev_( -80. )
  , tau_r_( 0.2 )
  , tau_d_( 10.0 )
{
  syn_idx = syn_index;

  double tp = ( tau_r_ * tau_d_ ) / ( tau_d_ - tau_r_ ) * std::log( tau_d_ / tau_r_ );
  g_norm_ = 1. / ( -std::exp( -tp / tau_r_ ) + std::exp( -tp / tau_d_ ) );
}
nest::GABA::GABA( const long syn_index, const DictionaryDatum& receptor_params )
  // initialization state variables
  : g_r_GABA_( 0.0 )
  , g_d_GABA_( 0.0 )
  // initialization parameters
  , e_rev_( -80. )
  , tau_r_( 0.2 )
  , tau_d_( 10.0 )
{
  syn_idx = syn_index;

  // update GABA receptor parameters
  if ( receptor_params->known( "e_GABA" ) )
  {
    e_rev_ = getValue< double >( receptor_params, "e_GABA" );
  }
  if ( receptor_params->known( "tau_r_GABA" ) )
  {
    tau_r_ = getValue< double >( receptor_params, "tau_r_GABA" );
  }
  if ( receptor_params->known( "tau_d_GABA" ) )
  {
    tau_d_ = getValue< double >( receptor_params, "tau_d_GABA" );
  }

  double tp = ( tau_r_ * tau_d_ ) / ( tau_d_ - tau_r_ ) * std::log( tau_d_ / tau_r_ );
  g_norm_ = 1. / ( -std::exp( -tp / tau_r_ ) + std::exp( -tp / tau_d_ ) );
}

void
nest::GABA::append_recordables( std::map< Name, double* >* recordables )
{
  ( *recordables )[ Name( "g_r_GABA_" + std::to_string( syn_idx ) ) ] = &g_r_GABA_;
  ( *recordables )[ Name( "g_d_GABA_" + std::to_string( syn_idx ) ) ] = &g_d_GABA_;
}

std::pair< double, double >
nest::GABA::f_numstep( const double v_comp, const long lag )
{
  // update conductance
  g_r_GABA_ *= prop_r_;
  g_d_GABA_ *= prop_d_;

  // add spikes
  double s_val = b_spikes_->get_value( lag ) * g_norm_;
  g_r_GABA_ -= s_val;
  g_d_GABA_ += s_val;

  // compute synaptic conductance
  double g_GABA = g_r_GABA_ + g_d_GABA_;

  // total current
  double i_tot = g_GABA * ( e_rev_ - v_comp );
  // voltage derivative of total current
  double d_i_tot_dv = -g_GABA;

  // for numerical integration
  double g_val = -d_i_tot_dv / 2.;
  double i_val = i_tot + g_val * v_comp;

  if (debug)
    printf("i_GABA     = %f\n",i_val);

  return std::make_pair( g_val, i_val );
}


nest::NMDA::NMDA( const long syn_index )
  // initialization state variables
  : g_r_NMDA_( 0.0 )
  , g_d_NMDA_( 0.0 )
  , i_NMDA_(0.0)
  // initialization parameters
  , e_rev_( 0. )
  , tau_r_( 0.2 )
  , tau_d_( 43.0 )
{
  syn_idx = syn_index;

  double tp = ( tau_r_ * tau_d_ ) / ( tau_d_ - tau_r_ ) * std::log( tau_d_ / tau_r_ );
  g_norm_ = 1. / ( -std::exp( -tp / tau_r_ ) + std::exp( -tp / tau_d_ ) );
}
nest::NMDA::NMDA( const long syn_index, const DictionaryDatum& receptor_params )
  // initialization state variables
  : g_r_NMDA_( 0.0 )
  , g_d_NMDA_( 0.0 )
  , i_NMDA_(0.0)
  // initialization parameters
  , e_rev_( 0. )
  , tau_r_( 0.2 )
  , tau_d_( 43.0 )
{
  syn_idx = syn_index;

  // update NMDA receptor parameters
  if ( receptor_params->known( "e_NMDA" ) )
  {
    e_rev_ = getValue< double >( receptor_params, "e_NMDA" );
  }
  if ( receptor_params->known( "tau_r_NMDA" ) )
  {
    tau_r_ = getValue< double >( receptor_params, "tau_r_NMDA" );
  }
  if ( receptor_params->known( "tau_d_NMDA" ) )
  {
    tau_d_ = getValue< double >( receptor_params, "tau_d_NMDA" );
  }

  double tp = ( tau_r_ * tau_d_ ) / ( tau_d_ - tau_r_ ) * std::log( tau_d_ / tau_r_ );
  g_norm_ = 1. / ( -std::exp( -tp / tau_r_ ) + std::exp( -tp / tau_d_ ) );
}

void
nest::NMDA::append_recordables( std::map< Name, double* >* recordables )
{
  ( *recordables )[ Name( "g_r_NMDA_" + std::to_string( syn_idx ) ) ] = &g_r_NMDA_;
  ( *recordables )[ Name( "g_d_NMDA_" + std::to_string( syn_idx ) ) ] = &g_d_NMDA_;
  ( *recordables )[ Name( "i_NMDA_" + std::to_string( syn_idx ) ) ] = &i_NMDA_;
}

std::pair< double, double >
nest::NMDA::f_numstep( const double v_comp, const long lag )
{
  // update conductance
  g_r_NMDA_ *= prop_r_;
  g_d_NMDA_ *= prop_d_;

  // add spikes
  double s_val = b_spikes_->get_value( lag ) * g_norm_;
  g_r_NMDA_ -= s_val;
  g_d_NMDA_ += s_val;

  // compute conductance window
  double g_NMDA = g_r_NMDA_ + g_d_NMDA_;

  // auxiliary variables
  std::pair< double, double > NMDA_sigmoid = NMDA_sigmoid__and__d_NMDAsigmoid_dv( v_comp );

  // total current
  double i_tot = g_NMDA * NMDA_sigmoid.first * ( e_rev_ - v_comp );
  // voltage derivative of total current
  double d_i_tot_dv = g_NMDA * ( NMDA_sigmoid.second * ( e_rev_ - v_comp ) - NMDA_sigmoid.first );

  // for numerical integration
  double g_val = -d_i_tot_dv / 2.;
  double i_val = i_tot + g_val * v_comp;

  i_NMDA_ = i_val;

  if (debug)
    printf("i_NMDA     = %f\n",i_val);

  return std::make_pair( g_val, i_val );
}


nest::AMPA_NMDA::AMPA_NMDA( const long syn_index )
  // initialization state variables
  : g_r_AN_AMPA_( 0.0 )
  , g_d_AN_AMPA_( 0.0 )
  , g_r_AN_NMDA_( 0.0 )
  , g_d_AN_NMDA_( 0.0 )
  // initialization parameters
  , e_rev_( 0. )
  , tau_r_AMPA_( 0.2 )
  , tau_d_AMPA_( 3.0 )
  , tau_r_NMDA_( 0.2 )
  , tau_d_NMDA_( 43.0 )
  , NMDA_ratio_( 2.0 )
{
  syn_idx = syn_index;

  // AMPA normalization constant
  double tp = ( tau_r_AMPA_ * tau_d_AMPA_ ) / ( tau_d_AMPA_ - tau_r_AMPA_ ) * std::log( tau_d_AMPA_ / tau_r_AMPA_ );
  g_norm_AMPA_ = 1. / ( -std::exp( -tp / tau_r_AMPA_ ) + std::exp( -tp / tau_d_AMPA_ ) );
  // NMDA normalization constant
  tp = ( tau_r_NMDA_ * tau_d_NMDA_ ) / ( tau_d_NMDA_ - tau_r_NMDA_ ) * std::log( tau_d_NMDA_ / tau_r_NMDA_ );
  g_norm_NMDA_ = 1. / ( -std::exp( -tp / tau_r_NMDA_ ) + std::exp( -tp / tau_d_NMDA_ ) );
}
nest::AMPA_NMDA::AMPA_NMDA( const long syn_index, const DictionaryDatum& receptor_params )
  // initialization state variables
  : g_r_AN_AMPA_( 0.0 )
  , g_d_AN_AMPA_( 0.0 )
  , g_r_AN_NMDA_( 0.0 )
  , g_d_AN_NMDA_( 0.0 )
  // initialization parameters
  , e_rev_( 0. )
  , tau_r_AMPA_( 0.2 )
  , tau_d_AMPA_( 3.0 )
  , tau_r_NMDA_( 0.2 )
  , tau_d_NMDA_( 43.0 )
  , NMDA_ratio_( 2.0 )
{
  syn_idx = syn_index;

  // update AMPA+NMDA receptor parameters
  if ( receptor_params->known( "e_AMPA_NMDA" ) )
  {
    e_rev_ = getValue< double >( receptor_params, "e_AMPA_NMDA" );
  }
  if ( receptor_params->known( "tau_r_AMPA" ) )
  {
    tau_r_AMPA_ = getValue< double >( receptor_params, "tau_r_AMPA" );
  }
  if ( receptor_params->known( "tau_d_AMPA" ) )
  {
    tau_d_AMPA_ = getValue< double >( receptor_params, "tau_d_AMPA" );
  }
  if ( receptor_params->known( "tau_r_NMDA" ) )
  {
    tau_r_NMDA_ = getValue< double >( receptor_params, "tau_r_NMDA" );
  }
  if ( receptor_params->known( "tau_d_NMDA" ) )
  {
    tau_d_NMDA_ = getValue< double >( receptor_params, "tau_d_NMDA" );
  }
  if ( receptor_params->known( "NMDA_ratio" ) )
  {
    NMDA_ratio_ = getValue< double >( receptor_params, "NMDA_ratio" );
  }

  // AMPA normalization constant
  double tp = ( tau_r_AMPA_ * tau_d_AMPA_ ) / ( tau_d_AMPA_ - tau_r_AMPA_ ) * std::log( tau_d_AMPA_ / tau_r_AMPA_ );
  g_norm_AMPA_ = 1. / ( -std::exp( -tp / tau_r_AMPA_ ) + std::exp( -tp / tau_d_AMPA_ ) );
  // NMDA normalization constant
  tp = ( tau_r_NMDA_ * tau_d_NMDA_ ) / ( tau_d_NMDA_ - tau_r_NMDA_ ) * std::log( tau_d_NMDA_ / tau_r_NMDA_ );
  g_norm_NMDA_ = 1. / ( -std::exp( -tp / tau_r_NMDA_ ) + std::exp( -tp / tau_d_NMDA_ ) );
}

void
nest::AMPA_NMDA::append_recordables( std::map< Name, double* >* recordables )
{
  ( *recordables )[ Name( "g_r_AN_AMPA_" + std::to_string( syn_idx ) ) ] = &g_r_AN_AMPA_;
  ( *recordables )[ Name( "g_d_AN_AMPA_" + std::to_string( syn_idx ) ) ] = &g_d_AN_AMPA_;
  ( *recordables )[ Name( "g_r_AN_NMDA_" + std::to_string( syn_idx ) ) ] = &g_r_AN_NMDA_;
  ( *recordables )[ Name( "g_d_AN_NMDA_" + std::to_string( syn_idx ) ) ] = &g_d_AN_NMDA_;
}

std::pair< double, double >
nest::AMPA_NMDA::f_numstep( const double v_comp, const long lag )
{
  // update conductance
  g_r_AN_AMPA_ *= prop_r_AMPA_;
  g_d_AN_AMPA_ *= prop_d_AMPA_;
  g_r_AN_NMDA_ *= prop_r_NMDA_;
  g_d_AN_NMDA_ *= prop_d_NMDA_;

  // add spikes
  double s_val_ = b_spikes_->get_value( lag );
  double s_val = s_val_ * g_norm_AMPA_;
  g_r_AN_AMPA_ -= s_val;
  g_d_AN_AMPA_ += s_val;
  s_val = s_val_ * g_norm_NMDA_;
  g_r_AN_NMDA_ -= s_val;
  g_d_AN_NMDA_ += s_val;

  // compute conductance window
  double g_AMPA = g_r_AN_AMPA_ + g_d_AN_AMPA_;
  double g_NMDA = g_r_AN_NMDA_ + g_d_AN_NMDA_;

  // auxiliary variable
  std::pair< double, double > NMDA_sigmoid = NMDA_sigmoid__and__d_NMDAsigmoid_dv( v_comp );

  // total current
  double i_tot = ( g_AMPA + NMDA_ratio_ * g_NMDA * NMDA_sigmoid.first ) * ( e_rev_ - v_comp );
  // voltage derivative of total current
  double d_i_tot_dv =
    -g_AMPA + NMDA_ratio_ * g_NMDA * ( NMDA_sigmoid.second * ( e_rev_ - v_comp ) - NMDA_sigmoid.first );

  // for numerical integration
  double g_val = -d_i_tot_dv / 2.;
  double i_val = i_tot + g_val * v_comp;

  if (debug)
    printf("i_AMPA_NMDA     = %f\n",i_val);

  return std::make_pair( g_val, i_val );
}




////////////
// Back_I start
////////////

nest::Back_I::Back_I( const long syn_index )
  // initialization state variables
  : g_Back_I_( 0.0 )
  // initialization parameters
  , e_rev_( 0.0 )
  , g_max_( 0.0 )
  , tau_( 1.0 )
{
  syn_idx = syn_index;

  g_norm_ = g_max_ / tau_ * std::exp( 1 );
}
nest::Back_I::Back_I( const long syn_index, const DictionaryDatum& receptor_params )
  // initialization state variables
  : g_Back_I_( 0.0 )
  // initialization parameters
  , e_rev_( 0.0 )
  , g_max_( 0.0 )
  , tau_( 1.0 )
{
  syn_idx = syn_index;

  // update Back_I receptor parameters
  if ( receptor_params->known( "e_Back_I" ) )
  {
    e_rev_ = getValue< double >( receptor_params, "e_Back_I" );
  }
  if ( receptor_params->known( "g_max" ) )
  {
    g_max_ = getValue< double >( receptor_params, "g_max" );
  }
  if ( receptor_params->known( "tau_Back_I" ) )
  {
    tau_ = getValue< double >( receptor_params, "tau_Back_I" );
  }

  g_norm_ = g_max_ / tau_ * std::exp( 1 );
}

void
nest::Back_I::append_recordables( std::map< Name, double* >* recordables )
{
  ( *recordables )[ Name( "g_Back_I_" + std::to_string( syn_idx ) ) ] = &g_Back_I_;
}

std::pair< double, double >
nest::Back_I::f_numstep( const double v_comp, const long lag )
{
  // update conductance
  g_Back_I_ *= prop_;

  // add spikes
  double s_val = b_spikes_->get_value( lag ) * g_norm_;
  g_Back_I_ -= s_val;

  // compute synaptic conductance
  double g_Back_I = g_Back_I_;

  // total current
  double i_tot = g_Back_I * ( e_rev_ - v_comp );
  // voltage derivative of total current
  double d_i_tot_dv = -g_Back_I;

  // for numerical integration
  double g_val = -d_i_tot_dv / 2.;
  double i_val = i_tot + g_val * v_comp;

  if (debug)
    printf("i_Back (not used!)     = %f\n",i_val);

  return std::make_pair( g_val, i_val );
}

////////////
// Back_I end
////////////


nest::DC::DC( const long syn_index )
  // initialization state variables
  : g_r_DC_( 0.0 )
  , g_d_DC_( 0.0 )
  // initialization parameters
  , e_rev_( 0.0 )
  , tau_r_( 0.2 )
  , tau_d_( 3.0 )
{
  syn_idx = syn_index;

  double tp = ( tau_r_ * tau_d_ ) / ( tau_d_ - tau_r_ ) * std::log( tau_d_ / tau_r_ );
  g_norm_ = 1. / ( -std::exp( -tp / tau_r_ ) + std::exp( -tp / tau_d_ ) );
  g_norm_ = 1.;
}
nest::DC::DC( const long syn_index, const DictionaryDatum& receptor_params )
  // initialization state variables
  : g_r_DC_( 0.0 )
  , g_d_DC_( 0.0 )
  // initialization parameters
  , e_rev_( 0.0 )
  , tau_r_( 0.2 )
  , tau_d_( 3.0 )
{
  syn_idx = syn_index;

  // update DC receptor parameters
  if ( receptor_params->known( "e_DC" ) )
  {
    e_rev_ = getValue< double >( receptor_params, "e_DC" );
  }
  if ( receptor_params->known( "tau_r_DC" ) )
  {
    tau_r_ = getValue< double >( receptor_params, "tau_r_DC" );
  }
  if ( receptor_params->known( "tau_d_DC" ) )
  {
    tau_d_ = getValue< double >( receptor_params, "tau_d_DC" );
  }

  double tp = ( tau_r_ * tau_d_ ) / ( tau_d_ - tau_r_ ) * std::log( tau_d_ / tau_r_ );
  g_norm_ = 1. / ( -std::exp( -tp / tau_r_ ) + std::exp( -tp / tau_d_ ) );
  g_norm_ = 1.;
}

void
nest::DC::append_recordables( std::map< Name, double* >* recordables )
{
  ( *recordables )[ Name( "g_r_DC_" + std::to_string( syn_idx ) ) ] = &g_r_DC_;
  ( *recordables )[ Name( "g_d_DC_" + std::to_string( syn_idx ) ) ] = &g_d_DC_;
}

std::pair< double, double >
nest::DC::f_numstep( const double v_comp, const long lag )
{
  double g_val = 0.;
  double i_val = 0.;

  // update conductance
  //g_r_DC_ *= prop_r_;
  //g_d_DC_ *= prop_d_;

  // add spikes
  double s_val = b_spikes_->get_value( lag ) * g_norm_;
  const double dt = Time::get_resolution().get_ms();

  if(s_val != 0.)
    {
      up_time_ = 0.1;
      dc_permanent_ = s_val/2;
    }
  else
    {
      if(up_time_ != 0.)
	{
	  up_time_ += dt;
	  if(up_time_ > 5.)
	    {
	      up_time_ = 0.;
	      dc_permanent_ = 0.;
	    }
	}
    }

  i_val = dc_permanent_;

  //printf("up_time = %f     i_val = %f\n",up_time_,dc_permanent_);

  //g_r_DC_ -= s_val;
  //g_d_DC_ += s_val;

  // compute synaptic conductance
  //double g_DC = g_r_DC_ + g_d_DC_;

  // total current
  //double i_tot = g_DC * ( e_rev_ - v_comp );
  // voltage derivative of total current
  //double d_i_tot_dv = -g_DC;

  // for numerical integration
  //double g_val = -d_i_tot_dv / 2.;
  //double i_val = i_tot + g_val * v_comp;

  if (debug)
    printf("i_DC (not used!)     = %f\n",i_val);

  return std::make_pair( g_val, i_val );
}





nest::AHP::AHP( const long syn_index )
  // initialization state state variables
  : Ca_conc_AHP_( 0.0 )
  // initialization parameters
  , g_AHP_( 0.0 )
  , e_K_AHP_( -90.0 )
  , tau_K_AHP_( 80.0 )
{
  syn_idx = syn_index;

}
nest::AHP::AHP( const long syn_index, const DictionaryDatum& receptor_params )
  // initialization state state variables
  : Ca_conc_AHP_( 0.0 )
  // initialization parameters
  , g_AHP_( 0.0 )
  , e_K_AHP_( -90.0 )
  , tau_K_AHP_( 80.0 )
{
  syn_idx = syn_index;

  // update AHP receptor parameters
  if ( receptor_params->known( "g_AHP" ) )
  {
    g_AHP_ = getValue< double >( receptor_params, "g_AHP" );
  }
  if ( receptor_params->known( "e_K_AHP" ) )
  {
    e_K_AHP_ = getValue< double >( receptor_params, "e_K_AHP" );
  }
  if ( receptor_params->known( "tau_K_AHP" ) )
  {
    tau_K_AHP_ = getValue< double >( receptor_params, "tau_K_AHP" );
  }
}

void
nest::AHP::append_recordables( std::map< Name, double* >* recordables )
{
  ( *recordables )[ Name( "Ca_conc_AHP_" + std::to_string( syn_idx ) ) ] = &Ca_conc_AHP_;
}

std::pair< double, double >
nest::AHP::f_numstep( const double v_comp, const long lag )
{
  double g_val = 0.;
  double i_val = 0.;

  double e_L = -70.6;       // [mV]
  double a_ = 0.0;          // [1/mV]

  // add spikes
  double s_val = b_spikes_->get_value( lag );
  const double dt = Time::get_resolution().get_ms();

  // advance state variable 'Ca_conc' one timestamp  (Larkum-Senn 2004)
  Ca_conc_AHP_ *= std::exp( -dt/ tau_K_AHP_);
  if (s_val != 0.)
    Ca_conc_AHP_ += s_val;
  // additional term "looking for" compatibility with adex
  Ca_conc_AHP_ += dt/tau_K_AHP_ * a_ * (e_L - v_comp);

  // add to variables for numerical integration
  g_val += Ca_conc_AHP_ * g_AHP_ / 2.;
  i_val += Ca_conc_AHP_ * g_AHP_ * ( e_K_AHP_ - v_comp / 2.);

  if (debug)
    printf("i_AHP    = %f     Ca_conc_AHP_ = %f\n",i_val,Ca_conc_AHP_);

  return std::make_pair( g_val, i_val );
}




nest::REFRACT::REFRACT( const long syn_index )
  // initialization state variables
  : t_ref_( 0.0 )
  , V_reset_(0.0)
  ,g_refract_(100000)
  // initialization parameters
  , t_ref_residual_( 0.0 )
  , V_refract_( 0.0 )
{
  syn_idx = syn_index;
}

nest::REFRACT::REFRACT( const long syn_index, const DictionaryDatum& receptor_params )
  // initialization state variables
  : t_ref_( 0.0 )
  , V_reset_(0.0)
  , g_refract_(100000)
  // initialization parameters
  , t_ref_residual_( 0.0 )
  , V_refract_( 0.0 )
{
  syn_idx = syn_index;

  // update REFRACT receptor parameters
  if ( receptor_params->known( "t_ref" ) )
  {
    t_ref_ = getValue< double >( receptor_params, "t_ref" );
  }
  if ( receptor_params->known( "V_reset" ) )
  {
    V_reset_ = getValue< double >( receptor_params, "V_reset" );
  }
  if ( receptor_params->known( "g_refract" ) )
  {
    g_refract_ = getValue< double >( receptor_params, "g_refract" );
  }
}

void
nest::REFRACT::append_recordables( std::map< Name, double* >* recordables )
{
  ( *recordables )[ Name( "t_ref_residual_" + std::to_string( syn_idx ) ) ] = &t_ref_residual_;
}

std::pair< double, double >
nest::REFRACT::f_numstep( const double v_comp, const long lag )
{
  double g_val = 0.0;
  double i_val = 0.0;

  // get spike value
  double s_val = b_spikes_->get_value( lag );
  const double dt = Time::get_resolution().get_ms();

  assert( !(s_val != 0.0 &&  t_ref_residual_ > 0.0) );

  if (s_val != 0.)
    {
      t_ref_residual_ = t_ref_;
      V_refract_ = v_comp;
    }
  if (t_ref_residual_ > 0)
    {
      t_ref_residual_ = t_ref_residual_ - dt;
      g_val = g_refract_/2;
      //i_val = g_refract_ * (V_reset_ - v_comp/2);
      i_val = g_refract_ * (V_refract_ - v_comp/2);
    }

  if (debug)
    printf("i_REFRACT     = %f\n",i_val);

  return std::make_pair( g_val, i_val );
}






nest::ADAPT::ADAPT( const long syn_index )
  // initialization state state variables
  : w_( 0.0 )
  // initialization parameters
  , tau_w_( 80.0 )
  , e_w_( -90.0 )
  , a_( 0.0 )
  , b_( 0.0 )
{
  syn_idx = syn_index;

}
nest::ADAPT::ADAPT( const long syn_index, const DictionaryDatum& receptor_params )
  // initialization state state variables
  : w_( 0.0 )
  // initialization parameters
  , tau_w_( 80.0 )
  , e_w_( -90.0 )
  , a_( 0.0 )
  , b_( 0.0 )
{
  syn_idx = syn_index;

  // update AHP receptor parameters
  if ( receptor_params->known( "tau_w" ) )
  {
    tau_w_ = getValue< double >( receptor_params, "tau_w" );
  }
  if ( receptor_params->known( "e_w" ) )
  {
    e_w_ = getValue< double >( receptor_params, "e_w" );
  }
  if ( receptor_params->known( "a" ) )
  {
    a_ = getValue< double >( receptor_params, "a" );
  }
  if ( receptor_params->known( "b" ) )
  {
    b_ = getValue< double >( receptor_params, "b" );
  }
}

void
nest::ADAPT::append_recordables( std::map< Name, double* >* recordables )
{
  ( *recordables )[ Name( "w_" + std::to_string( syn_idx ) ) ] = &w_;
}

std::pair< double, double >
nest::ADAPT::f_numstep( const double v_comp, const long lag )
{
  double g_val = 0.;
  double i_val = 0.;
  
  double E_K_ = -200;

  // add spikes
  double s_val = b_spikes_->get_value( lag );
  const double dt = Time::get_resolution().get_ms();
  
  // adaptation current
  w_ += dt/tau_w_ * (a_ * (v_comp - e_w_) - w_);
  if (s_val != 0.)
    w_ += b_;
  
  if( (v_comp-E_K_) > 0.001)
    {
      g_val += w_ / (2 * fabs(v_comp - E_K_)) ;
      i_val += w_ / fabs(v_comp - E_K_) * (E_K_ - v_comp/2);
    }
  else
    {
      g_val += 0.;
      i_val += 0.;
    }


  if (debug)
    printf("i_ADAPT    = %f\n",i_val);
  
  return std::make_pair( g_val, i_val );
}







nest::BETA::BETA( const long syn_index )
  // initialization state variables
  : i_r_BETA_( 0.0 )
  , i_d_BETA_( 0.0 )
  // initialization parameters
  , i_input_( 0.0 )
  , tau_r_( 0.2 )
  , tau_d_( 3.0 )
{
  syn_idx = syn_index;

  //double tp = ( tau_r_ * tau_d_ ) / ( tau_d_ - tau_r_ ) * std::log( tau_d_ / tau_r_ );
  //i_norm_ = 1. / ( -std::exp( -tp / tau_r_ ) + std::exp( -tp / tau_d_ ) );
  tau_c_ = tau_r_ * tau_d_ / (tau_r_ + tau_d_);
  double t_max_ = tau_d_*tau_c_/(tau_c_-tau_d_)*std::log(tau_c_/tau_d_);
  beta_norm_ = std::exp(-t_max_/tau_d_)-std::exp(-t_max_/tau_c_);
}
nest::BETA::BETA( const long syn_index, const DictionaryDatum& receptor_params )
  // initialization state variables
  : i_r_BETA_( 0.0 )
  , i_d_BETA_( 0.0 )
  // initialization parameters
  , i_input_( 0.0 )
  , tau_r_( 0.2 )
  , tau_d_( 3.0 )
{
  syn_idx = syn_index;

  // update BETA receptor parameters
  //if ( receptor_params->known( "i_input" ) )
  //{
  //  i_input_ = getValue< double >( receptor_params, "i_input" );
  //}
  if ( receptor_params->known( "tau_r_BETA" ) )
  {
    tau_r_ = getValue< double >( receptor_params, "tau_r_BETA" );
  }
  if ( receptor_params->known( "tau_d_BETA" ) )
  {
    tau_d_ = getValue< double >( receptor_params, "tau_d_BETA" );
  }

  //double tp = ( tau_r_ * tau_d_ ) / ( tau_d_ - tau_r_ ) * std::log( tau_d_ / tau_r_ );
  //i_norm_ = 1. / ( -std::exp( -tp / tau_r_ ) + std::exp( -tp / tau_d_ ) );
  tau_c_ = tau_r_ * tau_d_ / (tau_r_ + tau_d_);
  double t_max_ = tau_d_*tau_c_/(tau_c_-tau_d_)*std::log(tau_c_/tau_d_);
  beta_norm_ = std::exp(-t_max_/tau_d_)-std::exp(-t_max_/tau_c_);
}

void
nest::BETA::append_recordables( std::map< Name, double* >* recordables )
{
  ( *recordables )[ Name( "i_r_BETA_" + std::to_string( syn_idx ) ) ] = &i_r_BETA_;
  ( *recordables )[ Name( "i_d_BETA_" + std::to_string( syn_idx ) ) ] = &i_d_BETA_;
  ( *recordables )[ Name( "i_BETA_" + std::to_string( syn_idx ) ) ] = &i_BETA_;
}

std::pair< double, double >
nest::BETA::f_numstep( const double v_comp, const long lag )
{
  // update conductance
  i_r_BETA_ *= prop_r_;
  i_d_BETA_ *= prop_d_;

  // add input currents
  //double s_val = b_spikes_->get_value( lag ) * i_norm_;
  double s_val = b_spikes_->get_value( lag )/beta_norm_;
  i_r_BETA_ -= s_val;
  i_d_BETA_ += s_val;

  // compute total input current
  i_BETA_ = i_r_BETA_ + i_d_BETA_;

  // for numerical integration
  double g_val = i_BETA_ / (2 * fabs(e_rev_ - v_comp));
  double i_val = i_BETA_ / (2 * fabs(e_rev_ - v_comp)) * (e_rev_ - v_comp/2);

  if (debug){
    //printf("i_r_BETA = %f     i_d_BETA = %f     i_BETA = %f\n",i_r_BETA_,i_d_BETA_,i_BETA_);
    //printf("i_tot = %f     e_rev = %f     v_comp = %f\n",i_tot,e_rev_,v_comp);
    //printf("g_val     = %f\n",g_val);
    if (tau_r_ == 5.)
      printf("i_Beta   = %f\n",i_val);
    else
      printf("i_BAP/BETA    = %f\n",i_val);
    }
  return std::make_pair( g_val, i_val );
}
