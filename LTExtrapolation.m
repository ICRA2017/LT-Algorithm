%========================================================================
%   LTExtrapolation
%   version 1.0 - January 18th, 2017
%   
%   This script computes the Laminar-Turbulent (LT) extrapolation as presented in the
%   article  "Probabilistic Air Flow Modelling Using Turbulent and Laminar 
%   Characteristics for Ground and Aerial Robots", according to equations 5a and 5c 
%   in the article. 
%   
%   Inputs:
%   dir_mesh: a 2D matrices with the wind speed values 
%   for each state in Gamma (state space)
%   dir_states, speed_states: vectors that contain the possible wind speed/wind  
%   direction values of the discrete joint PDFs
%   sensor_model: object that defines the sensor model according to
%   equation 2
%   dir_std,speed_std: Parameters of the sensor model (equation 2)
%   dir_bin_size,speed_bin_size: Discretization parameters for the state
%   space.
%   l_pdf: A matrix of size A X B, where A X B is
%   the size of the state space. l_pdf corresponds to the laminar
%   extrapolation 
%   t_pdf: A matrix of size A X B, where A X B is
%   the size of the state space. l_pdf corresponds to the turbulent
%   extrapolation 
%   slope, turn point: growth rate and mid point parameters of the logistic
%   function used to estimate the turbulece index (kappa in equation 5c)
%========================================================================



function [weight_pdf,turbulence_idx]=LTExtrapolation(l_pdf,t_pdf,dir_mesh,slope,turn_point)

% computes the turbulence index
turbulence_idx=dispersion_index(dir_mesh,turn_point,slope,t_pdf);
weight_pdf=turbulence_idx.*t_pdf+(1-turbulence_idx).*l_pdf;
weight_pdf=weight_pdf./sum(weight_pdf(:));