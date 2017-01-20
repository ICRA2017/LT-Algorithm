%========================================================================
%   dispersion_index
%   version 1.0 - January 18th, 2017
%
%   Computes an ad-hoc dispersion measurement based on a logistic function
%   Inputs:
%   dir_mesh (varargin{1}): a 2D matrices with the wind speed values 
%   turn_point(varargin{2}): mid point of the logistic function
%   slope(varargin{3}): growth rate of the logistic function
%   t_pdf (varargin{4}): Turbulent extrapolation matrix (equation 5b)
%   
%   NOTE: This script uses the circStat toolbox developed by Philipp
%   Berens and co-authors. The circStat toolbox used in this work was 
%   downloaded from:  
%   https://se.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox--directional-statistics-
%
%========================================================================

function disp_meas=dispersion_index(varargin)

% if length(varargin)>3
%     weight_vector=varargin{4};
%     weight_vector=weight_vector(:);
%     weight_vector=weight_vector./sum(weight_vector);
% end

angles=varargin{1};
angles=angles(:);
angles=wrapTo360(angles);
Ro=varargin{2};
K=varargin{3};
weight_vector=varargin{4};
weight_vector=weight_vector(:);
% S=sum(sind(angles).*weight_vector);
% C=sum(cosd(angles).*weight_vector);
R=circ_var(deg2rad(angles),weight_vector);
%R=1-sqrt(S*S+C*C);
disp_meas=1/(1+exp(-K*(R-Ro)));


