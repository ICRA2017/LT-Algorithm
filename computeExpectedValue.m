%========================================================================
%   dispersion_index
%   version 1.0 - January 18th, 2017
%
%   Computes the expected value according to the input PDF
%   Inputs:
%   dir_mesh,speed_mesh 2D matrices with the wind speed/direction values
%   pdf_mat: input pdf matrix
%
%   NOTE: This script uses the circStat toolbox developed by Philipp
%   Berens and co-authors. The circStat toolbox used in this work was 
%   downloaded from:  
%   https://se.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox--directional-statistics-
%
%========================================================================
function [expected_dir expected_speed]=computeExpectedValue(dir_mesh,speed_mesh,pdf_mat)

posterior_vector=pdf_mat(:);
dir_vector=dir_mesh(:);
speed_vector=speed_mesh(:);

expected_speed=sum(speed_vector.*posterior_vector);
expected_dir=wrapTo360(rad2deg(circ_mean(deg2rad(dir_vector),posterior_vector)));
