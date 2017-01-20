%========================================================================
%   generateMatrices
%   version 1.0 - January 18th, 2017
%   
%   Initializes the matrices and vectors to be used in the computation of
%   the wind speed/wind direction joint PDFs and in the computation of the
%   LT algorihtm. 
%   input parameters:
%       dir_bin_size - wind direction bin size (degrees)
%       speed_bin_size - speed bin size (m/s)
%       speed_range - min and max values of the wind speed to be considered in the modeled PDF
% 
%   outputs:
%       dir_mesh, speed_mesh - rectangular grids that contain the state
%       tuples (wind direction,wind speed) for each state of the discrete joint PDF 
%       dir_states, speed_states - vectors that contain the possible wind speed/wind
%       direction values of the discrete joint PDFs
%========================================================================


function [dir_mesh,speed_mesh,dir_states,speed_states]=generateMatrices(dir_bin_size,speed_bin_size,speed_range)

n_dir=round((359.9)/dir_bin_size);
dir_states=linspace(0,359.9,n_dir);

n_speed=round((speed_range(2)-speed_range(1))/speed_bin_size);
speed_states=linspace(speed_range(1),speed_range(2),n_speed);

[dir_mesh,speed_mesh]=meshgrid(dir_states,speed_states);


