%========================================================================
%   normalizeAngles
%   version 1.0 - January 18th, 2017
%   
%   Function that takes a vector of n given angular values, wrap them 
%   between [0,360] and returns them in a vector of nx1
%========================================================================

function output_values=normalizeAngles(input_values)
output_values=input_values(:);
output_values=wrapTo360(output_values);
output_values(output_values==360)=0;