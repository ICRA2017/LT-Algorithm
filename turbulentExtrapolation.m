%========================================================================
%   turbulentExtrapolation
%   version 1.0 - January 18th, 2017
%   
%   This script computes the turbulent extrapolation as presented in the
%   article  "Probabilistic Air Flow Modelling Using Turbulent and Laminar 
%   Characteristics for Ground and Aerial Robots", according to equation 5b
%   and 5c in the article. 
%   Inputs:
%   pdf_vectors: A matrix of size A X B X Nm, where A X B is
%   the size of the state space and Nm is the number of
%   measurement locations. pdf_vectors corresponds to the computed joint
%   pdfs at each measurement location visited by the mobile sensors (i.e.
%   robot) 
%   pdf_locations: this is a Nm X 2 vector that contains the (x,y)
%   coordinates of the positions at which the joint pdfs in pdf_vectors
%   were computed.
%   query_location: x,y location of the query point (in m)
%   kernel_size: size (in m) of the Gaussian kernel for extrapolation (equations 4d) 
%   n_samples: Number of artificial samples used in the extrapolation process (Equation 5b) 
%========================================================================


function [mixed_pdf]=turbulentExtrapolation(pdf_vector_in,pdf_locations,query_location,kernel_size)
pos_x_in=pdf_locations(:,1);
pos_y_in=pdf_locations(:,2);

sel_pos_x=query_location(1);
sel_pos_y=query_location(2);

[n_rows,n_cols]=size(squeeze(pdf_vector_in(:,:,1)));
mixed_pdf=zeros(n_rows,n_cols);
    sum_weights=0;


for i=1:length(pos_y_in)
    x_pos=pos_x_in(i);
    y_pos=pos_y_in(i);
    
    distance=sqrt((sel_pos_x-x_pos).^2+(sel_pos_y-y_pos).^2);
    norm_fact=1/(kernel_size*sqrt(2*pi));
    exp_fact=1/(2*kernel_size^2);
    kernel_coef=norm_fact.*exp(-distance.^2*exp_fact);
    mixed_pdf=mixed_pdf+pdf_vector_in(:,:,i).*kernel_coef;
    sum_weights=sum_weights+kernel_coef;
end
% mixed_pdf=mixed_pdf./sum_weights;
mixed_pdf=mixed_pdf./sum(mixed_pdf(:));