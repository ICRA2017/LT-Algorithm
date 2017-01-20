%========================================================================
%   laminarExtrapolation
%   version 1.0 - January 18th, 2017
%   
%   This script computes the laminar extrapolation as presented in the
%   article  "Probabilistic Air Flow Modelling Using Turbulent and Laminar 
%   Characteristics for Ground and Aerial Robots", according to equation 4
%   in the article. The laminar extrapolation is based on a Nadaraya-Watson
%   estimatior.
%   
%   Inputs:
%   dir_mesh,speed_mesh 2D matrices with the wind speed/direction values 
%   for each state in Gamma (state space)
%   dir_states, speed_states: vectors that contain the possible wind speed/wind  
%   direction values of the discrete joint PDFs
%   sensor_model: object that defines the sensor model according to
%   equation 2
%   dir_std,speed_std: Parameters of the sensor model (equation 2)
%   dir_bin_size,speed_bin_size: Discretization parameters for the state
%   space.
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


function [out_pdf_vector]=laminarExtrapolation(dir_mesh,speed_mesh,dir_states,sensor_model,speed_states,dir_std,speed_std,dir_bin_size,speed_bin_size,pdf_vectors,pdf_locations,query_location,kernel_size,n_samples)

[n_rows,n_cols,n_pdfs]=size(pdf_vectors);


x_locations=pdf_locations(:,1);
y_locations=pdf_locations(:,2);

extrapolated_dir_vector=zeros(n_samples,1);
extrapolated_speed_vector=zeros(n_samples,1);

dir_samp_vector=zeros(n_pdfs,n_samples);
speed_samp_vector=zeros(n_pdfs,n_samples);

% Acquires artificial samples from the joint PDFS
for i=1:n_pdfs
    [dir_v,speed_v]=samplePdf(dir_states,speed_states,dir_bin_size,speed_bin_size,pdf_vectors(:,:,i),n_samples);
    dir_samp_vector(i,:)=dir_v;
    speed_samp_vector(i,:)=speed_v;
end


% Performs Naradaya-Watson extrapolation using the artificial samples
% (equations 4b to 4d)
for j=1:n_samples
    query_x=query_location(1);
    query_y=query_location(2);
    
    norm_fact=1/(kernel_size*sqrt(2*pi));
    exp_fact=1/(2*kernel_size^2);
    
    distance_vector=sqrt((query_x-x_locations).^2+(query_y-y_locations).^2);
    kernel_coefs=norm_fact.*exp(-distance_vector.^2.*exp_fact);
    
    weight_vector=sum(kernel_coefs);
    
    extrapolated_speed=sum(speed_samp_vector(:,j).*kernel_coefs)/weight_vector;
    extrapolated_dir=wrapTo360(rad2deg(circ_mean(deg2rad(dir_samp_vector(:,j)),kernel_coefs)));
    extrapolated_dir_vector(j)=extrapolated_dir;
    extrapolated_speed_vector(j)=extrapolated_speed;
end

wind_data=[extrapolated_dir_vector,extrapolated_speed_vector];

% Using the extrapolated artificial samples, a histogram that considers the
% sensor model is computed (equation 4a)
[out_pdf_vector]=HistogramFilterSensorOnly(dir_mesh,speed_mesh,sensor_model,dir_std,speed_std,dir_bin_size,speed_bin_size,wind_data);
