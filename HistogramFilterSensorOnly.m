%========================================================================
%   HistogramFilterSensorOnly
%   version 1.0 - January 18th, 2017
%
%   Computes a histogram filter that considers only the sensor model
%   Inputs: 
%   dir_mesh,speed_mesh 2D matrices with the wind speed/direction values 
%   for each state in Gamma (state space)
%   sensor_model: object that defines the sensor model according to
%   equation 2
%   dir_std,speed_std: Parameters of the sensor model (equation 2)
%   dir_bin_size,speed_bin_size: Discretization parameters for the state
%   space.
%   wind data: Nx2 Measurement matrix that contains values of wind direction 
%   (degrees) in the first columd and wind speed (m/s) in the second column
%========================================================================


function [sens_posteriors]=HistogramFilterSensorOnly(dir_mesh,speed_mesh,sensor_model,dir_std,speed_std,dir_bin_size,speed_bin_size,wind_data)


%---------------------------------------
% Data initialization
%--------------------------------------
n_data=length(wind_data);
[n_rows,n_cols]=size(dir_mesh);

state_accummulator_sens=zeros(n_rows,n_cols);

for i=1:n_data
  
 
    
    dir_val=wind_data(i,1);
    speed_val=wind_data(i,2);
    
    dir_vector=reshape(dir_mesh,[n_rows*n_cols,1]);
    speed_vector=reshape(speed_mesh,[n_rows*n_cols,1]);
    
    distance_vector_dir=abs(rad2deg(circ_dist(deg2rad(dir_vector),deg2rad(dir_val))));
    distance_vector_speed=abs(speed_val-speed_vector);

    distance_matrix_dir=reshape(distance_vector_dir,n_rows,n_cols);
    distance_matrix_speed=reshape(distance_vector_speed,n_rows,n_cols);
    
   idx_dir_mat=distance_matrix_dir<=3*dir_std;
   idx_speed_mat=distance_matrix_speed<=3*speed_std; 
   idx_mat=idx_dir_mat.*idx_speed_mat+0;
   [idx_r, idx_c]=find(idx_mat>0);  
   
 %  PX=zeros(n_rows,n_cols)*(1/(n_rows*n_cols));
   for k=1:length(idx_r)
       i_r=idx_r(k);
       i_c=idx_c(k);
       bin_center_dir=dir_mesh(i_r,i_c);
       bin_center_speed=speed_mesh(i_r,i_c);
       dir_prior=sensor_model.dirPrior(dir_val,wrapTo360(bin_center_dir-dir_bin_size/2),wrapTo360(bin_center_dir+dir_bin_size/2));
       speed_prior=sensor_model.speedPrior(speed_val,bin_center_speed-speed_bin_size/2,bin_center_speed+speed_bin_size/2);
       state_accummulator_sens(i_r,i_c)=state_accummulator_sens(i_r,i_c)+dir_prior*speed_prior;
   end
   %sens_posteriors=state_accummulator_sens/sum(state_accummulator_sens(:));
    
end
sens_posteriors=state_accummulator_sens/sum(state_accummulator_sens(:));




