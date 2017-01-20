function [state_posteriors,sens_posteriors,histogram_out]=jointPDFBayes(dir_mesh,speed_mesh,dir_states,speed_states,sensor_model,dir_std,speed_std,dir_bin_size,speed_bin_size,wind_data)

%========================================================================
%   jointPDFBayes
%   version 1.0 - January 18th, 2017
%   
%   This script computes the joint PDFs of the airflow (wind speed/wind
%   direction) at a given location. For the computations of the posteriors,
%   the script uses equations 1 to 3 from the article "Probabilistic Air 
%   Flow Modelling Using Turbulent and Laminar Characteristics for Ground 
%   and Aerial Robots".
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
%   wind data: Nx2 Measurement matrix that contains values of wind direction 
%   (degrees) in the first columd and wind speed (m/s) in the second column
%   Outputs:
%   histogram_out: Simple 2D histogram computed with the raw input
%   state_posteriors: Joint PDFs computed with the Bayesian filter
%   sens_posteriors: Joint PDFs computed using the sensor model ONLY
%   
%   NOTE: This script uses the circStat toolbox developed by Philipp
%   Berens and co-authors. The circStat toolbox used in this work was 
%   downloaded from:  
%   https://se.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox--directional-statistics-
%
%========================================================================




%---------------------------------------
% Data initialization
%--------------------------------------
n_data=length(wind_data);
[n_rows,n_cols]=size(dir_mesh);
state_posteriors=(1/(n_rows*n_cols)).*ones(n_rows,n_cols);
state_accummulator=zeros(n_rows,n_cols);
state_accummulator_sens=zeros(n_rows,n_cols);

centers={dir_states,speed_states};
[N,C]=hist3(wind_data,centers);
histogram_out=N'./sum(N(:));

PX=ones(n_rows,n_cols).*(1/(n_rows*n_cols));
transition_matrices=zeros(n_rows,n_cols,n_rows,n_cols);

dir_vector=reshape(dir_mesh,[n_rows*n_cols,1]);
speed_vector=reshape(speed_mesh,[n_rows*n_cols,1]);


for i=1:n_data
    %---------------------------------------
    % Updates the transition matrix
    %--------------------------------------
    
    if i==1
        prev_dir_val=wind_data(1,1);
        prev_speed_val=wind_data(1,2);
        
        distance_dir_vector=abs(rad2deg(circ_dist(deg2rad(dir_states),deg2rad(prev_dir_val))));
        distance_speed_vector=abs(prev_speed_val-speed_states);
        
        dir_idx=find(distance_dir_vector==min(distance_dir_vector));
        previous_dir_idx=dir_idx(1);
        
        speed_idx=find(distance_speed_vector==min(distance_speed_vector));
        previous_speed_idx=speed_idx(1);
    else
        dir_val=wind_data(i,1);
        speed_val=wind_data(i,2);
        
        distance_vector_dir=abs(rad2deg(circ_dist(deg2rad(dir_vector),deg2rad(dir_val))));
        distance_vector_speed=abs(speed_val-speed_vector);
        
        distance_matrix_dir=reshape(distance_vector_dir,n_rows,n_cols);
        distance_matrix_speed=reshape(distance_vector_speed,n_rows,n_cols);
        
        idx_dir_mat=distance_matrix_dir<=3*dir_std;
        idx_speed_mat=distance_matrix_speed<=3*speed_std;
        idx_mat=idx_dir_mat.*idx_speed_mat;
        
        [idx_r, idx_c]=find(idx_mat>0);
        
        for k=1:length(idx_r)
            i_r=idx_r(k);
            i_c=idx_c(k);
            bin_center_dir=dir_mesh(i_r,i_c);
            bin_center_speed=speed_mesh(i_r,i_c);
            dir_prior=sensor_model.dirPrior(dir_val,wrapTo360(bin_center_dir-dir_bin_size/2),wrapTo360(bin_center_dir+dir_bin_size/2));
            speed_prior=sensor_model.speedPrior(speed_val,bin_center_speed-speed_bin_size/2,bin_center_speed+speed_bin_size/2);
            transition_matrices(i_r,i_c,previous_speed_idx,previous_dir_idx)=transition_matrices(i_r,i_c,previous_speed_idx,previous_dir_idx)+dir_prior*speed_prior;
        end
        
        distance_dir_vector=abs(rad2deg(circ_dist(deg2rad(dir_states),deg2rad(dir_val))));
        distance_speed_vector=abs(speed_val-speed_states);
        
        dir_idx=find(distance_dir_vector==min(distance_dir_vector));
        previous_dir_idx=dir_idx(1);
        
        speed_idx=find(distance_speed_vector==min(distance_speed_vector));
        previous_speed_idx=speed_idx(1);
        
    end
    
    
    %---------------------------------------
    % Prediction Stage
    %--------------------------------------    
    
    
    prediction_matrix=zeros(n_rows,n_cols);
    
    
    for k=1:n_rows
        for l=1:n_cols
            trans_mat=transition_matrices(:,:,k,l);
            trans_mat=squeeze(trans_mat);
            prediction_matrix = prediction_matrix + PX(k,l).*trans_mat;
        end
    end

    
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
   

    %---------------------------------------
    % Correction Stage
    %--------------------------------------    
   
   PX=zeros(n_rows,n_cols)*(1/(n_rows*n_cols));
   for k=1:length(idx_r)
       i_r=idx_r(k);
       i_c=idx_c(k);
       bin_center_dir=dir_mesh(i_r,i_c);
       bin_center_speed=speed_mesh(i_r,i_c);
       dir_prior=sensor_model.dirPrior(dir_val,wrapTo360(bin_center_dir-dir_bin_size/2),wrapTo360(bin_center_dir+dir_bin_size/2));
       speed_prior=sensor_model.speedPrior(speed_val,bin_center_speed-speed_bin_size/2,bin_center_speed+speed_bin_size/2);
       corrected_value=dir_prior*speed_prior*prediction_matrix(i_r,i_c);
       state_accummulator(i_r,i_c)=state_accummulator(i_r,i_c)+corrected_value;
       state_accummulator_sens(i_r,i_c)=state_accummulator_sens(i_r,i_c)+dir_prior*speed_prior;
       PX(i_r,i_c)=corrected_value;
   end
 
 if (sum(PX(:))==0)
     PX=(1/(n_rows*n_cols)).*ones(n_rows,n_cols);
 else
       PX=PX./sum(PX(:));
 end
   state_posteriors=state_accummulator/sum(state_accummulator(:));
   sens_posteriors=state_accummulator_sens/sum(state_accummulator_sens(:));
    
end





