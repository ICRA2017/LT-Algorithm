%========================================================================
%   samplePdf
%   version 1.0 - January 18th, 2017
%
%   Acquires artificial samples from a joint PDF
%
%   inputs:
%   dir_states, speed_states: vectors that contain the possible wind speed/wind  
%   dir_bin_size,speed_bin_size: Discretization parameters for the state
%   space.
%   pdf_in: joint pdf to be sampled from
%   n_meas: number of samples
%   outputs:
%   dir_vector,speed_vector: artificial samples
%========================================================================


function [dir_vector,speed_vector]=samplePdf(dir_states,speed_states,dir_bin_size,speed_bin_size,pdf_in,n_meas)
dir_vector=[];
speed_vector=[];
[r_idx,c_idx]=find(pdf_in~=0);
meas_per_state=length(c_idx)*100;
for i=1:length(r_idx)
    posterior=pdf_in(r_idx(i),c_idx(i));
    n_samp=round(posterior*meas_per_state);
    curr_dir=dir_states(c_idx(i));
    curr_spd=speed_states(r_idx(i));
    spd_values=curr_spd+speed_bin_size.*rand(n_samp,1);
    dir_values=curr_dir+dir_bin_size.*rand(n_samp,1);
    dir_vector=[dir_vector;dir_values];
    speed_vector=[speed_vector;spd_values];
end
dir_vector=normalizeAngles(dir_vector);
rand_idx=randi(length(dir_vector),[n_meas,1]);
dir_vector=dir_vector(rand_idx);
speed_vector=speed_vector(rand_idx);

end