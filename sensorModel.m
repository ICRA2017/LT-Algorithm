%========================================================================
%   sensorModel
%   version 1.0 - January 18th, 2017
%   
%   According to what is stated in the paper "Probabilistic Air Flow Modelling 
%   Using Turbulent and Laminar Characteristics for Ground and Aerial 
%   Robots", this script creates an object of the type "sensorModel". This
%   object is used in the computation of the joint PDFs. "sensorModel"
%   implements a 2D gaussian sensor model according to equation 2 in the
%   paper. 
%   The sensor model requires two initialization parameters as follows
%   dirSTD - wind direction standard deviation
%   spdSTD - wind speed standard deviation
%========================================================================



classdef sensorModel < handle
    properties
       
        dirSTD;
        dirPDF;
        dirSupport;
        
        spdSTD
        spdPDF
        spdSupport;

        max_spd=60;
    end
    
    methods
         function obj=sensorModel(dir_std,spd_std)
             
             
             % Direction model
     
             obj.dirSTD=normalizeAngles(dir_std);
             
             
             delta_dir=0.1;
             obj.dirSupport=0:delta_dir:360;
             
             dirSupport_rad=deg2rad(obj.dirSupport);
             
             std_in=normalizeAngles(obj.dirSTD);
             std_in_rad=deg2rad(std_in);
             mu_in=normalizeAngles(180);
             mu_rad=deg2rad(mu_in);
             
             norm_factor=1/(std_in_rad*sqrt(2*pi));
             K=-1:1;
             
             exp_factor=zeros(length(K),1);
             obj.dirPDF=zeros(length(obj.dirSupport),1);
             
             for i=1:length(dirSupport_rad)
                 curr_suppport_rad=dirSupport_rad(i);
                 for j=1:length(K)
                     efactor= ( curr_suppport_rad-mu_rad+2*pi*K(j) ).^2 / (2*std_in_rad*std_in_rad);
                     exp_factor(j)=exp(-efactor);
                 end
                 obj.dirPDF(i)=norm_factor*sum(exp_factor);
                 
             end
             
             area=trapz(obj.dirPDF);
            obj.dirPDF=obj.dirPDF./area;
            
            % speed model
            
            delta_spd=0.05;
            obj.spdSupport=0:delta_spd:obj.max_spd;
            obj.spdSTD=spd_std;
            obj.spdPDF=normpdf(obj.spdSupport,obj.max_spd/2,spd_std);
            area=trapz(obj.spdPDF);
            obj.spdPDF=obj.spdPDF./area;
         end
         
         function [p_out]=speedPrior(obj,sensor_read,min_state,max_state)
             
             delta_spd=abs(obj.max_spd/2-sensor_read);
             
             
             
             min_state=min_state+delta_spd;
    
             
             max_state=max_state+delta_spd;

             
             [val_x,min_idx]=min(abs(obj.spdSupport-min_state));
             [val_x,max_idx]=min(abs(obj.spdSupport-max_state));
             
             area_region=sum(obj.spdPDF(min_idx:max_idx));
             p_out=area_region;
             
             
             
         end
         
         
         function [p_out]=dirPrior(obj,sensor_read,min_state,max_state)
             delta_loc=180-sensor_read;
             
             min_state=min_state+delta_loc;
             min_state=normalizeAngles(min_state);
             
             max_state=max_state+delta_loc;
             max_state=normalizeAngles(max_state);
             
             [val_x,min_idx]=min(abs(obj.dirSupport-min_state));
             [val_x,max_idx]=min(abs(obj.dirSupport-max_state));
             
             area_region=sum(obj.dirPDF(min_idx:max_idx));
             p_out=area_region; 
         end
    end
end

