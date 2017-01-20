%========================================================================
%   example_02_airflow_modelings
%   version 1.0 - January 18th, 2017
%   
%  In this example, an airflow map is created from sample measurements
%  collected with a quadcopter. The airflow map is computed using the
%  implementation of the LT Algorithm presented in the paper  "Probabilistic 
%   Air Flow Modelling Using Turbulent and Laminar Characteristics for 
%   Ground and Aerial Robots". 
%========================================================================

clear;
clc;
close all;
%------------------------------------------------------------
%   Algorithm Parameters
%------------------------------------------------------------

%----  Joint PDF computation parameters (Bayesian Filter) Equations 1 to 3
dir_bin_size=10;        % wind direction bin size (degrees)
dir_std=14;             % wind direction standard deviation for the sensor model (equation 2)


speed_bin_size=0.25;    % speed bin size (m/s)
speed_std=0.36;         % wind speed standard deviation for the sensor model (equation 2)
speed_range=[0 4];      % min and max values of the wind speed to be considered in the modeled PDF


%----  LT Algorithm Parameters
kernel_size=0.5;        % size (in m) of the Gaussian kernel for extrapolation (equations 4d) 
cell_size=0.5;          % To create the maps, the explored area is discretized in a grid of cells of size equal to cell_size (in m)
n_samples=50;           % Number of artificial samples used in the extrapolation process (Equation 5b) 

% Parameters for the computation for the turbulence indicator. In this
% particular implementation, the turbulence indicator (equation 5c) is
% computed using a logistic function
slope=10;               % growth rate (kappa in equation 5c)
turn_point=0.5;         % mid point (set as 0.5 in the paper)



%------------------------------------------------------------
%   Initialization
%------------------------------------------------------------
% Matrix initialization
[dir_mesh,speed_mesh,dir_states,speed_states]=generateMatrices(dir_bin_size,speed_bin_size,speed_range);
% Sensor model initialization
sensor_model=sensorModel(dir_std,speed_std);


%------------------------------------------------------------
%   Data loading
%------------------------------------------------------------
load('sample_data');

dir_vector=wrapTo360(wind_direction);
speed_vector=wind_speed;
position_index_vector=way_point_vector;
position_x_vector=command_x;
position_y_vector=command_y;

% minimum and maximum values in x,y are defined in the explored area
min_x=min(position_x_vector)-cell_size/2;
max_x=max(position_x_vector)+cell_size/2;
min_y=min(position_y_vector)-cell_size/2;
max_y=max(position_y_vector)+cell_size/2;



% 
% 
% %------------------------------------------------------------
% %------------------------------------------------------------
% % C) Generates the PDFs
% %------------------------------------------------------------
% %------------------------------------------------------------
% 
% 
% 
[n_rows,n_cols]=size(dir_mesh);
n_locations=max(position_index_vector);
pos_vector_x=zeros(n_locations,1);
pos_vector_y=zeros(n_locations,1);
pdf_vector=zeros(n_rows,n_cols,n_locations);



parfor selected_location=1:n_locations
    
    selected_idx=position_index_vector==selected_location;
    current_wind_speed=speed_vector(selected_idx);
    current_wind_speed=current_wind_speed(end:-1:1);
    current_wind_dir=dir_vector(selected_idx);
    current_wind_dir=current_wind_dir(end:-1:1);
    
    wind_data=[current_wind_dir(:) current_wind_speed(:)];
    
    
   
    
    [state_posteriors,sens_posteriorsx,histogram_outx]=jointPDFBayes(dir_mesh,speed_mesh,dir_states,speed_states,sensor_model,dir_std,speed_std,dir_bin_size,speed_bin_size,wind_data);
     pdf_vector(:,:,selected_location)=state_posteriors;

  
    pos_vector_x(selected_location)=mean(position_x_vector(selected_idx));
    pos_vector_y(selected_location)=mean(position_y_vector(selected_idx));
    
    fprintf('\n%d.- Joint PDF computed at location: %f, %f',selected_location,mean(position_x_vector(selected_idx)),mean(position_y_vector(selected_idx)));
end

% 

%------------------------------------------------------------
%------------------------------------------------------------
% E) Extrapolates the PDFs
%------------------------------------------------------------
%------------------------------------------------------------

xv=linspace(min_x,max_x,round( (max_x-min_x)/cell_size ));
yv=linspace(min_y,max_y,round( (max_y-min_y)/cell_size ));


[x_mesh,y_mesh]=meshgrid(xv,yv);

[nr,nc]=size(x_mesh);

pdf_vector_map=zeros(n_rows,n_cols,nr,nc);
turbulence_map=zeros(nr,nc);
expected_dir_map=zeros(nr,nc);
expected_speed_map=zeros(nr,nc);
x_loc_map=zeros(nr*nc,1);
y_loc_map=zeros(nr*nc,1);
pdf_locations=[pos_vector_x pos_vector_y];


fprintf('\n Running the LT Algorithm on query locations');
parfor i=1:nr
    fprintf('\n');
    for j=1:nc
        current_x=x_mesh(i,j);
        current_y=y_mesh(i,j);
        query_location=[current_x current_y];
        [L_pdf]=laminarExtrapolation(dir_mesh,speed_mesh,dir_states,sensor_model,speed_states,dir_std,speed_std,dir_bin_size,speed_bin_size,pdf_vector,pdf_locations,query_location,kernel_size,n_samples);
        [T_pdf]=turbulentExtrapolation(pdf_vector,pdf_locations,query_location,kernel_size);
        [LT_pdf,t_indicator]=LTExtrapolation(L_pdf,T_pdf,dir_mesh,slope,turn_point);
        [expected_dir,expected_speed]=computeExpectedValue(dir_mesh,speed_mesh,LT_pdf);
        pdf_vector_map(:,:,i,j)=LT_pdf;
        expected_dir_map(i,j)=expected_dir;
        expected_speed_map(i,j)=expected_speed;
        turbulence_map(i,j)=t_indicator;
        fprintf('.');
    end
end

% 
% ------------------------------------------------------------
% ------------------------------------------------------------
%   Plots the wind map
% ------------------------------------------------------------
% ------------------------------------------------------------
uMat=expected_speed_map.*cosd(expected_dir_map);
vMat=expected_speed_map.*sind(expected_dir_map);
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
plotWindMap(x_mesh,y_mesh,uMat,vMat,'jet','k',[1 1]);
plot(pos_vector_x,pos_vector_y,'wo','markersize',12,'linewidth',3);
axis equal;
%title('Expected map');
set(gca,'FontSize',25);
set(gca,'FontWeight','bold');
ylabel('Y(m)');
xlabel('X(m)');
caxis([0 2.3]);
title('Airflow map');


% 
% ------------------------------------------------------------
% ------------------------------------------------------------
%   Plots the turbulence map
% ------------------------------------------------------------
% ------------------------------------------------------------

figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
[hC hC]=contourf(x_mesh,y_mesh,turbulence_map);
set(hC,'LineStyle','none');
plot(pos_vector_x,pos_vector_y,'wo','markersize',12,'linewidth',3);
axis equal;
%title('Turbulence map');
colormap hot;
colorbar;
caxis([0 1]);
set(gca,'FontSize',25);
set(gca,'FontWeight','bold');
ylabel('Y(m)');
xlabel('X(m)');
title('Turbulence map')

