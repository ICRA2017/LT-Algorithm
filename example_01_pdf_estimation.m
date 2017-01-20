%========================================================================
%   example_01_pdf_estimation
%   version 1.0 - January 18th, 2017
%   
%   This example estimates Probability Density Functions (PDFs) of
%   simulated airflow data. The pdf estimation is performed according to
%   the equations 1 to 3 in the paper "Probabilistic Air Flow Modelling 
%   Using Turbulent and Laminar Characteristics for Ground and Aerial 
%   Robots". 
%========================================================================


clear;
clc;
close all;



%------------------------------------------------------------
%       Joint PDF estimation parameters
%------------------------------------------------------------
dir_bin_size=10;        % wind direction bin size (degrees)
dir_std=4;              % wind direction standard deviation for the sensor model (equation 2)

speed_bin_size=0.025;   % speed bin size (m/s)
speed_std=0.08;          % wind speed standard deviation for the sensor model (equation 2)
speed_range=[0 0.5];    % min and max values of the wind speed to be considered in the modeled PDF

%------------------------------------------------------------
%       Simulated data parameters
%------------------------------------------------------------
% The generated test measurements simulate two different wind flow
% conditions, one laminar flow and one turbulent flow. In this example,
% both are simulated as random gaussan data, centered at the same wind
% direction/speed means. The turbulent flow however, is simulated with a
% bigger STD.

n_meas=1500;            % number of TOTAL measurements
k_turbulence=0.35;      % percentage of the measurements belong to the turbulent flow

lam_mu_dir=100;         % mean wind direction (degrees)of the laminar flow
lam_std_dir=20;            % direction std of the laminar flow

lam_mu_spd=0.25;        % mean wind speed (m/s)of the laminar flow
lam_std_spd=0.05;       % sdt wind speed (m/s) of the laminar flow

turb_mu_dir=lam_mu_dir; % turbulent flow mean speed
turb_std_dir=1400;  % turbulent flow std

turb_mu_spd=lam_mu_spd;        % mean wind speed (m/s)of the laminar flow
turb_std_spd=0.1;       % sdt wind speed (m/s) of the laminar flow



%-- Simulated data generation
n_laminar=n_meas*(1-k_turbulence);
n_laminar=round(n_laminar);
n_turbulent=n_meas-n_laminar;

dir_laminar=wrapTo360(normrnd(lam_mu_dir,lam_std_dir,[n_laminar 1]));
dir_turbulent=wrapTo360(normrnd(turb_mu_dir,turb_std_dir,[n_turbulent 1]));

spd_laminar=normrnd(lam_mu_spd,lam_std_spd,[n_laminar,1]);
spd_turbulent=normrnd(turb_mu_spd,turb_std_spd,[n_turbulent 1]);

dir_vector=[dir_laminar(:);dir_turbulent(:)];
spd_vector=[spd_laminar(:);spd_turbulent(:)];
wind_data=[dir_vector spd_vector];

%--- Runs the Bayesian Filter and the computation of the Joint PDF
% Matrix initialization
[dir_mesh,speed_mesh,dir_states,speed_states]=generateMatrices(dir_bin_size,speed_bin_size,speed_range);

% Sensor model initialization
sensor_model=sensorModel(dir_std,speed_std);


% ---> Computes the joint PDFs using a Bayesian filter
[state_posteriors,sens_posteriors,histogram_out]=jointPDFBayes(dir_mesh,...
    speed_mesh,dir_states,speed_states,sensor_model,dir_std,speed_std,dir_bin_size,...
    speed_bin_size,wind_data);

%-- Plots the results
figure('units','normalized','outerposition',[0 0 1 1]);
subplot(1,3,1);

% Raw data
plot(wind_data(:,1),wind_data(:,2),'ro','markerfacecolor','r','markersize',3);
xlim([0 359]);
ylim(speed_range);
set(gca,'FontSize',20);
set(gca,'FontWeight','bold')
grid on;
xlabel('Wind direction ^\circ');
ylabel('Wind speed (m/s)');
title('Raw data');
% Raw histogram
subplot(1,3,2);

h=pcolor(dir_mesh,speed_mesh,histogram_out);
set(h, 'EdgeColor', 'none');
set(gca,'FontSize',20);
set(gca,'FontWeight','bold')
grid on;
xlabel('Wind direction ^\circ');
ylabel('Wind speed (m/s)');
title('Simple Histogram');

% Filter output
subplot(1,3,3);

h=pcolor(dir_mesh,speed_mesh,state_posteriors);
set(h, 'EdgeColor', 'none');
set(gca,'FontSize',20);
set(gca,'FontWeight','bold')
grid on;
xlabel('Wind direction ^\circ');
ylabel('Wind speed (m/s)');
title('\Gamma(J|z_{1:t})');











