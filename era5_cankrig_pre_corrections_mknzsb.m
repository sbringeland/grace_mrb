%% Extend CanKrig using ERA5 Precipitation
close all
clear
clc

% Import
cankrigtime = (datetime(2002,1,15):calmonths:datetime(2019,09,15))';
era5tptime = ncread("era5_e_tp_updated2024-04-05.nc",'time');
era5tptime = datetime(1900,1,1,era5tptime,0,0);era5tptime = era5tptime(1:end-2);
era5_mknz = importdata("era5_tp_mknz_av.mat");
cankrig_mknz = importdata("cankrig_mknz_av.mat");

diff_mean_mknz = mean(era5_mknz(1:size(cankrigtime,1)))/mean(cankrig_mknz);

% Apply correction to end of ERA5
era5ext_ind = find(era5tptime>cankrigtime(end));
era5ext_time = era5tptime(era5ext_ind);
era5ext_mknz = era5_mknz(era5ext_ind)/diff_mean_mknz;

%% Plot the results
figure
plot(cankrigtime,cankrig_mknz);
hold on
plot(era5tptime,era5_mknz);
hold on
plot(era5ext_time,era5ext_mknz)
legend('CanKrig','ERA5','ERA5 Adjusted')
title('mknz')





