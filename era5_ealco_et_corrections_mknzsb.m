%% Extend EALCO using ERA5 Precipitation
close all
clear
clc

% Import
ealcotime = (datetime(2002,1,15):calmonths:datetime(2016,12,15))';
era5ettime = ncread("era5_e_tp_updated2024-04-05.nc",'time');
era5ettime = datetime(1900,1,1,era5ettime,0,0);era5ettime = era5ettime(1:end-2);
era5_mknz = -1*importdata("era5_te_mknz_av.mat");
ealco_mknz = importdata("ealco_et_wf_mknz_av.mat");
ealco_e0_mknz = importdata("ealco_e0_wf_mknz_av.mat");

ealco_mknz = ealco_mknz + ealco_e0_mknz;
diff_mean_mknz = mean(era5_mknz(1:size(ealcotime,1)))/mean(ealco_mknz);

% Apply correction to end of ERA5
era5ext_ind = find(era5ettime>ealcotime(end));
era5ext_time = era5ettime(era5ext_ind);
era5ext_mknz = era5_mknz(era5ext_ind)/diff_mean_mknz;

%% Plot the results
figure
plot(ealcotime,ealco_mknz);
hold on
plot(era5ettime,era5_mknz);
hold on
plot(era5ext_time,era5ext_mknz)
legend('EALCO','ERA5','ERA5 Adjusted')
title('mknz')

