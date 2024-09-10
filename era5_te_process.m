%% ERA 5 Evapotranspiration
close all
clear
clc 
%% Import
te = ncread("era5_e_tp_updated2024-04-05.nc",'e');
lon = ncread("era5_e_tp_updated2024-04-05.nc",'longitude');
lat = ncread("era5_e_tp_updated2024-04-05.nc",'latitude');
time = ncread("era5_e_tp_updated2024-04-05.nc",'time');
time = datetime(1900,1,1,time,0,0);time = time(1:end-2);
[y,m,d] = ymd(time);
monthdays = eomday(y,m);
era5_dn_y = datenum(time)/365;

%%
te025deg = importdata("era5_et_025deg_new042024.mat");
te025deg = pagetranspose(te025deg);te025deg = te025deg(4:end-5,4:end-5,1:end-2);
% Convert to mm/month
te025deg=te025deg.*1000;te025deg = bsxfun(@times,te025deg,reshape(monthdays,1,1,size(monthdays,1)));
%%
mknz025 = importdata("mknz_extent_025_lakesin.txt");
canada025 = importdata("canada_extent_full_lakesin_025.txt");
mknzcan_qd = mknz025.*canada025;
nanval = max(mknzcan_qd,[],'all');
mknzcan_qd(mknzcan_qd==nanval)=nan;
mknzcan_qd = mknzcan_qd(:,~all(isnan(mknzcan_qd)));
mknzcan_qd = mknzcan_qd(~all(isnan(mknzcan_qd),2),:);
mknzcan_qd(mknzcan_qd<=0) = nan;
mknzcan_qd(mknzcan_qd>0) = 1;
canada025(canada025~=-9999)=1;
canada025(canada025==-9999)=0;

subbasins = readgeotable("mknz_outline_new_08_2023.shp");
lakes = readgeotable('worldlakes.shp');
rivers = readgeotable('worldrivers.shp');

for j = 1:size(te025deg,3)
    era5_page = te025deg(:,:,j);
    te_mknz_page = era5_page.*mknzcan_qd;
    te_mknz_page = te_mknz_page(:,~all(isnan(te_mknz_page)));
    te_mknz_page = te_mknz_page(~all(isnan(te_mknz_page),2),:);
    te_mknz(:,:,j) = te_mknz_page;
end

te_mknz_av = squeeze(mean(te_mknz,[1,2],'omitnan'));


