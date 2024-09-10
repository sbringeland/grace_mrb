%% ERA 5 Precipitation
close all
clear
clc 
%% Import
time = ncread("era5_e_tp_updated2024-04-05.nc",'time');
time = datetime(1900,1,1,time,0,0);time = time(1:end-2);
[y,m,d] = ymd(time);
monthdays = eomday(y,m);
era5_dn_y = datenum(time)/365;
tp025deg = importdata("era5_pre_025deg_04_2024.mat");
tp025deg = pagetranspose(tp025deg);tp025deg = tp025deg(4:end-5,4:end-5,1:end-2);
% Transform to mm/month
tp025deg=tp025deg.*1000;tp025deg = bsxfun(@times,tp025deg,reshape(monthdays,1,1,size(monthdays,1)));
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


for j = 1:size(tp025deg,3)
    era5_page = tp025deg(:,:,j);
    tp_mknz_page = era5_page.*mknzcan_qd;
    tp_mknz_page = tp_mknz_page(:,~all(isnan(tp_mknz_page)));
    tp_mknz_page = tp_mknz_page(~all(isnan(tp_mknz_page),2),:);
    tp_mknz(:,:,j) = tp_mknz_page;
end

tp_mknz_av = squeeze(mean(tp_mknz,[1,2],'omitnan'));
