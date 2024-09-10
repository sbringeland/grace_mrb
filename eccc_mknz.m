%% Flow Rate for Mackenzie River Basin (& subbasins)
close all
clear 
clc

% Import
mknz_flow_mean_all = importdata("10KA001_Monthly_MeanFlow_ts.csv");
mknz_flow_mean = mknz_flow_mean_all.data;

%% Separate date and flow
mknz_dt = string(mknz_flow_mean_all.textdata(2:end-2,3));
mknz_dt = datetime(mknz_dt,'InputFormat',"MM--uuuu");mknz_dt.Day = 15;

%% Remove seasonal average
mknz_ind = find(mknz_dt>=datetime(2002,01,01));
mknz_dt = mknz_dt(mknz_ind);
mknz_flow_mean = mknz_flow_mean(mknz_ind);
[y,m,~] = ymd(mknz_dt);
eomday_mknz = eomday(y,m);
mknz_dn_y = datenum(mknz_dt)/365;
mknz_flow_mean = 3600*24*mknz_flow_mean.*eomday_mknz;
