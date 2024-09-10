%% Plotting evapotranspiration trends from EALCO and ERA5
close all
clear
clc

%% Import data
ealco_time = (datetime(2002,1,31):calmonths:datetime(2016,12,31))';
ealco_et = importdata("ealco_et_wf.mat");
coords_et = importdata("Canada_ET/coordinates_ealco.txt");
lon_list = coords_et.data(:,4);
lat_list = coords_et.data(:,5);
lon_map_et = reshape(lon_list,[1140,960])';
lat_map_et = reshape(lat_list,[1140,960])';

mknz025 = importdata("mknz_extent_025_lakesin.txt");
canada025 = importdata("canada_extent_full_lakesin_025.txt");
mknzcan_qd = mknz025.*canada025;
nanval = max(mknzcan_qd,[],'all');
mknzcan_qd(mknzcan_qd==nanval)=nan;
mknzcan_qd = mknzcan_qd(:,~all(isnan(mknzcan_qd)));
mknzcan_qd = mknzcan_qd(~all(isnan(mknzcan_qd),2),:);
mknzcan_qd(mknzcan_qd<=0) = nan;
mknzcan_qd(mknzcan_qd>0) = 1;
mknzcan_qd = mknzcan_qd(5:end,:);

%% EALCO Data
canada025(canada025>=0) = 1;
canada025(canada025<0) = nan;
canada = canada025(:,~all(isnan(canada025)));
canada = canada(~all(isnan(canada),2),:);

[row,col] = find(canada==1);
lat_025 = flip(42.125:0.25:83.125);
latlim = [min(lat_025) max(lat_025)];
lon_025 = -141:0.25:-53;
lonlim = [min(lon_025) max(lon_025)];
can_lon_range = sort([lon_025(min(col)) lon_025(max(col))]);
can_lat_range = sort([lat_025(min(row)) lat_025(max(row))]);
can_lat_range = [can_lat_range(1),can_lat_range(2)];
xq = can_lon_range(1):0.25:can_lon_range(2);
yq = can_lat_range(1):0.25:can_lat_range(2);
[canXQ_025,canYQ_025] = meshgrid(xq,yq);
canYQ_025 = flip(canYQ_025);
canXQ_025 = canXQ_025(5:end,:);
canYQ_025 = canYQ_025(5:end,:);

% ET
ealco_et_sb = nan(size(ealco_et,[1 2]));
for i = 1:size(ealco_et,1)
    for j = 1:size(ealco_et,2)
        ealco_lat = lat_map_et(i,j);
        ealco_lon = lon_map_et(i,j);
        [~,ind_lat] = min(abs(canYQ_025 - ealco_lat),[],1);
        [~,ind_lon] = min(abs(canXQ_025 - ealco_lon),[],2);
    
        if isnan(mknzcan_qd(ind_lat,ind_lon))==true
            ealco_et_sb(i,j) = 0;
        elseif mknzcan_qd(ind_lat,ind_lon) == 1
            ealco_et_sb(i,j) = 1;
        end
    end                                                                                                                                                                                                                                                                                                                                                    
end

[row_mknz,col_mknz] = find(ealco_et_sb==1);

% Mackenzie
for i = 1:size(row_mknz,1)
    ealco_mknz(:,i) = squeeze(ealco_et(row_mknz(i),col_mknz(i),:));
end

ealco_et_mknz_av = squeeze(mean(ealco_mknz,2,'omitnan'));
