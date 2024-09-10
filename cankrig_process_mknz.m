%% Plotting precipitation trends from NRCan stations
close all
clear
clc

%% Import data
time = (datetime(2002,1,31):calmonths:datetime(2019,09,30))';
[y,m,d] = ymd(time);
cankrig_dn_y = datenum(time)/365;
cankrig = importdata("cankrig_2002on.mat");
cankrig = cankrig{:,:};
cankrig(cankrig<0) = nan;
cankrig = cankrig';
statdetails = importdata("cankrig_latlon.mat");
subbasins = readgeotable("mknz_outline_new_08_2023.shp");
lakes = readgeotable('worldlakes.shp');
rivers = readgeotable('worldrivers.shp');

canada025 = importdata("canada_extent_full_lakesin_025.txt");
canada025(canada025~=-9999)=1;
canada025(canada025==-9999)=0;

mknz025 = importdata("mknz_extent_025_lakesin.txt");
canada025 = importdata("canada_extent_full_lakesin_025.txt");
mknzcan_qd = mknz025.*canada025;
nanval = max(mknzcan_qd,[],'all');
mknzcan_qd(mknzcan_qd==nanval)=nan;
mknzcan_qd = mknzcan_qd(:,~all(isnan(mknzcan_qd)));
mknzcan_qd = mknzcan_qd(~all(isnan(mknzcan_qd),2),:);
mknzcan_qd(mknzcan_qd<=0) = nan;
mknzcan_qd(mknzcan_qd>0) = 1;


%% Station Data
canada025(canada025>=0) = 1;
canada025(canada025<0) = nan;
canada = canada025(:,~all(isnan(canada025)));
canada = canada(~all(isnan(canada),2),:);

[row,col] = find(canada==1);
lat_025 = flip(42.125:0.25:83.125);
lon_025 = -141:0.25:-53;
can_lon_range = sort([lon_025(min(col)) lon_025(max(col))]);
can_lat_range = sort([lat_025(min(row)) lat_025(max(row))]);
can_lat_range = [can_lat_range(1),can_lat_range(2)];
xq = can_lon_range(1):0.25:can_lon_range(2);
yq = can_lat_range(1):0.25:can_lat_range(2);
[canXQ_025,canYQ_025] = meshgrid(xq,yq);
canYQ_025 = flip(canYQ_025);

for i = 1:size(statdetails,1)
    stat_lat = statdetails{i,1};
    stat_lon = statdetails{i,2};
    [~,ind_lat] = min(abs(canYQ_025 - stat_lat),[],1);
    [~,ind_lon] = min(abs(canXQ_025 - stat_lon),[],2);

    if isnan(mknzcan_qd(ind_lat,ind_lon))==true
        stat_sb(i,1) = 0;
    elseif mknzcan_qd(ind_lat,ind_lon) == 1
        stat_sb(i,1) = 1;
    end
end

ind_mknz = find(stat_sb==1);
cankrig_mknz = cankrig(:,ind_mknz);
cankrig_mknz_av = squeeze(mean(cankrig_mknz,2,'omitnan'));
