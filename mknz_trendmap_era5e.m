close all
clear
clc

%% ERA5 Evaporation + Evapotranspiration
csr_can = importdata("csr_can_025.mat");
ealcotime = (datetime(2002,1,15):calmonths:datetime(2016,12,15))';
mknz025 = importdata("mknz_extent_025_lakesin.txt");
canada025 = importdata("canada_extent_full_lakesin_025.txt");
mknzcan_qd = mknz025.*canada025;
nanval = max(mknzcan_qd,[],'all');
mknzcan_qd(mknzcan_qd==nanval)=nan;
mknzcan_qd = mknzcan_qd(:,~all(isnan(mknzcan_qd)));
mknzcan_qd = mknzcan_qd(~all(isnan(mknzcan_qd),2),:);
mknzcan_qd(mknzcan_qd<=0) = nan;
mknzcan_qd(mknzcan_qd>0) = 1;
% Linear trends
% ERA5 ET
era5_et = ncread("era5_e_tp_updated2024-04-05.nc",'e');
era5ettime = ncread("era5_e_tp_updated2024-04-05.nc",'time');
era5ettime = datetime(1900,1,1,era5ettime,0,0);era5ettime = era5ettime(1:end-2);
[y,m,~] = ymd(era5ettime);
monthdays = eomday(y,m);
era5_et = squeeze(era5_et(:,:,1:end-2));
era5_et = era5_et.*-1000;
era5_et = pagetranspose(bsxfun(@times,era5_et,reshape(monthdays,1,1,size(monthdays,1))));

era5_et_mean_ind = isbetween(era5ettime,datetime(2004,1,1),datetime(2009,12,31));
era5ettime_ealcoidx = isbetween(era5ettime,ealcotime(1),ealcotime(end));
era5_et_dn_y = datenum(era5ettime)./365;

% ERA5 ET Trend Map
era5_e_cs_map = nan(size(era5_et));
for i = 1:size(era5_et,1)
    for j = 1:size(era5_et,2)
        et_tmp = squeeze(era5_et(i,j,:));
        if isnan(et_tmp(1)) == false
            era5e_mean = mean(et_tmp(era5ettime_ealcoidx));
            era5_e_mr = et_tmp-era5e_mean;
            era5_e_cs = cumsum(era5_e_mr);
            era5_e_mean_ind_mean = mean(era5_e_cs(era5_et_mean_ind));
            era5_e_cs = era5_e_cs-era5_e_mean_ind_mean;
            era5_e_cs_map(i,j,:) = era5_e_cs;
        end
    end
end

% ERA5 ET lat/lon
lon_era5 = double(ncread("era5_e_tp_updated2024-04-05.nc",'longitude'));
lat_era5 = double(ncread("era5_e_tp_updated2024-04-05.nc",'latitude'));
[lon_era5,lat_era5] = meshgrid(lon_era5,lat_era5);
% CSR lat/lon
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
% Create a CSR mascon mask
unique_csr = unique(csr_can(:,:,1));
csr_mascon_mask = csr_can(:,:,1);
for i = 1:size(unique_csr)
    unique_val = unique_csr(i);
    csr_mascon_mask(csr_mascon_mask==unique_val) = i;
end
% Create an ERA5 ET mask with same resolution as CSR
era5_et_csr_mask = nan(size(lon_era5));
era5_et_csr_mask_mknz = nan(size(lon_era5));
for i = 1:size(lon_era5,1)
    for j = 1:size(lon_era5,2)
        if isnan(era5_et(i,j)) == false
            lat_era5_tmp = lat_era5(i,j);
            lon_era5_tmp = lon_era5(i,j);
            
            [~,ind_lat] = min(abs(canYQ_025 - lat_era5_tmp),[],1);
            [~,ind_lon] = min(abs(canXQ_025 - lon_era5_tmp),[],2);
    
            era5_et_csr_mask(i,j) = csr_mascon_mask(ind_lat(1),ind_lon(1));
            if mknzcan_qd(ind_lat(1),ind_lon(1)) == 1
                era5_et_csr_mask_mknz(i,j) = csr_mascon_mask(ind_lat(1),ind_lon(1));
            end
        end
    end
end

% Apply Mackenzie mask to data
unq_era5et_csr_mknz = unique(era5_et_csr_mask_mknz(isnan(era5_et_csr_mask_mknz)==false));

era5_trend_csr_mknz = nan(size(era5_et,[1 2])); % Initialize trend map
era5_trend_csr_mknz_2016 = nan(size(era5_et,[1 2])); % Initialize trend map
era5_trend_csrmknz_pval = nan(size(era5_et,[1 2])); % Initialize trend map
era5_trend_csrmknz_pval_2016 = nan(size(era5_et,[1 2])); % Initialize trend map

% Initialize variables for trend removal
t = era5_et_dn_y;
Fc_1yr = 2*pi;
a1 = ones(size(t));
a2 = t;
a3 = sin(t*Fc_1yr);
a4 = cos(t*Fc_1yr);
A = [a1 a2 a3 a4];
Fc_6mo = 2*pi/0.5;
b1 = ones(size(t));
b2 = t;
b3 = sin(t*Fc_6mo);
b4 = cos(t*Fc_6mo);
B = [b1 b2 b3 b4];

for i = 1:length(unq_era5et_csr_mknz)
    unq_val = unq_era5et_csr_mknz(i);
    [row,col] = find(era5_et_csr_mask_mknz==unq_val);
    av_et_cs = squeeze(mean(era5_e_cs_map(row,col,:),[1 2]));

    % Least Squares
    A_soln = (A'*A)\(A'*av_et_cs);
    sin_A = A_soln(3)*sin(t*Fc_1yr)+A_soln(4)*cos(t*Fc_1yr);
   
    B_soln = (B'*B)\(B'*av_et_cs);    
    sin_B = B_soln(3)*sin(t*Fc_6mo)+B_soln(4)*cos(t*Fc_6mo);
    era5e_ann_6moSigRem = av_et_cs - sin_A - sin_B;

    [e_trend,~,~,~,stats] = regress(era5e_ann_6moSigRem,[era5_et_dn_y ones(size(era5_et_dn_y))]);
    era5_trend_csr_mknz(row,col) = e_trend(1);
    era5_trend_csrmknz_pval(row,col) = stats(3);
    [e_trend_2016,~,~,~,stats] = regress(era5e_ann_6moSigRem(era5ettime_ealcoidx),...
        [era5_et_dn_y(era5ettime_ealcoidx) ones(size(era5_et_dn_y(era5ettime_ealcoidx)))]);
    era5_trend_csr_mknz_2016(row,col) = e_trend_2016(1);
    era5_trend_csrmknz_pval_2016(row,col) = stats(3);
end

%% ERA5 Trend Mackenzie Basin Map
cmap = importdata("cmap_blueyellowred_symmetrical.mat");
subbasins = readgeotable("mknz_outline_new_08_2023.shp");
lakes = readgeotable('worldlakes.shp');
rivers = readgeotable('worldrivers.shp');

[row,col] = find(mknzcan_qd==1);
lat_025 = flip(42.125:0.25:83.125);
latlim = [min(lat_025) max(lat_025)];
lon_025 = -141:0.25:-53;
lonlim = [min(lon_025) max(lon_025)];
mknz_lon_range = sort([lon_025(min(col)) lon_025(max(col))]);
mknz_lat_range = sort([lat_025(min(row)) lat_025(max(row))]);
mknz_lat_range = [mknz_lat_range(1),mknz_lat_range(2)];
plot_latlim = [mknz_lat_range(1)-0.25 mknz_lat_range(2)+0.75];
plot_lonlim = [mknz_lon_range(1)-0.25 mknz_lon_range(2)+1.25];

%% Remove insignificant trends
era5_trend_csrmknz_pval(era5_trend_csrmknz_pval>=0.05) = 1;
era5_trend_csrmknz_pval(era5_trend_csrmknz_pval<0.05) = 0;
era5_trend_csrmknz_pval_2016(era5_trend_csrmknz_pval_2016>=0.05) = 1;
era5_trend_csrmknz_pval_2016(era5_trend_csrmknz_pval_2016<0.05) = 0;

unsig_ind = find(era5_trend_csrmknz_pval==1);

era5_trend_csr_mknz(unsig_ind) = nan;
figure
worldmap(plot_latlim,plot_lonlim)
setm(gca,'FLineWidth',1)
setm(gca,'FontSize',15)
pcolorm(lat_era5+0.0625,lon_era5+0.0625,era5_trend_csr_mknz*0.1) % mm to cm
clim([-2.5 2.5])
colormap(flipud(cmap))
geoshow(lakes,"FaceColor",'white','FaceAlpha',.1,'EdgeColor','black')
geoshow(subbasins,"FaceColor",'white','FaceAlpha',0.01,'EdgeColor','black')
geoshow(rivers,"Color","black")
bordersm('Canada')

unsig_ind = find(era5_trend_csrmknz_pval_2016==1);

era5_trend_csr_mknz_2016(unsig_ind) = nan;
figure
worldmap(plot_latlim,plot_lonlim)
setm(gca,'FLineWidth',1)
setm(gca,'FontSize',15)
pcolorm(lat_era5+0.0625,lon_era5+0.0625,era5_trend_csr_mknz_2016*0.1) % mm to cm
clim([-2.5 2.5])
colormap(flipud(cmap))
geoshow(lakes,"FaceColor",'white','FaceAlpha',.1,'EdgeColor','black')
geoshow(subbasins,"FaceColor",'white','FaceAlpha',0.01,'EdgeColor','black')
geoshow(rivers,"Color","black")
bordersm('Canada')
