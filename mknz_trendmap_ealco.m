%% Plotting trends in the MRB
close all
clear
clc

csr_can = importdata("csr_can_025.mat");
ealcotime = (datetime(2002,1,15):calmonths:datetime(2016,12,15))';

subbasins = readgeotable("mknz_outline_new_08_2023.shp");
lakes = readgeotable('worldlakes.shp');
rivers = readgeotable('worldrivers.shp');

mknz025 = importdata("mknz_extent_025_lakesin.txt");
canada025 = importdata("canada_extent_full_lakesin_025.txt");
mknzcan_qd = mknz025.*canada025;
nanval = max(mknzcan_qd,[],'all');
mknzcan_qd(mknzcan_qd==nanval)=nan;
mknzcan_qd = mknzcan_qd(:,~all(isnan(mknzcan_qd)));
mknzcan_qd = mknzcan_qd(~all(isnan(mknzcan_qd),2),:);
mknzcan_qd(mknzcan_qd<=0) = nan;
mknzcan_qd(mknzcan_qd>0) = 1;
[row,col] = find(mknzcan_qd==1);
lat_025 = flip(42.125:0.25:83.125);
lon_025 = -141:0.25:-53;
mknz_lon_range = sort([lon_025(min(col)) lon_025(max(col))]);
mknz_lat_range = sort([lat_025(min(row)) lat_025(max(row))]);
mknz_lat_range = [mknz_lat_range(1),mknz_lat_range(2)];
plot_latlim = [mknz_lat_range(1)-0.25 mknz_lat_range(2)+0.75];
plot_lonlim = [mknz_lon_range(1)-0.25 mknz_lon_range(2)+1.25];


%% EALCO Evaporation + Evapotranspiration
ealco_et = importdata("ealco_et_wf.mat"); % mm
ealco_e0 = importdata("ealco_e0_wf.mat"); % mm
ealco_e_mean_ind = isbetween(ealcotime,datetime(2004,1,1),datetime(2009,12,31));
ealco_dn_y = datenum(ealcotime)./365;

%% EALCO ET lat/lon
coords = importdata("Canada_ET/coordinates_ealco.txt");
lon_ealco_et = coords.data(:,4);
lat_ealco_et = coords.data(:,5);
lon_ealco_et = reshape(lon_ealco_et,[1140,960])';
lat_ealco_et = reshape(lat_ealco_et,[1140,960])';
% EALCO E0 lat/lon
xind = 1:2:size(lon_ealco_et,1)-1;
yind = 1:2:size(lon_ealco_et,2)-1;
lon_ealco_e0 = lon_ealco_et(xind,yind);
lat_ealco_e0 = lat_ealco_et(xind,yind);
% CSR lat/lon
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
% Create a CSR mascon mask
unique_csr = unique(csr_can(:,:,1));
csr_mascon_mask = csr_can(:,:,1);
for i = 1:size(unique_csr)
    unique_val = unique_csr(i);
    csr_mascon_mask(csr_mascon_mask==unique_val) = i;
end
% Create an EALCO ET mask with same resolution as CSR
ealco_et_csr_mask = nan(size(lon_ealco_et));
ealco_et_csr_mask_mknz = nan(size(lon_ealco_et));
for i = 1:size(lon_ealco_et,1)
    for j = 1:size(lon_ealco_et,2)
        if isnan(ealco_et(i,j)) == false
            lat_ealco_et_tmp = lat_ealco_et(i,j);
            lon_ealco_et_tmp = lon_ealco_et(i,j);
            
            [~,ind_lat] = min(abs(canYQ_025 - lat_ealco_et_tmp),[],1);
            [~,ind_lon] = min(abs(canXQ_025 - lon_ealco_et_tmp),[],2);
    
            ealco_et_csr_mask(i,j) = csr_mascon_mask(ind_lat(1),ind_lon(1));
            if mknzcan_qd(ind_lat(1),ind_lon(1)) == 1
                ealco_et_csr_mask_mknz(i,j) = csr_mascon_mask(ind_lat(1),ind_lon(1));
            end
        end
    end
end
% Create an EALCO E0 mask with same resolution as CSR
ealco_e0_csr_mask = nan(size(lon_ealco_e0));
ealco_e0_csr_mask_mknz = nan(size(lon_ealco_e0));
for i = 1:size(lon_ealco_e0,1)
    for j = 1:size(lon_ealco_e0,2)
        if isnan(ealco_e0(i,j)) == false
            lat_ealco_e0_tmp = lat_ealco_e0(i,j);
            lon_ealco_e0_tmp = lon_ealco_e0(i,j);
            
            [~,ind_lat] = min(abs(canYQ_025 - lat_ealco_e0_tmp),[],1);
            [~,ind_lon] = min(abs(canXQ_025 - lon_ealco_e0_tmp),[],2);
    
            ealco_e0_csr_mask(i,j) = csr_mascon_mask(ind_lat(1),ind_lon(1));
            if mknzcan_qd(ind_lat(1),ind_lon(1)) == 1
                ealco_e0_csr_mask_mknz(i,j) = csr_mascon_mask(ind_lat(1),ind_lon(1));
            end
        end
    end
end
%% Apply mask to data, and add ET to E0
unq_ealcoet_csr = unique(ealco_et_csr_mask(isnan(ealco_et_csr_mask)==false));
unq_ealcoe0_csr = unique(ealco_e0_csr_mask(isnan(ealco_e0_csr_mask)==false));
if length(unq_ealcoe0_csr)~=length(unq_ealcoet_csr)
    error('The number of unique values in EALCO ET and EALCO E0 are different')
end

%% Mackenzie Basin only
unq_ealcoet_csr_mknz = unique(ealco_et_csr_mask_mknz(isnan(ealco_et_csr_mask_mknz)==false));
unq_ealcoe0_csr_mknz = unique(ealco_e0_csr_mask_mknz(isnan(ealco_e0_csr_mask_mknz)==false));
if length(unq_ealcoe0_csr_mknz)~=length(unq_ealcoe0_csr_mknz)
    error('The number of unique values in EALCO ET and EALCO E0 are different')
end
ealco_e_trendmap_mknz = nan(size(ealco_e0,[1 2])); % Initialize trend map
ealco_trend_csr_pval_mknz = nan(size(ealco_e0,[1 2])); % Initialize trend map

% Initialize variables for removing harmonic trends
t = ealco_dn_y;
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

for i = 1:length(unq_ealcoet_csr_mknz)
% Add E0 and ET to get total E
    unq_val = unq_ealcoet_csr_mknz(i);
    [row_et,col_et] = find(ealco_et_csr_mask_mknz==unq_val);
    ind_et = find(ealco_et_csr_mask_mknz==unq_val);
    [row_e0,col_e0] = find(ealco_e0_csr_mask_mknz==unq_val);
    ind_e0 = find(ealco_e0_csr_mask_mknz==unq_val);
    for rc = 1:length(row_et)
        av_et_all(rc,:) = squeeze(ealco_et(row_et(rc),col_et(rc),:));
    end
    av_et = squeeze(mean(av_et_all,1))';
    for rc = 1:length(row_e0)
        av_e0_all(rc,:) = squeeze(ealco_e0(row_e0(rc),col_e0(rc),:));
    end
    av_e0 = squeeze(mean(av_e0_all,1))';
    av_e = av_et + av_e0;
    % Compute cumulative sum of E
    ealco_e_mean = mean(av_e);
    ealco_e_mr = av_e-ealco_e_mean;
    ealco_e_cs = cumsum(ealco_e_mr);
    ealco_e_mean_ind_mean = mean(ealco_e_cs(ealco_e_mean_ind));
    ealco_e_cs = ealco_e_cs-ealco_e_mean_ind_mean;

    % Least Squares - Annual
    A_soln = (A'*A)\(A'*ealco_e_cs);   
    sin_A = A_soln(3)*sin(t*Fc_1yr)+A_soln(4)*cos(t*Fc_1yr);    
    B_soln = (B'*B)\(B'*ealco_e_cs);
    sin_B = B_soln(3)*sin(t*Fc_6mo)+B_soln(4)*cos(t*Fc_6mo);
    ealco_ann_6moSigRem = ealco_e_cs - sin_A - sin_B;

    [trend,~,~,~,stats] = regress(ealco_ann_6moSigRem,[ealco_dn_y ones(size(ealco_dn_y))]);
    ealco_e_trendmap_mknz(ind_e0) = trend(1);
    ealco_trend_csr_pval_mknz(ind_e0) = stats(3);
end

%% EALCO Trend Mackenzie Basin Map
cmap = importdata("cmap_blueyellowred_symmetrical.mat");

% Remove insignificant trends
ealco_trend_csr_pval_mknz(ealco_trend_csr_pval_mknz>=0.05) = 1;
ealco_trend_csr_pval_mknz(ealco_trend_csr_pval_mknz<0.05) = 0;

unsig_ind = find(ealco_trend_csr_pval_mknz==1);

ealco_e_trendmap_mknz(unsig_ind) = nan;
figure
worldmap([plot_latlim(1)+0.25,plot_latlim(2)+0.625],plot_lonlim)
setm(gca,'FLineWidth',1)
setm(gca,'FontSize',15)
pcolorm(lat_ealco_e0,lon_ealco_e0,ealco_e_trendmap_mknz*0.1) % mm to cm
clim([-2.5 2.5])
colormap(flipud(cmap))
geoshow(lakes,"FaceColor",'white','FaceAlpha',.1,'EdgeColor','black')
geoshow(subbasins,"FaceColor",'white','FaceAlpha',0.01,'EdgeColor','black')
geoshow(rivers,"Color","black")
bordersm('Canada')
