close all
clear
clc

%% CanKrig Precipitation
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

cankrigtime = (datetime(2002,1,31):calmonths:datetime(2019,09,30))';
[y,m,d] = ymd(cankrigtime);
cankrig_dn_y = datenum(cankrigtime)/365;
cankrig_p = importdata("cankrig_2002on.mat");
cankrig_p = cankrig_p{:,:};
cankrig_p(cankrig_p<-1) = nan;
cankrig_p(cankrig_p<0) = 0;
cankrig_p = reshape(cankrig_p,[470,621,213]);

cankrig_p_mean_ind = isbetween(cankrigtime,datetime(2004,1,1),datetime(2009,12,31));
cankrigtime_ealcoidx = isbetween(cankrigtime,ealcotime(1),ealcotime(end));
cankrig_dn_y = datenum(cankrigtime)./365;

% CanKrig-P Map
cankrig_p_cs_2016 = nan(size(cankrig_p));
for i = 1:size(cankrig_p,1)
    for j = 1:size(cankrig_p,2)
        tp_tmp = squeeze(cankrig_p(i,j,:));
        if isnan(tp_tmp(1)) == false
            cankrig_p_mean = mean(tp_tmp(cankrigtime_ealcoidx));
            cankrig_p_mr = tp_tmp-cankrig_p_mean;
            cankrig_p_cs = cumsum(cankrig_p_mr);
            cankrig_p_mean_ind_mean = mean(cankrig_p_cs(cankrig_p_mean_ind));
            cankrig_p_cs = cankrig_p_cs-cankrig_p_mean_ind_mean;
            cankrig_p_cs_2016(i,j,:) = cankrig_p_cs;
        end
    end
end
% CanKrig lat/lon
cankrig_latlon = importdata("cankrig_latlon.mat");
lat_cankrig = cankrig_latlon{:,1};
lon_cankrig = cankrig_latlon{:,2};
lat_cankrig = reshape(lat_cankrig,[470,621]);
lon_cankrig = reshape(lon_cankrig,[470,621]);
lon_cankrig(lon_cankrig>0) = lon_cankrig(lon_cankrig>0)*-1;
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
% Create an mask with same resolution as CSR
cankrig_p_csr_mask = nan(size(lon_cankrig));
cankrig_p_csr_mask_mknz = nan(size(lon_cankrig));
for i = 1:size(lon_cankrig,1)
    for j = 1:size(lon_cankrig,2)
        if isnan(cankrig_p(i,j)) == false
            lat_cankrig_tmp = lat_cankrig(i,j);
            lon_cankrig_tmp = lon_cankrig(i,j);
            
            [~,ind_lat] = min(abs(canYQ_025 - lat_cankrig_tmp),[],1);
            [~,ind_lon] = min(abs(canXQ_025 - lon_cankrig_tmp),[],2);
    
            cankrig_p_csr_mask(i,j) = csr_mascon_mask(ind_lat(1),ind_lon(1));
            if mknzcan_qd(ind_lat(1),ind_lon(1)) == 1
                cankrig_p_csr_mask_mknz(i,j) = csr_mascon_mask(ind_lat(1),ind_lon(1));
            end
        end
    end
end

% Apply Mackenzie mask to data
unq_cankrig_csr_mknz = unique(cankrig_p_csr_mask_mknz(isnan(cankrig_p_csr_mask_mknz)==false));
cankrig_trend_csr_mknz = nan(size(cankrig_p,[1 2])); % Initialize trend map
cankrig_trend_csr_mknz_pval = nan(size(cankrig_p,[1 2])); % Initialize pval map

cankrig_trend_csr_mknz_2016 = nan(size(cankrig_p,[1 2])); % Initialize trend map
cankrig_trend_csr_mknz_pval_2016 = nan(size(cankrig_p,[1 2])); % Initialize pval map

% Initialize variables for sine trend removal
t = cankrig_dn_y;
Fc_1yr = 2*pi;
a1 = ones(size(t));
a2 = t;
a3 = sin(t*Fc_1yr);
a4 = cos(t*Fc_1yr);
A = [a1 a2 a3 a4];

Fc_5yr = 2*pi/5.4;
b1 = ones(size(t));
b2 = t;
b3 = sin(t*Fc_5yr);
b4 = cos(t*Fc_5yr);
B = [b1 b2 b3 b4];
for i = 1:length(unq_cankrig_csr_mknz)
    unq_val = unq_cankrig_csr_mknz(i);
    [row_p,col_p] = find(cankrig_p_csr_mask_mknz==unq_val);
    ind_p = find(cankrig_p_csr_mask_mknz==unq_val);
    for rc = 1:length(row_p)
        tp_cs_all(rc,:) = squeeze(cankrig_p_cs_2016(row_p(rc),col_p(rc),:));
    end
    av_p_cs = squeeze(mean(tp_cs_all,1))';

    % Least Squares - Annual
    A_soln = (A'*A)\(A'*av_p_cs);    
    sin_A = A_soln(3)*sin(t*Fc_1yr)+A_soln(4)*cos(t*Fc_1yr);
        
    % Least Squares - 5 years
    B_soln = (B'*B)\(B'*av_p_cs);
    sin_B = B_soln(3)*sin(t*Fc_5yr)+B_soln(4)*cos(t*Fc_5yr);
    cankrig_ann_5yrSigRem = av_p_cs - sin_A - sin_B;

    [trend,~,~,~,stats] = regress(cankrig_ann_5yrSigRem,[cankrig_dn_y ones(size(cankrig_dn_y))]);
    cankrig_trend_csr_mknz(ind_p) = trend(1);
    cankrig_trend_csr_mknz_pval(ind_p) = stats(3);

    [trend_2016,~,~,~,stats] = ...
        regress(cankrig_ann_5yrSigRem(cankrigtime_ealcoidx),...
        [cankrig_dn_y(cankrigtime_ealcoidx) ones(size(cankrig_dn_y(cankrigtime_ealcoidx)))]);
    cankrig_trend_csr_mknz_2016(ind_p) = trend_2016(1);
    cankrig_trend_csr_mknz_pval_2016(ind_p) = stats(3);

    clear tp_cs_all
end

%% CanKrig Trend Mackenzie Basin Map
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

%% Add stippling for non-significant trends
cankrig_trend_csr_mknz_pval(cankrig_trend_csr_mknz_pval>=0.05) = 1;
cankrig_trend_csr_mknz_pval(cankrig_trend_csr_mknz_pval<0.05) = 0;
cankrig_trend_csr_mknz_pval_2016(cankrig_trend_csr_mknz_pval_2016>=0.05) = 1;
cankrig_trend_csr_mknz_pval_2016(cankrig_trend_csr_mknz_pval_2016<0.05) = 0;

unsig_ind = find(cankrig_trend_csr_mknz_pval==1);

cankrig_trend_csr_mknz(unsig_ind) = nan;
figure
worldmap(plot_latlim,plot_lonlim)
setm(gca,'FLineWidth',1)
setm(gca,'FontSize',15)
pcolorm(lat_cankrig,lon_cankrig,cankrig_trend_csr_mknz*0.1) % mm to cm
clim([-5 5])
colormap(cmap)
geoshow(lakes,"FaceColor",'white','FaceAlpha',.1,'EdgeColor','black')
geoshow(subbasins,"FaceColor",'white','FaceAlpha',0.01,'EdgeColor','black')
geoshow(rivers,"Color","black")
bordersm('Canada')

unsig_ind = find(cankrig_trend_csr_mknz_pval_2016==1);

cankrig_trend_csr_mknz_2016(unsig_ind) = nan;
figure
worldmap(plot_latlim,plot_lonlim)
setm(gca,'FLineWidth',1)
setm(gca,'FontSize',15)
pcolorm(lat_cankrig,lon_cankrig,cankrig_trend_csr_mknz_2016*0.1) % mm to cm
clim([-5 5])
colormap(cmap)
geoshow(lakes,"FaceColor",'white','FaceAlpha',.1,'EdgeColor','black')
geoshow(subbasins,"FaceColor",'white','FaceAlpha',0.01,'EdgeColor','black')
geoshow(rivers,"Color","black")
bordersm('Canada')

