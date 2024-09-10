%% Plotting trends in the MRB
close all
clear
clc

%% Import data
csr_can = importdata("csr_can_025.mat");
jpl_can = importdata("jpl_can_050.mat");
csrtime = ncread('CSR_GRACE_GRACE-FO_RL06_Mascons_all-corrections_v02.nc','time');
start_date = daysact('1-jan-2002');
csrtime = start_date + csrtime;
csrtime = datetime(csrtime,'ConvertFrom','datenum');
jpltime = ncread('GRCTellus.JPL.200204_202306.GLO.RL06.1M.MSCNv03CRI.nc','time');
start_date = daysact('1-jan-2002');
jpltime = start_date + jpltime;
jpltime = datetime(jpltime,'ConvertFrom','datenum');
ealcotime = (datetime(2002,1,15):calmonths:datetime(2016,12,15))';

csr_dn_y = datenum(csrtime)./365;
jpl_dn_y = datenum(jpltime)./365;

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

mknz050 = importdata("mknz_extent_050_lakesin.txt");
canada050 = importdata("canada_extent_full_lakesin_050.txt");
mknzcan_hd = mknz050.*canada050;
nanval = max(mknzcan_hd,[],'all');
mknzcan_hd(mknzcan_hd==nanval)=nan;
mknzcan_hd = mknzcan_hd(:,~all(isnan(mknzcan_hd)));
mknzcan_hd = mknzcan_hd(~all(isnan(mknzcan_hd),2),:);
mknzcan_hd(mknzcan_hd<=0) = nan;
mknzcan_hd(mknzcan_hd>0) = 1;

%% Get CSR linear trend for each pixel in the MRB
csrtime_ealcoidx = isbetween(csrtime,ealcotime(1),ealcotime(end));
csrtime_era5idx = isbetween(csrtime,ealcotime(1),datetime(2022,10,31));

for i = 1:size(csr_can,3)
    csr_can_page = csr_can(:,:,i);
    page = csr_can_page.*mknzcan_qd;
    page = page(:,~all(isnan(page)));
    page = page(~all(isnan(page),2),:);
    csr_mknz(:,:,i) = page;
end

csr_mknz_trend = nan(size(csr_mknz,[1 2]));
csr_mknz_trend_2016 = nan(size(csr_mknz,[1 2]));
csr_mknz_trend_2022 = nan(size(csr_mknz,[1 2]));
csr_mknz_trend_pval = nan(size(csr_mknz,[1 2]));
csr_mknz_trend_pval_2016 = nan(size(csr_mknz,[1 2]));
csr_mknz_trend_pval_2022 = nan(size(csr_mknz,[1 2]));

Fc_1yr = 2*pi;
a1 = ones(size(csr_dn_y));
a2 = csr_dn_y;
a3 = sin(csr_dn_y*Fc_1yr);
a4 = cos(csr_dn_y*Fc_1yr);
A = [a1 a2 a3 a4];

Fc_4yr = 2*pi/3.87;
b1 = ones(size(csr_dn_y));
b2 = csr_dn_y;
b3 = sin(csr_dn_y*Fc_4yr);
b4 = cos(csr_dn_y*Fc_4yr);
B = [b1 b2 b3 b4];
for i = 1:size(csr_mknz,1)
    for j = 1:size(csr_mknz,2)
        data = squeeze(csr_mknz(i,j,:));
        if isnan(data(1)) == true
            csr_mknz_trend(i,j) = nan;
            csr_mknz_trend_2016(i,j) = nan;
            csr_mknz_trend_2022(i,j) = nan;
        else
            % Least Squares - Annual            
            A_soln = (A'*A)\(A'*data);            
            sinA = A_soln(3)*sin(csr_dn_y*Fc_1yr)+A_soln(4)*cos(csr_dn_y*Fc_1yr);           
            % Least Squares - 3.87 yr
            B_soln = (B'*B)\(B'*data);            
            sinB = B_soln(3)*sin(csr_dn_y*Fc_4yr)+B_soln(4)*cos(csr_dn_y*Fc_4yr);
            csr_ann_4yrSigRem = data - sinA - sinB;

            [p,~,~,~,stats] = regress(csr_ann_4yrSigRem,[csr_dn_y ones(size(csr_dn_y))]);
            csr_mknz_trend(i,j) = p(1);
            csr_mknz_trend_pval(i,j) = stats(3);
            [p,~,~,~,stats] = regress(csr_ann_4yrSigRem(csrtime_ealcoidx),...
                [csr_dn_y(csrtime_ealcoidx) ones(size(csr_dn_y(csrtime_ealcoidx)))]);
            csr_mknz_trend_2016(i,j) = p(1);
            csr_mknz_trend_pval_2016(i,j) = stats(3);
            [p,~,~,~,stats] = regress(csr_ann_4yrSigRem(csrtime_era5idx),...
                [csr_dn_y(csrtime_era5idx) ones(size(csr_dn_y(csrtime_era5idx)))]);
            csr_mknz_trend_2022(i,j) = p(1);
            csr_mknz_trend_pval_2022(i,j) = stats(3);
        end
    end
end

%% Get JPL linear trend for each pixel in the MRB
jpltime_ealcoidx = isbetween(jpltime,ealcotime(1),ealcotime(end));
jpltime_era5idx = isbetween(jpltime,ealcotime(1),datetime(2022,10,31));

for i = 1:size(jpl_can,3)
    jpl_can_page = jpl_can(:,:,i);
    page = jpl_can_page.*mknzcan_hd;
    page = page(:,~all(isnan(page)));
    page = page(~all(isnan(page),2),:);
    jpl_mknz(:,:,i) = page;
end

jpl_mknz_trend = nan(size(jpl_mknz,[1 2]));
jpl_mknz_trend_2016 = nan(size(jpl_mknz,[1 2]));
jpl_mknz_trend_2022 = nan(size(jpl_mknz,[1 2]));
jpl_mknz_trend_pval = nan(size(jpl_mknz,[1 2]));
jpl_mknz_trend_pval_2016 = nan(size(jpl_mknz,[1 2]));
jpl_mknz_trend_pval_2022 = nan(size(jpl_mknz,[1 2]));

Fc_1yr = 2*pi;
a1 = ones(size(jpl_dn_y));
a2 = jpl_dn_y;
a3 = sin(jpl_dn_y*Fc_1yr);
a4 = cos(jpl_dn_y*Fc_1yr);
A = [a1 a2 a3 a4];

Fc_4yr = 2*pi/3.87;
b1 = ones(size(jpl_dn_y));
b2 = jpl_dn_y;
b3 = sin(jpl_dn_y*Fc_4yr);
b4 = cos(jpl_dn_y*Fc_4yr);
B = [b1 b2 b3 b4];

for i = 1:size(jpl_mknz,1)
    for j = 1:size(jpl_mknz,2)
        data = squeeze(jpl_mknz(i,j,:));
        if isnan(data(1)) == true
            jpl_mknz_trend(i,j) = nan;
            jpl_mknz_trend_2016(i,j) = nan;
            jpl_mknz_trend_2022(i,j) = nan;
        else
            % Least Squares - Annual            
            A_soln = (A'*A)\(A'*data);            
            sinA = A_soln(3)*sin(jpl_dn_y*Fc_1yr)+A_soln(4)*cos(jpl_dn_y*Fc_1yr);           
            % Least Squares - 3.87 yr
            B_soln = (B'*B)\(B'*data);            
            sinB = B_soln(3)*sin(jpl_dn_y*Fc_4yr)+B_soln(4)*cos(jpl_dn_y*Fc_4yr);
            jpl_ann_4yrSigRem = data - sinA - sinB;

            [p,~,~,~,stats] = regress(jpl_ann_4yrSigRem,[jpl_dn_y ones(size(jpl_dn_y))]);
            jpl_mknz_trend(i,j) = p(1);
            jpl_mknz_trend_pval(i,j) = stats(3);
            [p,~,~,~,stats] = regress(jpl_ann_4yrSigRem(jpltime_ealcoidx),...
                [jpl_dn_y(jpltime_ealcoidx) ones(size(jpl_dn_y(jpltime_ealcoidx)))]);
            jpl_mknz_trend_2016(i,j) = p(1);
            jpl_mknz_trend_pval_2016(i,j) = stats(3);
            [p,~,~,~,stats] = regress(jpl_ann_4yrSigRem(jpltime_era5idx),...
                [jpl_dn_y(jpltime_era5idx) ones(size(jpl_dn_y(jpltime_era5idx)))]);
            jpl_mknz_trend_2022(i,j) = p(1);
            jpl_mknz_trend_pval_2022(i,j) = stats(3);
        end
    end
end

% Plot within the MRB
% CSR
cmap = importdata("cmap_blueyellowred_symmetrical.mat");
[row,col] = find(mknzcan_qd==1);
lat_025 = flip(42.125:0.25:83.125);
lon_025 = -141:0.25:-53;
mknz_lon_range = sort([lon_025(min(col)) lon_025(max(col))]);
mknz_lat_range = sort([lat_025(min(row)) lat_025(max(row))]);
mknz_lat_range = [mknz_lat_range(1),mknz_lat_range(2)];
xq = mknz_lon_range(1):0.25:mknz_lon_range(2);
yq = mknz_lat_range(1):0.25:mknz_lat_range(2);
[XQ_025,YQ_025] = meshgrid(xq,yq);
YQ_025 = flip(YQ_025);

% JPL
[row,col] = find(mknzcan_hd==1);
lat = flip(41.5:0.5:82.5);
latlim = [min(lat) max(lat)];
lon = -140.75:0.5:-52.75;
lonlim = [min(lon) max(lon)];
mknz_lon_range = sort([lon(min(col)) lon(max(col))]);
mknz_lat_range = sort([lat(min(row)) lat(max(row))]);
mknz_lat_range = [mknz_lat_range(1),mknz_lat_range(2)];
xq = mknz_lon_range(1):0.50:mknz_lon_range(2);
yq = mknz_lat_range(1):0.50:mknz_lat_range(2);
[XQ_050,YQ_050] = meshgrid(xq,yq);
plot_latlim = [mknz_lat_range(1)-0.25 mknz_lat_range(2)+0.75];
plot_lonlim = [mknz_lon_range(1)-0.25 mknz_lon_range(2)+1.25];
YQ_050 = flip(YQ_050);

%% Remove insignificant trends - CSR
csr_mknz_trend_pval_2016(csr_mknz_trend_pval_2016>=0.05) = 1;
csr_mknz_trend_pval_2016(csr_mknz_trend_pval_2016<0.05) = 0;
csr_mknz_trend_pval_2022(csr_mknz_trend_pval_2022>=0.05) = 1;
csr_mknz_trend_pval_2022(csr_mknz_trend_pval_2022<0.05) = 0;

unsig_ind = find(csr_mknz_trend_pval_2016==1);

csr_mknz_trend_2016(unsig_ind) = nan;
figure
worldmap([mknz_lat_range(1)-0.25,mknz_lat_range(2)+0.75],mknz_lon_range)
setm(gca,'FLineWidth',1)
setm(gca,'FontSize',15)
pcolorm(YQ_025+0.125,XQ_025,csr_mknz_trend_2016)
clim([-5 5])
colormap(cmap)
geoshow(lakes,"FaceColor",'white','FaceAlpha',.1,'EdgeColor','black')
geoshow(subbasins,"FaceColor",'white','FaceAlpha',0.01,'EdgeColor','black')
geoshow(rivers,"Color","black")
bordersm('Canada')
bordersm('Canada')

unsig_ind = find(csr_mknz_trend_pval_2022==1);

csr_mknz_trend_2022(unsig_ind) = nan;
figure
worldmap([mknz_lat_range(1)-0.25,mknz_lat_range(2)+0.75],mknz_lon_range)
setm(gca,'FLineWidth',1)
setm(gca,'FontSize',15)
pcolorm(YQ_025+0.125,XQ_025,csr_mknz_trend_2022)
clim([-5 5])
colormap(cmap)
geoshow(lakes,"FaceColor",'white','FaceAlpha',.1,'EdgeColor','black')
geoshow(subbasins,"FaceColor",'white','FaceAlpha',0.01,'EdgeColor','black')
geoshow(rivers,"Color","black")
bordersm('Canada')

%% Remove insignificant trends - JPL
jpl_mknz_trend_pval_2016(jpl_mknz_trend_pval_2016>=0.05) = 1;
jpl_mknz_trend_pval_2016(jpl_mknz_trend_pval_2016<0.05) = 0;
jpl_mknz_trend_pval_2022(jpl_mknz_trend_pval_2022>=0.05) = 1;
jpl_mknz_trend_pval_2022(jpl_mknz_trend_pval_2022<0.05) = 0;

unsig_ind = find(jpl_mknz_trend_pval_2016==1);

jpl_mknz_trend_2016(unsig_ind) = nan;
figure
worldmap([plot_latlim(1)+0.25,plot_latlim(2)+0.625],plot_lonlim)
setm(gca,'FLineWidth',1)
setm(gca,'FontSize',15)
pcolorm(YQ_050+0.50,XQ_050-0.25,jpl_mknz_trend_2016)
clim([-5 5])
colormap(cmap)
geoshow(lakes,"FaceColor",'white','FaceAlpha',.1,'EdgeColor','black')
geoshow(subbasins,"FaceColor",'white','FaceAlpha',0.01,'EdgeColor','black')
geoshow(rivers,"Color","black")
bordersm('Canada')

unsig_ind = find(jpl_mknz_trend_pval_2022==1);

jpl_mknz_trend_2022(unsig_ind) = nan;
figure
worldmap([plot_latlim(1)+0.25,plot_latlim(2)+0.625],plot_lonlim)
setm(gca,'FLineWidth',1)
setm(gca,'FontSize',15)
pcolorm(YQ_050+0.50,XQ_050-0.25,jpl_mknz_trend_2022)
clim([-5 5])
colormap(cmap)
geoshow(lakes,"FaceColor",'white','FaceAlpha',.1,'EdgeColor','black')
geoshow(subbasins,"FaceColor",'white','FaceAlpha',0.01,'EdgeColor','black')
geoshow(rivers,"Color","black")
bordersm('Canada')
