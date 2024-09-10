%% CSR GIA Correction Grids
close all
clear
clc

%% Import
info = ncinfo("CSR_GRACE_GRACE-FO_RL0602_Mascons_GIA-component.nc");
csr_gia_corr = pagetranspose(ncread('CSR_GRACE_GRACE-FO_RL0602_Mascons_GIA-component.nc','lwe_thickness'));
time = ncread('CSR_GRACE_GRACE-FO_RL0602_Mascons_GIA-component.nc','time');
start_date = daysact('1-jan-2002');
time = start_date + time;
csrtime = datetime(time,'ConvertFrom','datenum');
csr_dn_y = datenum(csrtime)./365;
ealcotime = (datetime(2002,1,15):calmonths:datetime(2016,12,15))';
csrtime_ealcoidx = isbetween(csrtime,ealcotime(1),ealcotime(end));
csrtime_era5idx = isbetween(csrtime,ealcotime(1),datetime(2022,10,31));
subbasins = readgeotable("mknz_outline_new_08_2023.shp");
lakes = readgeotable('worldlakes.shp');
rivers = readgeotable('worldrivers.shp');

%% Import land/ocean mask
landocean_msk = flip(transpose(ncread('CSR_GRACE_GRACE-FO_RL06_Mascons_v02_LandMask.nc','LO_val')));
landocean_msk = [landocean_msk(:,720:1440) landocean_msk(:,1:719)];
landocean_msk(landocean_msk==0) = nan;
gia_data_csr = ones(size(csr_gia_corr));
for i = 1:size(gia_data_csr,3)
    gia_page = csr_gia_corr(:,:,i);
    gia_page = [gia_page(:,720:1440) gia_page(:,1:719)];
    gia_page = flip(gia_page,1);
    gia_data_csr(:,:,i) = gia_page.*landocean_msk;
end
%% Canada mask
can_mask_qd = importdata("canada_extent_full_lakesin_025.txt");
can_mask_qd(can_mask_qd>=0) = 1;
can_mask_qd(can_mask_qd<0) = nan;

for k = 1:size(gia_data_csr,3)
    csr_can_page = gia_data_csr(:,:,k).*can_mask_qd;
    csr_can_page = csr_can_page(:,~all(isnan(csr_can_page)));
    csr_can_page = csr_can_page(~all(isnan(csr_can_page),2),:);
    gia_can(:,:,k) = csr_can_page;
end
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

for i = 1:size(gia_can,3)
    csr_can_page = gia_can(:,:,i);
    page = csr_can_page.*mknzcan_qd;
    page = page(:,~all(isnan(page)));
    page = page(~all(isnan(page),2),:);
    gia_mknz(:,:,i) = page;
end

%% Rate of change
gia_mknz_trend = nan(size(gia_mknz,[1 2]));
gia_mknz_trend_pval = nan(size(gia_mknz,[1 2]));

for i = 1:size(gia_mknz,1)
    for j = 1:size(gia_mknz,2)
        data = squeeze(gia_mknz(i,j,:));
        if isnan(data(1)) == true
            gia_mknz_trend(i,j) = nan;
        else
            [p,~,~,~,stats] = regress(data,[csr_dn_y ones(size(csr_dn_y))]);
            gia_mknz_trend(i,j) = p(1);
            gia_mknz_trend_pval(i,j) = stats(3);            
        end
    end
end

%% Map
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

gia_mknz_trend_pval(gia_mknz_trend_pval>=0.05) = 1;
gia_mknz_trend_pval(gia_mknz_trend_pval<0.05) = 0;

unsig_ind = find(gia_mknz_trend_pval==1);
lat_unsig = YQ_025(unsig_ind)+0.0125;
lon_unsig = XQ_025(unsig_ind);
latlon_unsig = [lat_unsig lon_unsig];
latlon_unsig = sortrows(latlon_unsig,1,"ascend");
latlon_unsig = latlon_unsig(1:2:end,:);
gia_mknz_trend(unsig_ind) = nan;

figure
worldmap([mknz_lat_range(1)-0.25,mknz_lat_range(2)+0.75],mknz_lon_range)
setm(gca,'FLineWidth',1)
setm(gca,'FontSize',12)
pcolorm(YQ_025+0.125,XQ_025,gia_mknz_trend)
clim([-5 5])
colormap(cmap)
geoshow(lakes,"FaceColor",'white','FaceAlpha',.1,'EdgeColor','black')
geoshow(subbasins,"FaceColor",'white','FaceAlpha',0.01,'EdgeColor','black')
geoshow(rivers,"Color","black")
bordersm('Canada')

