clear 
clc
close all
%% Import JPL L3 data
%Input NetCDF file and extract time and lwe thickness
jplraw = pagetranspose(ncread('GRCTellus.JPL.200204_202306.GLO.RL06.1M.MSCNv03CRI.nc','lwe_thickness'));
jpltime = ncread('GRCTellus.JPL.200204_202306.GLO.RL06.1M.MSCNv03CRI.nc','time');
start_date = daysact('1-jan-2002');
jpltime = start_date + jpltime;
jpltime = datetime(jpltime,'ConvertFrom','datenum');
%% Import CSR L3 data
%Input NetCDF file and extract time and lwe thickness
csrraw = pagetranspose(ncread('CSR_GRACE_GRACE-FO_RL06_Mascons_all-corrections_v02.nc','lwe_thickness'));
csrtime = ncread('CSR_GRACE_GRACE-FO_RL06_Mascons_all-corrections_v02.nc','time');
start_date = daysact('1-jan-2002');
csrtime = start_date + csrtime;
csrtime = datetime(csrtime,'ConvertFrom','datenum');

% Filter GRACE data using land/ocean mask - CSR
% Import land/ocean mask
landocean_msk = ncread("CSR_GRACE_GRACE-FO_RL06_Mascons_v02_LandMask.nc",'LO_val');
landocean_msk = flip(landocean_msk');
landocean_msk = [landocean_msk(:,720:1440) landocean_msk(:,1:719)];
csr = ones(size(csrraw));
for i = 1:size(csr,3)
    grace_page = csrraw(:,:,i);
    grace_page = [grace_page(:,720:1440) grace_page(:,1:719)];
    grace_page = flip(grace_page,1);
    csr(:,:,i) = grace_page.*landocean_msk;
end

%% Filter GRACE data using land/ocean mask - JPL
% Import land/ocean mask
landocean_msk = ncread("LAND_MASK.CRI.nc",'land_mask');
landocean_msk = flip(landocean_msk');
landocean_msk = [landocean_msk(:,360:720) landocean_msk(:,1:359)];
jpl = ones(size(jplraw));
for i = 1:size(jpl,3)
    grace_page = jplraw(:,:,i);
    grace_page = [grace_page(:,360:720) grace_page(:,1:359)];
    grace_page = flip(grace_page,1);
    jpl(:,:,i) = grace_page.*landocean_msk;
end

%% Load and apply basin masks
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

%% Canada mask
canada025(canada025<0)=nan;
canada025(canada025>=0)=1;
canada050(canada050<0)=nan;
canada050(canada050>=0)=1;
for  i = 1:size(csrtime,1)
    csr_page = csr(:,:,i);
    csr_page = csr_page.*canada025;
    csr_page = csr_page(:,~all(isnan(csr_page)));
    csr_page = csr_page(~all(isnan(csr_page),2),:);
    csr_can(:,:,i) = csr_page;
end
for  i = 1:size(jpltime,1)
    jpl_page = jpl(:,:,i);
    jpl_page = jpl_page.*canada050;
    jpl_page = jpl_page(:,~all(isnan(jpl_page)));
    jpl_page = jpl_page(~all(isnan(jpl_page),2),:);
    jpl_can(:,:,i) = jpl_page;
end
%% SubMackenzie masks
% Mackenzie basin
csr_mknz_av = spatavbasin_gen(mknzcan_qd,csr_can,size(csrtime,1));
jpl_mknz_av = spatavbasin_gen(mknzcan_hd,jpl_can,size(jpltime,1));

figure
plot(jpltime,jpl_mknz_av);
title('Mackenzie');
hold on
plot(csrtime,csr_mknz_av);
ylim([-30 30]);
legend('JPL','CSR')
