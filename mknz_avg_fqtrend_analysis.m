%% Mackenzie River Basin - Removal of Dominant Harmonic Signals
% The purpose of this script is to remove prominent harmonic trends from
% the respective precipitation, evaporation, runoff, and TWSA data within
% the Mackenzie River Basin in order to fit a long term linear trend
% without bias from seasonal or other non-linear trends. 
% Author: Stephanie Bringeland

close all
clear 
clc
 
%% Import 
basin_area = 1.80588e+12; %m^2
cankrig_p = importdata("cankrig_mknz_av.mat"); % mm
era5_p = importdata("era5_tp_mknz_av.mat"); % mm
era5ext_p = importdata("era5ext_tp_mknz.mat");
era5ext_et = importdata("era5ext_te_mknz.mat");
ealco_et = importdata("ealco_et_wf_mknz_av.mat"); % mm
ealco_e0 = importdata("ealco_e0_wf_mknz_av.mat"); % mm
era5_et = -1*importdata("era5_te_mknz_av.mat"); % mm
eccc_sw = importdata("eccc_flow_mknz_av.mat"); % m^3
csr = importdata("csr_mknz_av.mat"); % cm
jpl = importdata("jpl_mknz_av.mat"); % cm
jpl = jpl(1:end-1);
jpl_pred = importdata("automl_mknz_pred.mat"); % Gap-filled Data

% Define time for each parameter
cankrigtime = (datetime(2002,1,31):calmonths:datetime(2019,09,30))';
ealcotime = (datetime(2002,1,15):calmonths:datetime(2016,12,15))';
era5ptime = ncread("era5_e_tp_updated2024-04-05.nc",'time');
era5ptime = datetime(1900,1,1,era5ptime,0,0);era5ptime = era5ptime(1:end-2);
era5ettime = ncread("era5_e_tp_updated2024-04-05.nc",'time');
era5ettime = datetime(1900,1,1,era5ettime,0,0);era5ettime = era5ettime(1:end-2);
era5ext_p_time = importdata('era5ext_tp_time.mat');
era5ext_et_time = importdata('era5ext_te_time.mat');
eccc_swtime = importdata("eccc_mknz_dt.mat");
csrtime = ncread('CSR_GRACE_GRACE-FO_RL06_Mascons_all-corrections_v02.nc','time');
start_date = daysact('1-jan-2002');
csrtime = start_date + csrtime;
csrtime = datetime(csrtime,'ConvertFrom','datenum');
jpltime = ncread('GRCTellus.JPL.200204_202306.GLO.RL06.1M.MSCNv03CRI.nc','time');
start_date = daysact('1-jan-2002');
jpltime = start_date + jpltime;
jpltime = datetime(jpltime,'ConvertFrom','datenum');
jpltime = jpltime(1:end-1);
jplgaptime = importdata("JPL_gap_datetime.mat");

% Drop days from time (set to 15 for consistency)
[y,m,~] = ymd(cankrigtime);
cankrigtime = datetime(y,m,15);
[y,m,~] = ymd(ealcotime);
ealcotime = datetime(y,m,15);
[y,m,~] = ymd(era5ptime);
era5ptime = datetime(y,m,15);
[y,m,~] = ymd(era5ettime);
era5ettime = datetime(y,m,15);
[y,m,~] = ymd(era5ext_et_time);
era5ext_et_time = datetime(y,m,15);
[y,m,~] = ymd(era5ext_p_time);
era5ext_p_time = datetime(y,m,15);
[y,m,~] = ymd(eccc_swtime);
eccc_swtime = datetime(y,m,15);

%% Convert everything to cm EWH
cankrig_p = cankrig_p./10; % mm to cm
era5_p = era5_p./10; % mm to cm
ealco_et = ealco_et./10; % mm to cm
ealco_e0 = ealco_e0./10; % mm to cm
era5_e = era5_et./10; % mm to cm
eccc_sw = eccc_sw./basin_area; % m3 to m
eccc_sw = eccc_sw.*100; % m to cm
era5ext_et = era5ext_et./10; % mm to cm
era5ext_p = era5ext_p./10; % mm to cm
ealco_e = ealco_et + ealco_e0;

% Gapfilled JPL time series
jpl_gap_index = importdata("gap_index.txt");
jpl_counter = 1;
gap_counter = 1;
jpl_gapfilled = nan([length(jpl)+length(jpl_pred),1]);
for i = 1:length(jpl_gapfilled)
    if min(abs(jpl_gap_index-i)) == 0
        jpl_gapfilled(i,1) = jpl_pred(gap_counter);
        jpl_gapfilled_dt(i,1) = jplgaptime(gap_counter);
        gap_counter = gap_counter + 1;
    else 
        jpl_gapfilled(i,1) = jpl(jpl_counter);
        jpl_gapfilled_dt(i,1) = jpltime(jpl_counter);
        jpl_counter = jpl_counter + 1;
    end
end

%% Trend analysis
% CSR and JPL
csr_dn_y = datenum(csrtime)./365;
jpl_gf_dn_y = datenum(jpl_gapfilled_dt)./365;
jpl_dn_y = datenum(jpltime)./365;
csr_trend = polyfit(csr_dn_y,csr,1);
jpl_trend = polyfit(jpl_gf_dn_y,jpl_gapfilled,1);
[b,bint,r,rint,stats] = regress(jpl_gapfilled,[jpl_gf_dn_y ones(size(jpl_gf_dn_y,1))]);

% Least Squares - Annual
t = csr_dn_y;
Fc = 2*pi;
b1 = ones(size(t));
b2 = t;
b3 = sin(t*Fc);
b4 = cos(t*Fc);
B = [b1 b2 b3 b4];
x_soln = (B'*B)\(B'*csr);

sin_linear_csr = x_soln(1) + x_soln(2)*t +x_soln(3)*sin(t*Fc)+x_soln(4)*cos(t*Fc);
sin_csr = x_soln(3)*sin(t*Fc)+x_soln(4)*cos(t*Fc);
csr_annSigRem = csr - sin_csr;
csr_trend_annSigRem = polyfit(csr_dn_y,csr_annSigRem,1);

% Least Squares - 3.87 yr
t = csr_dn_y;
Fc = 2*pi/3.87;
b1 = ones(size(t));
b2 = t;
b3 = sin(t*Fc);
b4 = cos(t*Fc);
B = [b1 b2 b3 b4];
x_soln = (B'*B)\(B'*csr);

sin_linear_csr = x_soln(1) + x_soln(2)*t +x_soln(3)*sin(t*Fc)+x_soln(4)*cos(t*Fc);
sin_csr = x_soln(3)*sin(t*Fc)+x_soln(4)*cos(t*Fc);
csr_ann_4yrSigRem = csr_annSigRem - sin_csr;
csr_trend_ann_4yrSigRem = polyfit(csr_dn_y,csr_ann_4yrSigRem,1);

%% Least Squares - Annual
t = jpl_gf_dn_y;
Fc = 2*pi;
b1 = ones(size(t));
b2 = t;
b3 = sin(t*Fc);
b4 = cos(t*Fc);
B = [b1 b2 b3 b4];
x_soln = (B'*B)\(B'*jpl_gapfilled);

sin_linear_jpl = x_soln(1) + x_soln(2)*t +x_soln(3)*sin(t*Fc)+x_soln(4)*cos(t*Fc);
sin_jpl = x_soln(3)*sin(t*Fc)+x_soln(4)*cos(t*Fc);
jpl_annSigRem = jpl_gapfilled - sin_jpl;
jpl_trend_annSigRem = polyfit(jpl_gf_dn_y,jpl_annSigRem,1);

% Least Squares - 3.87 yr
t = jpl_gf_dn_y;
Fc = 2*pi/3.87;
b1 = ones(size(t));
b2 = t;
b3 = sin(t*Fc);
b4 = cos(t*Fc);
B = [b1 b2 b3 b4];
x_soln = (B'*B)\(B'*jpl_gapfilled);

sin_linear_jpl = x_soln(1) + x_soln(2)*t +x_soln(3)*sin(t*Fc)+x_soln(4)*cos(t*Fc);
sin_jpl = x_soln(3)*sin(t*Fc)+x_soln(4)*cos(t*Fc);
jpl_ann_4yrSigRem = jpl_annSigRem - sin_jpl;
jpl_trend_ann_4yrSigRem = polyfit(jpl_gf_dn_y,jpl_ann_4yrSigRem,1);

%% Evap
era5_e_mean_ind = isbetween(era5ettime,datetime(2004,1,1),datetime(2009,12,31));
era5_e_ealco_ind = isbetween(era5ettime,ealcotime(1),ealcotime(end));
era5e_mean = mean(era5_e(era5_e_ealco_ind));
era5_e_mr = era5_e-era5e_mean;
era5_e_cs = cumsum(era5_e_mr);
era5_e_mean_ind_mean = mean(era5_e_cs(era5_e_mean_ind));
era5_e_cs = era5_e_cs-era5_e_mean_ind_mean;

ealco_e_mean_ind = isbetween(ealcotime,datetime(2004,1,1),datetime(2009,12,31));
ealcoe_mean = mean(ealco_e);
ealco_e_mr = ealco_e-ealcoe_mean;
ealco_e_cs = cumsum(ealco_e_mr);
ealco_e_mean_ind_mean = mean(ealco_e_cs(ealco_e_mean_ind));
ealco_e_cs = ealco_e_cs-ealco_e_mean_ind_mean;

ealco_dn_y = datenum(ealcotime)./365;
era5e_dn_y = datenum(era5ettime)./365;
ealco_trend = polyfit(ealco_dn_y,ealco_e_cs,1);
era5e_trend = polyfit(era5e_dn_y(1:size(ealcotime,1)),era5_e_cs(1:size(ealcotime,1)),1);

% Least Squares - Annual
t = ealco_dn_y;
Fc = 2*pi;
b1 = ones(size(t));
b2 = t;
b3 = sin(t*Fc);
b4 = cos(t*Fc);
B = [b1 b2 b3 b4];
x_soln = (B'*B)\(B'*ealco_e_cs);

sin_linear_ealco = x_soln(1) + x_soln(2)*t +x_soln(3)*sin(t*Fc)+x_soln(4)*cos(t*Fc);
sin_ealco = x_soln(3)*sin(t*Fc)+x_soln(4)*cos(t*Fc);
ealco_annSigRem = ealco_e_cs - sin_ealco;
ealco_trend_annSigRem = polyfit(ealco_dn_y,ealco_annSigRem,1);

% Least squares - 6 month
t = ealco_dn_y;
Fc = 2*pi/0.5;
b1 = ones(size(t));
b2 = t;
b3 = sin(t*Fc);
b4 = cos(t*Fc);
B = [b1 b2 b3 b4];
x_soln = (B'*B)\(B'*ealco_e_cs);

sin_linear_ealco = x_soln(1) + x_soln(2)*t +x_soln(3)*sin(t*Fc)+x_soln(4)*cos(t*Fc);
sin_ealco = x_soln(3)*sin(t*Fc)+x_soln(4)*cos(t*Fc);
ealco_ann_6moSigRem = ealco_annSigRem - sin_ealco;
ealco_trend_ann_6moSigRem = polyfit(ealco_dn_y,ealco_ann_6moSigRem,1);

% Least Squares
t = era5e_dn_y(1:size(ealcotime,1));
Fc = 2*pi;
b1 = ones(size(t));
b2 = t;
b3 = sin(t*Fc);
b4 = cos(t*Fc);
B = [b1 b2 b3 b4];
x_soln = (B'*B)\(B'*era5_e_cs(1:size(ealcotime,1)));

sin_linear_era5e = x_soln(1) + x_soln(2)*t +x_soln(3)*sin(t*Fc)+x_soln(4)*cos(t*Fc);
sin_era5e = x_soln(3)*sin(era5e_dn_y*Fc)+x_soln(4)*cos(era5e_dn_y*Fc);
era5e_annSigRem = era5_e_cs - sin_era5e;
era5e_trend_annSigRem = polyfit(era5e_dn_y(1:size(ealcotime,1)),era5e_annSigRem(1:size(ealcotime,1)),1);

% Least Squares - 6 months
t = era5e_dn_y(1:size(ealcotime,1));
Fc = 2*pi/0.5;
b1 = ones(size(t));
b2 = t;
b3 = sin(t*Fc);
b4 = cos(t*Fc);
B = [b1 b2 b3 b4];
x_soln = (B'*B)\(B'*era5_e_cs(1:size(ealcotime,1)));

sin_linear_era5e = x_soln(1) + x_soln(2)*t +x_soln(3)*sin(t*Fc)+x_soln(4)*cos(t*Fc);
sin_era5e = x_soln(3)*sin(era5e_dn_y*Fc)+x_soln(4)*cos(era5e_dn_y*Fc);
era5e_ann_6moSigRem = era5e_annSigRem - sin_era5e;
era5e_trend_ann_6moSigRem = polyfit(era5e_dn_y(1:size(ealcotime,1)),era5e_ann_6moSigRem(1:size(ealcotime,1)),1);

%% Precip
era5_p_mean_ind = isbetween(era5ptime,datetime(2004,1,1),datetime(2009,12,31));
era5_p_ck_ind = isbetween(era5ptime,cankrigtime(1),cankrigtime(end));
era5p_mean = mean(era5_p(era5_p_ck_ind));
era5_p_mr = era5_p-era5p_mean;
era5_p_cs = cumsum(era5_p_mr);
era5_p_mean_ind_mean = mean(era5_p_cs(era5_p_mean_ind));
era5_p_cs = era5_p_cs-era5_p_mean_ind_mean;

cankrig_p_mean_ind = isbetween(cankrigtime,datetime(2004,1,1),datetime(2009,12,31));
cankrig_p_mean = mean(cankrig_p);
cankrig_p_mr = cankrig_p-cankrig_p_mean;
cankrig_p_cs = cumsum(cankrig_p_mr);
cankrig_p_mean_ind_mean = mean(cankrig_p_cs(cankrig_p_mean_ind));
cankrig_p_cs = cankrig_p_cs-cankrig_p_mean_ind_mean;

cankrig_dn_y = datenum(cankrigtime)./365;
era5p_dn_y = datenum(era5ptime)./365;
cankrig_trend = polyfit(cankrig_dn_y,cankrig_p_cs,1);
era5p_trend = polyfit(era5p_dn_y(1:size(cankrigtime,1)),era5_p_cs(1:size(cankrigtime,1)),1);

% Least Squares - Annual
t = cankrig_dn_y;
Fc = 2*pi;
b1 = ones(size(t));
b2 = t;
b3 = sin(t*Fc);
b4 = cos(t*Fc);
B = [b1 b2 b3 b4];
x_soln = (B'*B)\(B'*cankrig_p_cs);

sin_linear_cankrig = x_soln(1) + x_soln(2)*t +x_soln(3)*sin(t*Fc)+x_soln(4)*cos(t*Fc);
sin_cankrig = x_soln(3)*sin(t*Fc)+x_soln(4)*cos(t*Fc);
cankrig_annSigRem = cankrig_p_cs - sin_cankrig;
cankrig_trend_annSigRem = polyfit(cankrig_dn_y,cankrig_annSigRem,1);

% Least Squares - 5 years
t = cankrig_dn_y;
Fc = 2*pi/5.4;
b1 = ones(size(t));
b2 = t;
b3 = sin(t*Fc);
b4 = cos(t*Fc);
B = [b1 b2 b3 b4];
x_soln = (B'*B)\(B'*cankrig_p_cs);

sin_linear_cankrig = x_soln(1) + x_soln(2)*t +x_soln(3)*sin(t*Fc)+x_soln(4)*cos(t*Fc);
sin_cankrig = x_soln(3)*sin(t*Fc)+x_soln(4)*cos(t*Fc);
cankrig_ann_5yrSigRem = cankrig_annSigRem - sin_cankrig;
cankrig_trend_ann_5yrSigRem = polyfit(cankrig_dn_y,cankrig_ann_5yrSigRem,1);

% Least Squares - Annual
t = era5p_dn_y(1:size(cankrigtime,1));
Fc = 2*pi;
b1 = ones(size(t));
b2 = t;
b3 = sin(t*Fc);
b4 = cos(t*Fc);
B = [b1 b2 b3 b4];
x_soln = (B'*B)\(B'*era5_p_cs(1:size(cankrigtime,1)));

sin_linear_era5p = x_soln(1) + x_soln(2)*t +x_soln(3)*sin(t*Fc)+x_soln(4)*cos(t*Fc);
sin_era5p = x_soln(3)*sin(era5p_dn_y*Fc)+x_soln(4)*cos(era5p_dn_y*Fc);
era5p_annSigRem = era5_p_cs - sin_era5p;
era5p_trend_annSigRem = polyfit(era5p_dn_y(1:size(cankrigtime,1)),era5p_annSigRem(1:size(cankrigtime,1)),1);

% Least Squares - 5 year
t = era5p_dn_y(1:size(cankrigtime,1));
Fc = 2*pi/4.7;
b1 = ones(size(t));
b2 = t;
b3 = sin(t*Fc);
b4 = cos(t*Fc);
B = [b1 b2 b3 b4];
x_soln = (B'*B)\(B'*era5_p_cs(1:size(cankrigtime,1)));

sin_linear_era5p = x_soln(1) + x_soln(2)*t +x_soln(3)*sin(t*Fc)+x_soln(4)*cos(t*Fc);
sin_era5p = x_soln(3)*sin(era5p_dn_y*Fc)+x_soln(4)*cos(era5p_dn_y*Fc);
era5p_ann_5yrSigRem = era5p_annSigRem - sin_era5p;
era5p_trend_ann_5yrSigRem = polyfit(era5p_dn_y(1:size(cankrigtime,1)),era5p_ann_5yrSigRem(1:size(cankrigtime,1)),1);

%% Runoff
eccc_sw_mean_ind = isbetween(eccc_swtime,datetime(2004,1,1),datetime(2009,12,31));
eccc_sw_mean = mean(eccc_sw);
eccc_sw_mr = eccc_sw-eccc_sw_mean;
eccc_sw_cs = cumsum(eccc_sw_mr);
eccc_sw_mean_ind_mean = mean(eccc_sw_cs(eccc_sw_mean_ind));
eccc_sw_cs = eccc_sw_cs-eccc_sw_mean_ind_mean;

eccc_sw_dn_y = datenum(eccc_swtime)./365;
eccc_sw_trend = polyfit(eccc_sw_dn_y,eccc_sw_cs,1);

% Least Squares
t = eccc_sw_dn_y;
Fc = 2*pi;
b1 = ones(size(t));
b2 = t;
b3 = sin(t*Fc);
b4 = cos(t*Fc);
B = [b1 b2 b3 b4];
x_soln = (B'*B)\(B'*eccc_sw_cs);

sin_linear_eccc = x_soln(1) + x_soln(2)*t +x_soln(3)*sin(t*Fc)+x_soln(4)*cos(t*Fc);
sin_eccc = x_soln(3)*sin(t*Fc)+x_soln(4)*cos(t*Fc);
eccc_annSigRem = eccc_sw_cs - sin_eccc;
eccc_trend_annSigRem = polyfit(eccc_sw_dn_y,eccc_annSigRem,1);

%% Tiled Layout plot
figure
tiledlayout(2,2,"TileSpacing","tight");
nexttile
plot(csrtime(csrtime<datetime(2017,07,01)),csr(csrtime<datetime(2017,07,01)),'Color',"#d95f02");
hold on
plot(csrtime(csrtime>datetime(2018,05,01)),csr(csrtime>datetime(2018,05,01)),'Color',"#d95f02");
plot(csrtime,csr_trend(1)*csr_dn_y+csr_trend(2),'--','Color',"#d95f02");
plot(jpltime(jpltime<datetime(2017,07,01)),jpl(jpltime<datetime(2017,07,01)),'Color',"#1b9e77");
plot(jpltime(jpltime>datetime(2018,05,01)),jpl(jpltime>datetime(2018,05,01)),'Color',"#1b9e77");
plot(jpltime,jpl_trend(1)*jpl_dn_y+jpl_trend(2),'--','Color',"#1b9e77");
plot(jplgaptime,jpl_pred,'*','Color',"#1b9e77");
legend('CSR-M_{res}','','CSR-M_{res} Trend','JPL-M_{res}','','JPL-M_{res} Trend',...
    'JPL-M_{res} Gap Predictions','Location','southwest')
xlabel('Year')
ylabel('cm EWH')
xlim([datetime(2002,4,1) datetime(2023,3,31)])
ylim([-20 15])
grid on

nexttile
plot(ealcotime,ealco_e_cs,'Color',"#e41a1c");
hold on
plot(ealcotime,ealco_trend(1)*ealco_dn_y+ealco_trend(2),'--','Color',"#e41a1c");
hold on
plot(era5ettime,era5_e_cs,'Color',"#377eb8");
hold on
plot(era5ettime(1:size(ealcotime,1)),era5e_trend(1)*era5e_dn_y(1:size(ealcotime,1))+era5e_trend(2),'--','Color',"#377eb8");
legend('EALCO-E \DeltaE_{res}','EALCO-E \DeltaE_{res} Trend',...
    'ERA5-E \DeltaE_{res}','ERA5-E \DeltaE_{res} Trend','Location','southwest')
xlabel('Year')
ylabel('cm EWH')
xlim([datetime(2002,4,1) datetime(2023,3,31)])
ylim([-20 15])
grid on

nexttile
plot(cankrigtime,cankrig_p_cs,'Color',"#ff7f00");
hold on
plot(cankrigtime,cankrig_trend(1)*cankrig_dn_y+cankrig_trend(2),'--','Color',"#ff7f00");
hold on
plot(era5ptime,era5_p_cs,'Color',"#377eb8");
hold on
plot(era5ptime(1:size(cankrigtime,1)),era5p_trend(1)*era5p_dn_y(1:size(cankrigtime,1))+era5p_trend(2),'--','Color',"#377eb8");
legend('CanKrig-P \DeltaP_{res}','CanKrig-P \DeltaP_{res} Trend','ERA5-P \DeltaP_{res}',...
    'ERA5-P \DeltaP_{res} Trend','Location','southwest')
xlabel('Year')
ylabel('cm EWH')
xlim([datetime(2002,4,1) datetime(2023,3,31)])
ylim([-20 15])
grid on

nexttile
plot(eccc_swtime,eccc_sw_cs,'Color',"#756bb1");
hold on
plot(eccc_swtime,eccc_sw_trend(1)*eccc_sw_dn_y+eccc_sw_trend(2),'--','Color',"#756bb1");
legend('ECCC-R \DeltaR_{res}','ECCC-R \DeltaR_{res} Trend','Location','southwest')
xlabel('Year')
ylabel('cm EWH')
xlim([datetime(2002,4,1) datetime(2023,3,31)])
ylim([-20 15])
grid on

%% Tiled Layout plot
jpl_pred_trendrm = jpl_ann_4yrSigRem(jpl_gap_index);
figure
tiledlayout(2,2,"TileSpacing","tight");
nexttile
plot(csrtime(csrtime<datetime(2017,07,01)),csr_ann_4yrSigRem(csrtime<datetime(2017,07,01)),'Color',"#d95f02");
hold on
plot(csrtime(csrtime>datetime(2018,05,01)),csr_ann_4yrSigRem(csrtime>datetime(2018,05,01)),'Color',"#d95f02");
plot(csrtime,csr_trend_ann_4yrSigRem(1)*csr_dn_y+csr_trend_ann_4yrSigRem(2),'--','Color',"#d95f02");
plot(jpltime(jpltime<datetime(2017,07,01)),jpl_ann_4yrSigRem(jpltime<datetime(2017,07,01)),'Color',"#1b9e77");
plot(jpltime(jpltime>datetime(2018,05,01)),jpl_ann_4yrSigRem(jpltime>datetime(2018,05,01)),'Color',"#1b9e77");
plot(jpl_gapfilled_dt,jpl_trend_ann_4yrSigRem(1)*jpl_gf_dn_y+jpl_trend_ann_4yrSigRem(2),'--','Color',"#1b9e77");
plot(jplgaptime,jpl_pred_trendrm,'*','Color',"#1b9e77");
legend('CSR-M_{res}','','CSR-M_{res} Trend','JPL-M_{res}','','JPL-M_{res} Trend',...
    'JPL-M_{res} Gap Predictions','Location','southwest')
xlabel('Year')
ylabel('cm EWH')
xlim([datetime(2002,4,1) datetime(2023,3,31)])
ylim([-12 12])
grid on

nexttile
plot(ealcotime,ealco_ann_6moSigRem,'Color',"#e41a1c");
hold on
plot(ealcotime,ealco_trend_ann_6moSigRem(1)*ealco_dn_y+ealco_trend_ann_6moSigRem(2),'--','Color',"#e41a1c");
hold on
plot(era5ettime,era5e_ann_6moSigRem,'Color',"#377eb8");
hold on
plot(era5ettime(1:size(ealcotime,1)),era5e_trend_ann_6moSigRem(1)*era5e_dn_y(1:size(ealcotime,1))+era5e_trend_ann_6moSigRem(2),'--','Color',"#377eb8");
legend('EALCO-E \DeltaE_{res}','EALCO-E \DeltaE_{res} Trend',...
    'ERA5-E \DeltaE_{res}','ERA5-E \DeltaE_{res} Trend','Location','southwest')
xlabel('Year')
ylabel('cm EWH')
xlim([datetime(2002,4,1) datetime(2023,3,31)])
ylim([-12 12])
grid on

nexttile
plot(cankrigtime,cankrig_ann_5yrSigRem,'Color',"#ff7f00");
hold on
plot(cankrigtime,cankrig_trend_ann_5yrSigRem(1)*cankrig_dn_y+cankrig_trend_ann_5yrSigRem(2),'--','Color',"#ff7f00");
hold on
plot(era5ptime,era5p_ann_5yrSigRem,'Color',"#377eb8");
hold on
plot(era5ptime(1:size(cankrigtime,1)),era5p_trend_ann_5yrSigRem(1)*era5p_dn_y(1:size(cankrigtime,1))+era5p_trend_ann_5yrSigRem(2),'--','Color',"#377eb8");
legend('CanKrig-P \DeltaP_{res}','CanKrig-P \DeltaP_{res} Trend','ERA5-P \DeltaP_{res}',...
    'ERA5-P \DeltaP_{res} Trend','Location','southwest')
xlabel('Year')
ylabel('cm EWH')
xlim([datetime(2002,4,1) datetime(2023,3,31)])
ylim([-12 12])
grid on

nexttile
plot(eccc_swtime,eccc_annSigRem,'Color',"#756bb1");
hold on
plot(eccc_swtime,eccc_trend_annSigRem(1)*eccc_sw_dn_y+eccc_trend_annSigRem(2),'--','Color',"#756bb1");
legend('ECCC-R \DeltaR_{res}','ECCC-R \DeltaR_{res} Trend','Location','southwest')
xlabel('Year')
ylabel('cm EWH')
xlim([datetime(2002,4,1) datetime(2023,3,31)])
ylim([-12 12])
grid on

%% Tiled Layout - PSD 
twomonths = seconds(days(61));
figure
tiledlayout(4,2,'TileSpacing','tight');
nexttile
[pxx,fvec] = plomb(csr,csrtime);
semilogy(fvec*1e9,pxx,'-k')
ylim([100 11^10])
xlim([3 180])
% legend('CSR');
nexttile
[pxx,fvec] = plomb(jpl_gapfilled,jpl_gapfilled_dt); 
semilogy(fvec*1e9,pxx,'-k')
ylim([100 11^10])
xlim([3 180])
% legend('JPL');
nexttile
[pxx,fvec] = plomb(era5_e_cs,era5ettime); 
semilogy(fvec*1e9,pxx,'-k')
ylim([100 11^10])
xlim([3 180])
% legend('ERA5-E');
nexttile
[pxx,fvec] = plomb(ealco_e_cs,ealcotime); 
semilogy(fvec*1e9,pxx,'-k')
ylim([100 11^10])
xlim([3 180])
% legend('EALCO');
nexttile
[pxx,fvec] = plomb(era5_p_cs,era5ptime); 
semilogy(fvec*1e9,pxx,'-k')
ylim([100 11^10])
xlim([3 180])
% legend('ERA5-P');
nexttile
[pxx,fvec] = plomb(cankrig_p_cs,cankrigtime); 
semilogy(fvec*1e9,pxx,'-k')
ylim([100 11^10])
xlim([3 180])
% legend('CanKrig');
nexttile
[pxx,fvec] = plomb(eccc_sw_cs,eccc_swtime); 
semilogy(fvec*1e9,pxx,'-k')
ylim([100 11^10])
xlim([3 180])
% legend('ECCC');

