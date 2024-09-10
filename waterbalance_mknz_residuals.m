%% Water Balance for the Mackenzie River Basin
% The purpose of the script is to compute the water balance (WB = P - E -
% R) for the Mackenzie River Basin using different combinations of
% hydrological parameters.

close all
clear 
clc

%% Import 
basin_area = 1.80588e+12; %m^2
cankrig_tp = importdata("cankrig_mknz_av.mat"); % mm
era5_tp = importdata("era5_tp_mknz_av.mat"); % mm
era5ext_tp = importdata("era5ext_tp_mknz.mat");
era5ext_et = importdata("era5ext_te_mknz.mat");
ealco_et = importdata("ealco_et_wf_mknz_av.mat"); % mm
ealco_e0 = importdata("ealco_e0_wf_mknz_av.mat"); % mm
era5_e = -1*importdata("era5_te_mknz_av.mat"); % mm
eccc_sw = importdata("eccc_flow_mknz_av.mat"); % m^3
csr = importdata("csr_mknz_av.mat"); % cm
jpl = importdata("jpl_mknz_av.mat"); % cm
jpl_pred = importdata("automl_mknz_pred.mat"); % cm

% Define time for each parameter
cankrigtime = (datetime(2002,1,31):calmonths:datetime(2019,09,30))';
ealcotime = (datetime(2002,1,15):calmonths:datetime(2016,12,15))';
era5tptime = ncread("era5_e_tp_updated2024-04-05.nc",'time');
era5tptime = datetime(1900,1,1,era5tptime,0,0);era5tptime = era5tptime(1:end-2);
era5ettime = ncread("era5_e_tp_updated2024-04-05.nc",'time');
era5ettime = datetime(1900,1,1,era5ettime,0,0);era5ettime = era5ettime(1:end-2);
era5ext_tp_time = importdata('era5ext_tp_time.mat');
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
jpltime(112) = jpltime(112) - caldays(6);
jpltime(145) = jpltime(145) + caldays(6);
jplgaptime = importdata("JPL_gap_datetime.mat");
[y,m,~] = ymd(jplgaptime);
jplgaptime = datetime(y,m,15);

%% Gap-filled JPL time series
jpl_gap_index = importdata("gap_index.txt");
jpl_counter = 1;
gap_counter = 1;
jpl_gapfilled = nan([length(jpl)+length(jpl_pred),1]);
for i = 1:length(jpl_gapfilled)
    if min(abs(jpl_gap_index-i)) == 0
        jpl_gapfilled(i) = jpl_pred(gap_counter);
        jpl_gapfilled_dt(i) = jplgaptime(gap_counter);
        gap_counter = gap_counter + 1;
    else 
        jpl_gapfilled(i) = jpl(jpl_counter);
        jpl_gapfilled_dt(i) = jpltime(jpl_counter);
        jpl_counter = jpl_counter + 1;
    end
end

jpl_gapfilled_dt = jpl_gapfilled_dt';
jpl_gapfilled_tt = timetable(jpl_gapfilled_dt,jpl_gapfilled);

% Drop days from time (set to 15 for consistency)
[y,m,~] = ymd(cankrigtime);
cankrigtime = datetime(y,m,15);
[y,m,~] = ymd(ealcotime);
ealcotime = datetime(y,m,15);
[y,m,~] = ymd(era5tptime);
era5tptime = datetime(y,m,15);
[y,m,~] = ymd(era5ettime);
era5ettime = datetime(y,m,15);
[y,m,~] = ymd(era5ext_et_time);
era5ext_et_time = datetime(y,m,15);
[y,m,~] = ymd(era5ext_tp_time);
era5ext_tp_time = datetime(y,m,15);
[y,m,~] = ymd(eccc_swtime);
eccc_swtime = datetime(y,m,15);

%% Convert everything to cm EWH
cankrig_tp = cankrig_tp./10; % mm to cm
era5_tp = era5_tp./10; % mm to cm
ealco_et = ealco_et./10; % mm to cm
ealco_e0 = ealco_e0./10; % mm to cm
era5_e = era5_e./10; % mm to cm
eccc_sw = eccc_sw./basin_area; % m3 to m
eccc_sw = eccc_sw.*100; % m to cm
era5ext_et = era5ext_et./10; % mm to cm
era5ext_tp = era5ext_tp./10; % mm to cm
ealco_e = ealco_et + ealco_e0;

%% Water balance 1 - ERA5 TP, ERA5 ET
% Synchronize times
wb1a_tt = timetable(era5tptime,era5_tp,'VariableNames',{'dataA'});
wb1b_tt = timetable(era5ettime,era5_e,'VariableNames',{'dataB'});
wb1c_tt = timetable(eccc_swtime,eccc_sw,'VariableNames',{'dataC'});
wb1_tt = synchronize(wb1a_tt,wb1b_tt,wb1c_tt);
wb1_tt = rmmissing(wb1_tt);
wb1 = wb1_tt.dataA - wb1_tt.dataB - wb1_tt.dataC;
wb1_time = wb1_tt.Properties.RowTimes;
% Remove mean consistent with mean removed from GRACE mascon data
wb1_mean_ind = isbetween(wb1_time,datetime(2004,1,1),datetime(2009,12,31));
wb1_ealco_ind = isbetween(wb1_time,ealcotime(1),ealcotime(end));
wb1_mean = mean(wb1(wb1_ealco_ind));
wb1_mr = wb1-wb1_mean;
wb1_cs = cumsum(wb1_mr);
wb1_mean = mean(wb1_cs(wb1_mean_ind));
wb1_cs = wb1_cs-wb1_mean;
figure
plot(wb1_time,wb1_cs,'Color','blue');
hold on
plot(csrtime(csrtime<datetime(2017,07,01)),csr(csrtime<datetime(2017,07,01)),'Color',"#e41a1c");
hold on
plot(csrtime(csrtime>datetime(2018,05,01)),csr(csrtime>datetime(2018,05,01)),'Color',"#e41a1c");
hold on
plot(jpltime(jpltime<datetime(2017,07,01)),jpl(jpltime<datetime(2017,07,01)),'Color',"#4daf4a");
hold on
plot(jpltime(jpltime>datetime(2018,05,01)),jpl(jpltime>datetime(2018,05,01)),'Color',"#4daf4a");
hold on
plot(jplgaptime,jpl_pred,'*','Color',"#4daf4a");
legend('Water Budget TWSA','CSR TWSA','','JPL TWSA','','JPL TWSA Gap Predictions')
xlabel('Year')
ylabel('cm equivalent water height')
ylim([-25 15])

%% Water balance 2 - CanKrig + ERA5 TP, EALCO + ERA5 ET
% Synchronize times
tp = [cankrig_tp;era5ext_tp];
tp_time = [cankrigtime;era5ext_tp_time];
et = [ealco_e;era5ext_et];
et_time = [ealcotime;era5ext_et_time];
wb2a_tt = timetable(tp_time,tp,'VariableNames',{'dataA'});
wb2b_tt = timetable(et_time,et,'VariableNames',{'dataB'});
wb2c_tt = timetable(eccc_swtime,eccc_sw,'VariableNames',{'dataC'});
wb2_tt = synchronize(wb2a_tt,wb2b_tt,wb2c_tt);
wb2_tt = rmmissing(wb2_tt);
wb2 = wb2_tt.dataA - wb2_tt.dataB - wb2_tt.dataC;
wb2_time = wb2_tt.Properties.RowTimes;

% Remove mean consistent with mean removed from GRACE mascon data
wb2_mean_ind = isbetween(wb2_time,datetime(2004,1,1),datetime(2009,12,31));
wb2_ealco_ind = isbetween(wb2_time,ealcotime(1),ealcotime(end));
wb2_mean = mean(wb2(wb2_ealco_ind));
wb2_mr = wb2-wb2_mean;
wb2_cs = cumsum(wb2_mr);
wb2_mean = mean(wb2_cs(wb2_mean_ind));
wb2_cs = wb2_cs-wb2_mean;
figure
plot(wb2_time,wb2_cs,'Color',"blue");
hold on
plot(csrtime(csrtime<datetime(2017,07,01)),csr(csrtime<datetime(2017,07,01)),'Color',"#e41a1c");
hold on
plot(csrtime(csrtime>datetime(2018,05,01)),csr(csrtime>datetime(2018,05,01)),'Color',"#e41a1c");
hold on
plot(jpltime(jpltime<datetime(2017,07,01)),jpl(jpltime<datetime(2017,07,01)),'Color',"#4daf4a");
hold on
plot(jpltime(jpltime>datetime(2018,05,01)),jpl(jpltime>datetime(2018,05,01)),'Color',"#4daf4a");
hold on
plot(jplgaptime,jpl_pred,'*','Color',"#4daf4a")
legend('Water Budget TWSA','CSR TWSA','','JPL TWSA','','JPL TWSA Gap Predictions')
ylim([-25 15])

%% Water balance 3 - CanKrig TP, EALCO ET
% Synchronize times
wb3a_tt = timetable(cankrigtime,cankrig_tp,'VariableNames',{'dataA'});
wb3b_tt = timetable(ealcotime,ealco_e,'VariableNames',{'dataB'});
wb3c_tt = timetable(eccc_swtime,eccc_sw,'VariableNames',{'dataC'});
wb3_tt = synchronize(wb3a_tt,wb3b_tt,wb3c_tt);
wb3_tt = rmmissing(wb3_tt);
wb3 = wb3_tt.dataA - wb3_tt.dataB - wb3_tt.dataC;
wb3_time = wb3_tt.Properties.RowTimes;

% Remove mean consistent with mean removed from GRACE mascon data
wb3_mean_ind = isbetween(wb3_time,datetime(2004,1,1),datetime(2009,12,31));
wb3_ealco_ind = isbetween(wb3_time,ealcotime(1),ealcotime(end));
wb3_mean = mean(wb3(wb3_ealco_ind));
wb3_mr = wb3-wb3_mean;
wb3_cs = cumsum(wb3_mr);
wb3_mean = mean(wb3_cs(wb3_mean_ind));
wb3_cs = wb3_cs-wb3_mean;
figure
plot(wb3_time,wb3_cs,'Color',"blue");
hold on
plot(csrtime(csrtime<datetime(2017,07,01)),csr(csrtime<datetime(2017,07,01)),'Color',"#e41a1c");
hold on
plot(csrtime(csrtime>datetime(2018,05,01)),csr(csrtime>datetime(2018,05,01)),'Color',"#e41a1c");
hold on
plot(jpltime(jpltime<datetime(2017,07,01)),jpl(jpltime<datetime(2017,07,01)),'Color',"#4daf4a");
hold on
plot(jpltime(jpltime>datetime(2018,05,01)),jpl(jpltime>datetime(2018,05,01)),'Color',"#4daf4a");
hold on
plot(jplgaptime,jpl_pred,'*','Color',"#4daf4a")
legend('Water Budget TWSA','CSR TWSA','','JPL TWSA','','JPL TWSA Gap Predictions')
ylim([-25 15])

%% Water balance 4 - ERA5 TP, EALCO + ERA5 ET
% Synchronize times
wb4a_tt = timetable(era5tptime,era5_tp,'VariableNames',{'dataA'});
et = [ealco_e;era5ext_et];
et_time = [ealcotime;era5ext_et_time];
wb4b_tt = timetable(et_time,et,'VariableNames',{'dataB'});
wb4c_tt = timetable(eccc_swtime,eccc_sw,'VariableNames',{'dataC'});
wb4_tt = synchronize(wb4a_tt,wb4b_tt,wb4c_tt);
wb4_tt = rmmissing(wb4_tt);
wb4 = wb4_tt.dataA - wb4_tt.dataB - wb4_tt.dataC;
wb4_time = wb4_tt.Properties.RowTimes;

% Remove mean consistent with mean removed from GRACE mascon data
wb4_mean_ind = isbetween(wb4_time,datetime(2004,1,1),datetime(2009,12,31));
wb4_ealco_ind = isbetween(wb4_time,ealcotime(1),ealcotime(end));
wb4_mean = mean(wb4(wb4_ealco_ind));
wb4_mr = wb4-wb4_mean;
wb4_cs = cumsum(wb4_mr);
wb4_mean = mean(wb4_cs(wb4_mean_ind));
wb4_cs = wb4_cs-wb4_mean;
figure
plot(wb4_time,wb4_cs,'Color',"blue");
hold on
plot(csrtime(csrtime<datetime(2017,07,01)),csr(csrtime<datetime(2017,07,01)),'Color',"#e41a1c");
hold on
plot(csrtime(csrtime>datetime(2018,05,01)),csr(csrtime>datetime(2018,05,01)),'Color',"#e41a1c");
hold on
plot(jpltime(jpltime<datetime(2017,07,01)),jpl(jpltime<datetime(2017,07,01)),'Color',"#4daf4a");
hold on
plot(jpltime(jpltime>datetime(2018,05,01)),jpl(jpltime>datetime(2018,05,01)),'Color',"#4daf4a");
hold on
plot(jplgaptime,jpl_pred,'*','Color',"#4daf4a")
legend('Water Budget TWSA','CSR TWSA','','JPL TWSA','','JPL TWSA Gap Predictions')
ylim([-25 15])

%% Water balance 5 - ERA5 TP, EALCO ET
% Synchronize times
wb5a_tt = timetable(era5tptime,era5_tp,'VariableNames',{'dataA'});
wb5b_tt = timetable(ealcotime,ealco_e,'VariableNames',{'dataB'});
wb5c_tt = timetable(eccc_swtime,eccc_sw,'VariableNames',{'dataC'});
wb5_tt = synchronize(wb5a_tt,wb5b_tt,wb5c_tt);
wb5_tt = rmmissing(wb5_tt);
wb5 = wb5_tt.dataA - wb5_tt.dataB - wb5_tt.dataC;
wb5_time = wb5_tt.Properties.RowTimes;

% Remove mean consistent with mean removed from GRACE mascon data
wb5_mean_ind = isbetween(wb5_time,datetime(2004,1,1),datetime(2009,12,31));
wb5_ealco_ind = isbetween(wb5_time,ealcotime(1),ealcotime(end));
wb5_mean = mean(wb5(wb5_ealco_ind));
wb5_mr = wb5-wb5_mean;
wb5_cs = cumsum(wb5_mr);
wb5_mean = mean(wb5_cs(wb5_mean_ind));
wb5_cs = wb5_cs-wb5_mean;
figure
plot(wb5_time,wb5_cs,'Color',"blue");
hold on
plot(csrtime(csrtime<datetime(2017,07,01)),csr(csrtime<datetime(2017,07,01)),'Color',"#e41a1c");
hold on
plot(csrtime(csrtime>datetime(2018,05,01)),csr(csrtime>datetime(2018,05,01)),'Color',"#e41a1c");
hold on
plot(jpltime(jpltime<datetime(2017,07,01)),jpl(jpltime<datetime(2017,07,01)),'Color',"#4daf4a");
hold on
plot(jpltime(jpltime>datetime(2018,05,01)),jpl(jpltime>datetime(2018,05,01)),'Color',"#4daf4a");
hold on
plot(jplgaptime,jpl_pred,'*','Color',"#4daf4a")
legend('Water Budget TWSA','CSR TWSA','','JPL TWSA','','JPL TWSA Gap Predictions')
ylim([-25 15])

%% Water balance 6 - CanKrig TP, ERA5 ET
% Synchronize times
wb6a_tt = timetable(cankrigtime,cankrig_tp,'VariableNames',{'dataA'});
wb6b_tt = timetable(era5ettime,era5_e,'VariableNames',{'dataB'});
wb6c_tt = timetable(eccc_swtime,eccc_sw,'VariableNames',{'dataC'});
wb6_tt = synchronize(wb6a_tt,wb6b_tt,wb6c_tt);
wb6_tt = rmmissing(wb6_tt);
wb6 = wb6_tt.dataA - wb6_tt.dataB - wb6_tt.dataC;
wb6_time = wb6_tt.Properties.RowTimes;

% Remove mean consistent with mean removed from GRACE mascon data
wb6_mean_ind = isbetween(wb6_time,datetime(2004,1,1),datetime(2009,12,31));
wb6_ealco_ind = isbetween(wb6_time,ealcotime(1),ealcotime(end));
wb6_mean = mean(wb6(wb6_ealco_ind));
wb6_mr = wb6-wb6_mean;
wb6_cs = cumsum(wb6_mr);
wb6_mean = mean(wb6_cs(wb6_mean_ind));
wb6_cs = wb6_cs-wb6_mean;
figure
plot(wb6_time,wb6_cs,'Color',"blue");
hold on
plot(csrtime(csrtime<datetime(2017,07,01)),csr(csrtime<datetime(2017,07,01)),'Color',"#e41a1c");
hold on
plot(csrtime(csrtime>datetime(2018,05,01)),csr(csrtime>datetime(2018,05,01)),'Color',"#e41a1c");
hold on
plot(jpltime(jpltime<datetime(2017,07,01)),jpl(jpltime<datetime(2017,07,01)),'Color',"#4daf4a");
hold on
plot(jpltime(jpltime>datetime(2018,05,01)),jpl(jpltime>datetime(2018,05,01)),'Color',"#4daf4a");
hold on
plot(jplgaptime,jpl_pred,'*','Color',"#4daf4a")
legend('Water Budget TWSA','CSR TWSA','','JPL TWSA','','JPL TWSA Gap Predictions')
ylim([-25 15])

%% Water balance 7 - CanKrig + ERA5 TP, ERA5 ET
% Synchronize times
tp = [cankrig_tp;era5ext_tp];
tp_time = [cankrigtime;era5ext_tp_time];
wb7a_tt = timetable(tp_time,tp,'VariableNames',{'dataA'});
wb7b_tt = timetable(era5ettime,era5_e,'VariableNames',{'dataB'});
wb7c_tt = timetable(eccc_swtime,eccc_sw,'VariableNames',{'dataC'});
wb7_tt = synchronize(wb7a_tt,wb7b_tt,wb7c_tt);
wb7_tt = rmmissing(wb7_tt);
wb7 = wb7_tt.dataA - wb7_tt.dataB - wb7_tt.dataC;
wb7_time = wb7_tt.Properties.RowTimes;

% Remove mean consistent with mean removed from GRACE mascon data
wb7_mean_ind = isbetween(wb7_time,datetime(2004,1,1),datetime(2009,12,31));
wb7_ealco_ind = isbetween(wb7_time,ealcotime(1),ealcotime(end));
wb7_mean = mean(wb7(wb7_ealco_ind));
wb7_mr = wb7-wb7_mean;
wb7_cs = cumsum(wb7_mr);
wb7_mean = mean(wb7_cs(wb7_mean_ind));
wb7_cs = wb7_cs-wb7_mean;
figure
plot(wb7_time,wb7_cs,'Color',"blue");
hold on
plot(csrtime(csrtime<datetime(2017,07,01)),csr(csrtime<datetime(2017,07,01)),'Color',"#e41a1c");
hold on
plot(csrtime(csrtime>datetime(2018,05,01)),csr(csrtime>datetime(2018,05,01)),'Color',"#e41a1c");
hold on
plot(jpltime(jpltime<datetime(2017,07,01)),jpl(jpltime<datetime(2017,07,01)),'Color',"#4daf4a");
hold on
plot(jpltime(jpltime>datetime(2018,05,01)),jpl(jpltime>datetime(2018,05,01)),'Color',"#4daf4a");
hold on
plot(jplgaptime,jpl_pred,'*','Color',"#4daf4a")
legend('Water Budget TWSA','CSR TWSA','','JPL TWSA','','JPL TWSA Gap Predictions')
ylim([-25 15])

%% Water balance 8 - CanKrig TP, EALCO + ERA5 ET
% Synchronize times
wb8a_tt = timetable(cankrigtime,cankrig_tp,'VariableNames',{'dataA'});
et = [ealco_e;era5ext_et];
et_time = [ealcotime;era5ext_et_time];
wb8b_tt = timetable(et_time,et,'VariableNames',{'dataB'});
wb8c_tt = timetable(eccc_swtime,eccc_sw,'VariableNames',{'dataC'});
wb8_tt = synchronize(wb8a_tt,wb8b_tt,wb8c_tt);
wb8_tt = rmmissing(wb8_tt);
wb8 = wb8_tt.dataA - wb8_tt.dataB - wb8_tt.dataC;
wb8_time = wb8_tt.Properties.RowTimes;

% Remove mean consistent with mean removed from GRACE mascon data
wb8_mean_ind = isbetween(wb8_time,datetime(2004,1,1),datetime(2009,12,31));
wb8_ealco_ind = isbetween(wb8_time,ealcotime(1),ealcotime(end));
wb8_mean = mean(wb8(wb8_ealco_ind));
wb8_mr = wb8-wb8_mean;
wb8_cs = cumsum(wb8_mr);
wb8_mean = mean(wb8_cs(wb8_mean_ind));
wb8_cs = wb8_cs-wb8_mean;
figure
plot(wb8_time,wb8_cs,'Color',"blue");
hold on
plot(csrtime(csrtime<datetime(2017,07,01)),csr(csrtime<datetime(2017,07,01)),'Color',"#e41a1c");
hold on
plot(csrtime(csrtime>datetime(2018,05,01)),csr(csrtime>datetime(2018,05,01)),'Color',"#e41a1c");
hold on
plot(jpltime(jpltime<datetime(2017,07,01)),jpl(jpltime<datetime(2017,07,01)),'Color',"#4daf4a");
hold on
plot(jpltime(jpltime>datetime(2018,05,01)),jpl(jpltime>datetime(2018,05,01)),'Color',"#4daf4a");
hold on
plot(jplgaptime,jpl_pred,'*','Color',"#4daf4a")
legend('Water Budget TWSA','CSR TWSA','','JPL TWSA','','JPL TWSA Gap Predictions')
xlabel('Year')
ylabel('cm equivalent water height')
ylim([-25 15])

%% RMSE and Correlation Coefficient
wb1_cs_tt = timetable(wb1_time,wb1_cs);
wb1_cs_csr = retime(wb1_cs_tt,csrtime(csrtime<=wb1_time(end)),"linear");
wb1_cs_jpl = retime(wb1_cs_tt,jpl_gapfilled_dt(jpl_gapfilled_dt<=wb1_time(end)),"linear");
wb1_csr_rmse = rmse(wb1_cs_csr.wb1_cs,csr(csrtime<=wb1_time(end)));
wb1_jpl_rmse = rmse(wb1_cs_jpl.wb1_cs,jpl_gapfilled(jpl_gapfilled_dt<=wb1_time(end)));
rmse_wb1 = mean([wb1_jpl_rmse,wb1_csr_rmse]);
wb1_csr_cc = corrcoef(wb1_cs_csr.wb1_cs,csr(csrtime<=wb1_time(end)));
wb1_jpl_cc = corrcoef(wb1_cs_jpl.wb1_cs,jpl_gapfilled(jpl_gapfilled_dt<=wb1_time(end)));
cc_wb1 = mean([wb1_jpl_cc(1,2),wb1_csr_cc(1,2)]);

wb2_cs_tt = timetable(wb2_time,wb2_cs);
wb2_cs_csr = retime(wb2_cs_tt,csrtime(csrtime<=wb2_time(end)),"linear");
wb2_cs_jpl = retime(wb2_cs_tt,jpltime(jpltime<=wb2_time(end)),"linear");
wb2_csr_rmse = rmse(wb2_cs_csr.wb2_cs,csr(csrtime<=wb2_time(end)));
wb2_jpl_rmse = rmse(wb2_cs_jpl.wb2_cs,jpl(jpltime<=wb2_time(end)));
rmse_wb2 = mean([wb2_jpl_rmse,wb2_csr_rmse]);
wb2_csr_cc = corrcoef(wb2_cs_csr.wb2_cs,csr(csrtime<=wb2_time(end)));
wb2_jpl_cc = corrcoef(wb2_cs_jpl.wb2_cs,jpl(jpltime<=wb2_time(end)));
cc_wb2 = mean([wb2_jpl_cc(1,2),wb2_csr_cc(1,2)]);

wb3_cs_tt = timetable(wb3_time,wb3_cs);
wb3_cs_csr = retime(wb3_cs_tt,csrtime(csrtime<=wb3_time(end)),"linear");
wb3_cs_jpl = retime(wb3_cs_tt,jpltime(jpltime<=wb3_time(end)),"linear");
wb3_csr_rmse = rmse(wb3_cs_csr.wb3_cs,csr(csrtime<=wb3_time(end)));
wb3_jpl_rmse = rmse(wb3_cs_jpl.wb3_cs,jpl(jpltime<=wb3_time(end)));
rmse_wb3 = mean([wb3_jpl_rmse,wb3_csr_rmse]);
wb3_csr_cc = corrcoef(wb3_cs_csr.wb3_cs,csr(csrtime<=wb3_time(end)));
wb3_jpl_cc = corrcoef(wb3_cs_jpl.wb3_cs,jpl(jpltime<=wb3_time(end)));
cc_wb3 = mean([wb3_jpl_cc(1,2),wb3_csr_cc(1,2)]);

wb4_cs_tt = timetable(wb4_time,wb4_cs);
wb4_cs_csr = retime(wb4_cs_tt,csrtime(csrtime<=wb4_time(end)),"linear");
wb4_cs_jpl = retime(wb4_cs_tt,jpltime(jpltime<=wb4_time(end)),"linear");
wb4_csr_rmse = rmse(wb4_cs_csr.wb4_cs,csr(csrtime<=wb4_time(end)));
wb4_jpl_rmse = rmse(wb4_cs_jpl.wb4_cs,jpl(jpltime<=wb4_time(end)));
rmse_wb4 = mean([wb4_jpl_rmse,wb4_csr_rmse]);
wb4_csr_cc = corrcoef(wb4_cs_csr.wb4_cs,csr(csrtime<=wb4_time(end)));
wb4_jpl_cc = corrcoef(wb4_cs_jpl.wb4_cs,jpl(jpltime<=wb4_time(end)));
cc_wb4 = mean([wb4_jpl_cc(1,2),wb4_csr_cc(1,2)]);

wb5_cs_tt = timetable(wb5_time,wb5_cs);
wb5_cs_csr = retime(wb5_cs_tt,csrtime(csrtime<=wb5_time(end)),"linear");
wb5_cs_jpl = retime(wb5_cs_tt,jpltime(jpltime<=wb5_time(end)),"linear");
wb5_csr_rmse = rmse(wb5_cs_csr.wb5_cs,csr(csrtime<=wb5_time(end)));
wb5_jpl_rmse = rmse(wb5_cs_jpl.wb5_cs,jpl(jpltime<=wb5_time(end)));
rmse_wb5 = mean([wb5_jpl_rmse,wb5_csr_rmse]);
wb5_csr_cc = corrcoef(wb5_cs_csr.wb5_cs,csr(csrtime<=wb5_time(end)));
wb5_jpl_cc = corrcoef(wb5_cs_jpl.wb5_cs,jpl(jpltime<=wb5_time(end)));
cc_wb5 = mean([wb5_jpl_cc(1,2),wb5_csr_cc(1,2)]);

wb6_cs_tt = timetable(wb6_time,wb6_cs);
wb6_cs_csr = retime(wb6_cs_tt,csrtime(csrtime<=wb6_time(end)),"linear");
wb6_cs_jpl = retime(wb6_cs_tt,jpltime(jpltime<=wb6_time(end)),"linear");
wb6_csr_rmse = rmse(wb6_cs_csr.wb6_cs,csr(csrtime<=wb6_time(end)));
wb6_jpl_rmse = rmse(wb6_cs_jpl.wb6_cs,jpl(jpltime<=wb6_time(end)));
rmse_wb6 = mean([wb6_jpl_rmse,wb6_csr_rmse]);
wb6_csr_cc = corrcoef(wb6_cs_csr.wb6_cs,csr(csrtime<=wb6_time(end)));
wb6_jpl_cc = corrcoef(wb6_cs_jpl.wb6_cs,jpl(jpltime<=wb6_time(end)));
cc_wb6 = mean([wb6_jpl_cc(1,2),wb6_csr_cc(1,2)]);

wb7_cs_tt = timetable(wb7_time,wb7_cs);
wb7_cs_csr = retime(wb7_cs_tt,csrtime(csrtime<=wb7_time(end)),"linear");
wb7_cs_jpl = retime(wb7_cs_tt,jpltime(jpltime<=wb7_time(end)),"linear");
wb7_csr_rmse = rmse(wb7_cs_csr.wb7_cs,csr(csrtime<=wb7_time(end)));
wb7_jpl_rmse = rmse(wb7_cs_jpl.wb7_cs,jpl(jpltime<=wb7_time(end)));
rmse_wb7 = mean([wb7_jpl_rmse,wb7_csr_rmse]);
wb7_csr_cc = corrcoef(wb7_cs_csr.wb7_cs,csr(csrtime<=wb7_time(end)));
wb7_jpl_cc = corrcoef(wb7_cs_jpl.wb7_cs,jpl(jpltime<=wb7_time(end)));
cc_wb7 = mean([wb7_jpl_cc(1,2),wb7_csr_cc(1,2)]);

wb8_cs_tt = timetable(wb8_time,wb8_cs);
wb8_cs_csr = retime(wb8_cs_tt,csrtime(csrtime<=wb8_time(end)),"linear");
wb8_cs_jpl = retime(wb8_cs_tt,jpltime(jpltime<=wb8_time(end)),"linear");
wb8_csr_rmse = rmse(wb8_cs_csr.wb8_cs,csr(csrtime<=wb8_time(end)));
wb8_jpl_rmse = rmse(wb8_cs_jpl.wb8_cs,jpl(jpltime<=wb8_time(end)));
rmse_wb8 = mean([wb8_jpl_rmse,wb8_csr_rmse]);
wb8_csr_cc = corrcoef(wb8_cs_csr.wb8_cs,csr(csrtime<=wb8_time(end)));
wb8_jpl_cc = corrcoef(wb8_cs_jpl.wb8_cs,jpl(jpltime<=wb8_time(end)));
cc_wb8 = mean([wb8_jpl_cc(1,2),wb8_csr_cc(1,2)]);

%% Compute cross-correlation
% wb1
wbdaysdt = (wb1_cs_tt.wb1_time(1):caldays(1):wb1_cs_tt.wb1_time(end))';
wb1_days = interp1(wb1_cs_tt.wb1_time,wb1_cs_tt.wb1_cs,wbdaysdt);
jpldaysdt = (jpl_gapfilled_dt(1):caldays(1):jpl_gapfilled_dt(end))';
jpl_days = interp1(jpl_gapfilled_dt,jpl_gapfilled,jpldaysdt);

ind_wb1 = isbetween(wbdaysdt,jpldaysdt(1),wbdaysdt(end));
ind_jpl = isbetween(jpldaysdt,jpldaysdt(1),wbdaysdt(end));
wb1_days = wb1_days(ind_wb1);
jpl1_days = jpl_days(ind_jpl);
[xcc,lags] = xcorr(wb1_days,jpl1_days,182,'normalized');
[xcc_wb1,ind] = max(xcc);
lag_wb1 = lags(ind);

% wb2
wbdaysdt = (wb2_cs_tt.wb2_time(1):caldays(1):wb2_cs_tt.wb2_time(end))';
wb2_days = interp1(wb2_cs_tt.wb2_time,wb2_cs_tt.wb2_cs,wbdaysdt);

ind_wb2 = isbetween(wbdaysdt,jpldaysdt(1),wbdaysdt(end));
ind_jpl = isbetween(jpldaysdt,jpldaysdt(1),wbdaysdt(end));
wb2_days = wb2_days(ind_wb2);
jpl2_days = jpl_days(ind_jpl);
[xcc,lags] = xcorr(wb2_days,jpl2_days,182,'normalized');
[xcc_wb2,ind] = max(xcc);
lag_wb2 = lags(ind);

% wb3
wbdaysdt = (wb3_cs_tt.wb3_time(1):caldays(1):wb3_cs_tt.wb3_time(end))';
wb3_days = interp1(wb3_cs_tt.wb3_time,wb3_cs_tt.wb3_cs,wbdaysdt);

ind_wb3 = isbetween(wbdaysdt,jpldaysdt(1),wbdaysdt(end));
ind_jpl = isbetween(jpldaysdt,jpldaysdt(1),wbdaysdt(end));
wb3_days = wb3_days(ind_wb3);
jpl3_days = jpl_days(ind_jpl);
[xcc,lags] = xcorr(wb3_days,jpl3_days,182,'normalized');
[xcc_wb3,ind] = max(xcc);
lag_wb3 = lags(ind);

% wb4
wbdaysdt = (wb4_cs_tt.wb4_time(1):caldays(1):wb4_cs_tt.wb4_time(end))';
wb4_days = interp1(wb4_cs_tt.wb4_time,wb4_cs_tt.wb4_cs,wbdaysdt);

ind_wb4 = isbetween(wbdaysdt,jpldaysdt(1),wbdaysdt(end));
ind_jpl = isbetween(jpldaysdt,jpldaysdt(1),wbdaysdt(end));
wb4_days = wb4_days(ind_wb4);
jpl4_days = jpl_days(ind_jpl);
[xcc,lags] = xcorr(wb4_days,jpl4_days,182,'normalized');
[xcc_wb4,ind] = max(xcc);
lag_wb4 = lags(ind);

% wb5
wbdaysdt = (wb5_cs_tt.wb5_time(1):caldays(1):wb5_cs_tt.wb5_time(end))';
wb5_days = interp1(wb5_cs_tt.wb5_time,wb5_cs_tt.wb5_cs,wbdaysdt);

ind_wb5 = isbetween(wbdaysdt,jpldaysdt(1),wbdaysdt(end));
ind_jpl = isbetween(jpldaysdt,jpldaysdt(1),wbdaysdt(end));
wb5_days = wb5_days(ind_wb5);
jpl5_days = jpl_days(ind_jpl);
[xcc,lags] = xcorr(wb5_days,jpl5_days,182,'normalized');
[xcc_wb5,ind] = max(xcc);
lag_wb5 = lags(ind);

% wb6
wbdaysdt = (wb6_cs_tt.wb6_time(1):caldays(1):wb6_cs_tt.wb6_time(end))';
wb6_days = interp1(wb6_cs_tt.wb6_time,wb6_cs_tt.wb6_cs,wbdaysdt);

ind_wb6 = isbetween(wbdaysdt,jpldaysdt(1),wbdaysdt(end));
ind_jpl = isbetween(jpldaysdt,jpldaysdt(1),wbdaysdt(end));
wb6_days = wb6_days(ind_wb6);
jpl6_days = jpl_days(ind_jpl);
[xcc,lags] = xcorr(wb6_days,jpl6_days,182,'normalized');
[xcc_wb6,ind] = max(xcc);
lag_wb6 = lags(ind);

% wb7
wbdaysdt = (wb7_cs_tt.wb7_time(1):caldays(1):wb7_cs_tt.wb7_time(end))';
wb7_days = interp1(wb7_cs_tt.wb7_time,wb7_cs_tt.wb7_cs,wbdaysdt);

ind_wb7 = isbetween(wbdaysdt,jpldaysdt(1),wbdaysdt(end));
ind_jpl = isbetween(jpldaysdt,jpldaysdt(1),wbdaysdt(end));
wb7_days = wb7_days(ind_wb7);
jpl7_days = jpl_days(ind_jpl);
[xcc,lags] = xcorr(wb7_days,jpl7_days,182,'normalized');
[xcc_wb7,ind] = max(xcc);
lag_wb7 = lags(ind);

% wb8
wbdaysdt = (wb8_cs_tt.wb8_time(1):caldays(1):wb8_cs_tt.wb8_time(end))';
wb8_days = interp1(wb8_cs_tt.wb8_time,wb8_cs_tt.wb8_cs,wbdaysdt);

ind_wb8 = isbetween(wbdaysdt,jpldaysdt(1),wbdaysdt(end));
ind_jpl = isbetween(jpldaysdt,jpldaysdt(1),wbdaysdt(end));
wb8_days = wb8_days(ind_wb8);
jpl8_days = jpl_days(ind_jpl);
[xcc,lags] = xcorr(wb8_days,jpl8_days,182,'normalized');
[xcc_wb8,ind] = max(xcc);
lag_wb8 = lags(ind);

%% Compute and plot residuals
% Residuals
% WB1
wbdaysdt = (wb1_cs_tt.wb1_time(1):caldays(1):wb1_cs_tt.wb1_time(end))';
wb1_days = interp1(wb1_cs_tt.wb1_time,wb1_cs_tt.wb1_cs,wbdaysdt);
jpldaysdt = (jpl_gapfilled_dt(1):caldays(1):jpl_gapfilled_dt(end))';
jpl_days = interp1(jpl_gapfilled_dt,jpl_gapfilled,jpldaysdt);
jpldaysdt_lag1 = jpldaysdt+caldays(lag_wb1);
ind_wb1 = isbetween(wbdaysdt,jpldaysdt_lag1(1),wbdaysdt(end));
ind_jpl = isbetween(jpldaysdt_lag1,jpldaysdt_lag1(1),wbdaysdt(end));
wb1_days = wb1_days(ind_wb1);
wb1daysdt = wbdaysdt(ind_wb1);
jpl1_days = jpl_days(ind_jpl);
jpldaysdt_lag1 = jpldaysdt_lag1(ind_jpl);

wb1_cs_tt = timetable(wb1daysdt,wb1_days);
wb1_cs_jpl_tt = retime(wb1_cs_tt,jpldaysdt_lag1,"linear");

res_wb1 = jpl1_days - wb1_cs_jpl_tt.wb1_days;
[~,locs,~,p] = findpeaks(res_wb1);
wb1_peaks = wb1_cs_jpl_tt.wb1daysdt(locs);
wb1_peaks = wb1_peaks(p>2);
t = tiledlayout(3,3,"TileSpacing","tight");
nexttile
plot(wb1_cs_jpl_tt.wb1daysdt,res_wb1)
ylabel('cm EWH')
xlim([datetime(2002,04,01) datetime(2020,12,31)])
ylim([-20 10])

% WB5
wbdaysdt = (wb5_cs_tt.wb5_time(1):caldays(1):wb5_cs_tt.wb5_time(end))';
wb5_days = interp1(wb5_cs_tt.wb5_time,wb5_cs_tt.wb5_cs,wbdaysdt);
jpldaysdt = (jpl_gapfilled_dt(1):caldays(1):jpl_gapfilled_dt(end))';
jpl_days = interp1(jpl_gapfilled_dt,jpl_gapfilled,jpldaysdt);
jpldaysdt_lag5 = jpldaysdt+caldays(lag_wb5);
ind_wb5 = isbetween(wbdaysdt,jpldaysdt_lag5(1),wbdaysdt(end));
ind_jpl = isbetween(jpldaysdt_lag5,jpldaysdt_lag5(1),wbdaysdt(end));
wb5_days = wb5_days(ind_wb5);
wb5daysdt = wbdaysdt(ind_wb5);
jpl5_days = jpl_days(ind_jpl);
jpldaysdt_lag5 = jpldaysdt_lag5(ind_jpl);

wb5_cs_tt = timetable(wb5daysdt,wb5_days);
wb5_cs_jpl_tt = retime(wb5_cs_tt,jpldaysdt_lag5,"linear");

res_wb5 = jpl5_days - wb5_cs_jpl_tt.wb5_days;
[~,locs,~,p] = findpeaks(res_wb5);
wb5_peaks = wb5_cs_jpl_tt.wb5daysdt(locs);
wb5_peaks = wb5_peaks(p>2);
nexttile
plot(wb5_cs_jpl_tt.wb5daysdt,res_wb5)
ylabel('cm EWH')
xlim([datetime(2002,04,01) datetime(2020,12,31)])
ylim([-20 10])

% WB4
wbdaysdt = (wb4_cs_tt.wb4_time(1):caldays(1):wb4_cs_tt.wb4_time(end))';
wb4_days = interp1(wb4_cs_tt.wb4_time,wb4_cs_tt.wb4_cs,wbdaysdt);
jpldaysdt = (jpl_gapfilled_dt(1):caldays(1):jpl_gapfilled_dt(end))';
jpl_days = interp1(jpl_gapfilled_dt,jpl_gapfilled,jpldaysdt);
jpldaysdt_lag4 = jpldaysdt+caldays(lag_wb4);
ind_wb4 = isbetween(wbdaysdt,jpldaysdt_lag4(1),wbdaysdt(end));
ind_jpl = isbetween(jpldaysdt_lag4,jpldaysdt_lag4(1),wbdaysdt(end));
wb4_days = wb4_days(ind_wb4);
wb4daysdt = wbdaysdt(ind_wb4);
jpl4_days = jpl_days(ind_jpl);
jpldaysdt_lag4 = jpldaysdt_lag4(ind_jpl);

wb4_cs_tt = timetable(wb4daysdt,wb4_days);
wb4_cs_jpl_tt = retime(wb4_cs_tt,jpldaysdt_lag4,"linear");

res_wb4 = jpl4_days - wb4_cs_jpl_tt.wb4_days;
[~,locs,~,p] = findpeaks(res_wb4);
wb4_peaks = wb4_cs_jpl_tt.wb4daysdt(locs);
wb4_peaks = wb4_peaks(p>2);
nexttile
plot(wb4_cs_jpl_tt.wb4daysdt,res_wb4)
ylabel('cm EWH')
xlim([datetime(2002,04,01) datetime(2020,12,31)])
ylim([-20 10])

% WB6
wbdaysdt = (wb6_cs_tt.wb6_time(1):caldays(1):wb6_cs_tt.wb6_time(end))';
wb6_days = interp1(wb6_cs_tt.wb6_time,wb6_cs_tt.wb6_cs,wbdaysdt);
jpldaysdt = (jpl_gapfilled_dt(1):caldays(1):jpl_gapfilled_dt(end))';
jpl_days = interp1(jpl_gapfilled_dt,jpl_gapfilled,jpldaysdt);
jpldaysdt_lag6 = jpldaysdt+caldays(lag_wb6);
ind_wb6 = isbetween(wbdaysdt,jpldaysdt_lag6(1),wbdaysdt(end));
ind_jpl = isbetween(jpldaysdt_lag6,jpldaysdt_lag6(1),wbdaysdt(end));
wb6_days = wb6_days(ind_wb6);
wb6daysdt = wbdaysdt(ind_wb6);
jpl6_days = jpl_days(ind_jpl);
jpldaysdt_lag6 = jpldaysdt_lag6(ind_jpl);

wb6_cs_tt = timetable(wb6daysdt,wb6_days);
wb6_cs_jpl_tt = retime(wb6_cs_tt,jpldaysdt_lag6,"linear");

res_wb6 = jpl6_days - wb6_cs_jpl_tt.wb6_days;
[~,locs,~,p] = findpeaks(-res_wb6);
wb6_peaks = wb6_cs_jpl_tt.wb6daysdt(locs);
wb6_peaks = wb6_peaks(p>2);
nexttile
plot(wb6_cs_jpl_tt.wb6daysdt,res_wb6)
ylabel('cm EWH')
xlim([datetime(2002,04,01) datetime(2020,12,31)])
ylim([-20 10])

% WB3
wbdaysdt = (wb3_cs_tt.wb3_time(1):caldays(1):wb3_cs_tt.wb3_time(end))';
wb3_days = interp1(wb3_cs_tt.wb3_time,wb3_cs_tt.wb3_cs,wbdaysdt);
jpldaysdt = (jpl_gapfilled_dt(1):caldays(1):jpl_gapfilled_dt(end))';
jpl_days = interp1(jpl_gapfilled_dt,jpl_gapfilled,jpldaysdt);
jpldaysdt_lag3 = jpldaysdt+caldays(lag_wb3);
ind_wb3 = isbetween(wbdaysdt,jpldaysdt_lag3(1),wbdaysdt(end));
ind_jpl = isbetween(jpldaysdt_lag3,jpldaysdt_lag3(1),wbdaysdt(end));
wb3_days = wb3_days(ind_wb3);
wb3daysdt = wbdaysdt(ind_wb3);
jpl3_days = jpl_days(ind_jpl);
jpldaysdt_lag3 = jpldaysdt_lag3(ind_jpl);

wb3_cs_tt = timetable(wb3daysdt,wb3_days);
wb3_cs_jpl_tt = retime(wb3_cs_tt,jpldaysdt_lag3,"linear");

res_wb3 = jpl3_days - wb3_cs_jpl_tt.wb3_days;
[~,locs,~,p] = findpeaks(-res_wb3);
wb3_peaks = wb3_cs_jpl_tt.wb3daysdt(locs);
wb3_peaks = wb3_peaks(p>2);
nexttile
plot(wb3_cs_jpl_tt.wb3daysdt,res_wb3)
ylabel('cm EWH')
xlim([datetime(2002,04,01) datetime(2020,12,31)])
ylim([-20 10])

% WB8
wbdaysdt = (wb8_cs_tt.wb8_time(1):caldays(1):wb8_cs_tt.wb8_time(end))';
wb8_days = interp1(wb8_cs_tt.wb8_time,wb8_cs_tt.wb8_cs,wbdaysdt);
jpldaysdt = (jpl_gapfilled_dt(1):caldays(1):jpl_gapfilled_dt(end))';
jpl_days = interp1(jpl_gapfilled_dt,jpl_gapfilled,jpldaysdt);
jpldaysdt_lag8 = jpldaysdt+caldays(lag_wb8);
ind_wb8 = isbetween(wbdaysdt,jpldaysdt_lag8(1),wbdaysdt(end));
ind_jpl = isbetween(jpldaysdt_lag8,jpldaysdt_lag8(1),wbdaysdt(end));
wb8_days = wb8_days(ind_wb8);
wb8daysdt = wbdaysdt(ind_wb8);
jpl8_days = jpl_days(ind_jpl);
jpldaysdt_lag8 = jpldaysdt_lag8(ind_jpl);

wb8_cs_tt = timetable(wb8daysdt,wb8_days);
wb8_cs_jpl_tt = retime(wb8_cs_tt,jpldaysdt_lag8,"linear");

res_wb8 = jpl8_days - wb8_cs_jpl_tt.wb8_days;
[~,locs,~,p] = findpeaks(-res_wb8);
wb8_peaks = wb8_cs_jpl_tt.wb8daysdt(locs);
wb8_peaks = wb8_peaks(p>2);
nexttile
plot(wb8_cs_jpl_tt.wb8daysdt,res_wb8)
ylabel('cm EWH')
xlim([datetime(2002,04,01) datetime(2020,12,31)])
ylim([-20 10])


% WB7
wbdaysdt = (wb7_cs_tt.wb7_time(1):caldays(1):wb7_cs_tt.wb7_time(end))';
wb7_days = interp1(wb7_cs_tt.wb7_time,wb7_cs_tt.wb7_cs,wbdaysdt);
jpldaysdt = (jpl_gapfilled_dt(1):caldays(1):jpl_gapfilled_dt(end))';
jpl_days = interp1(jpl_gapfilled_dt,jpl_gapfilled,jpldaysdt);
jpldaysdt_lag7 = jpldaysdt+caldays(lag_wb7);
ind_wb7 = isbetween(wbdaysdt,jpldaysdt_lag7(1),wbdaysdt(end));
ind_jpl = isbetween(jpldaysdt_lag7,jpldaysdt_lag7(1),wbdaysdt(end));
wb7_days = wb7_days(ind_wb7);
wb7daysdt = wbdaysdt(ind_wb7);
jpl7_days = jpl_days(ind_jpl);
jpldaysdt_lag7 = jpldaysdt_lag7(ind_jpl);

wb7_cs_tt = timetable(wb7daysdt,wb7_days);
wb7_cs_jpl_tt = retime(wb7_cs_tt,jpldaysdt_lag7,"linear");

res_wb7 = jpl7_days - wb7_cs_jpl_tt.wb7_days;
[~,locs,~,p] = findpeaks(-res_wb7);
wb7_peaks = wb7_cs_jpl_tt.wb7daysdt(locs);
wb7_peaks = wb7_peaks(p>2);
nexttile
plot(wb7_cs_jpl_tt.wb7daysdt,res_wb7)
ylabel('cm EWH')
xlim([datetime(2002,04,01) datetime(2020,12,31)])
ylim([-20 10])

nexttile

% WB2
wbdaysdt = (wb2_cs_tt.wb2_time(1):caldays(1):wb2_cs_tt.wb2_time(end))';
wb2_days = interp1(wb2_cs_tt.wb2_time,wb2_cs_tt.wb2_cs,wbdaysdt);
jpldaysdt = (jpl_gapfilled_dt(1):caldays(1):jpl_gapfilled_dt(end))';
jpl_days = interp1(jpl_gapfilled_dt,jpl_gapfilled,jpldaysdt);
jpldaysdt_lag2 = jpldaysdt+caldays(lag_wb2);
ind_wb2 = isbetween(wbdaysdt,jpldaysdt_lag2(1),wbdaysdt(end));
ind_jpl = isbetween(jpldaysdt_lag2,jpldaysdt_lag2(1),wbdaysdt(end));
wb2_days = wb2_days(ind_wb2);
wb2daysdt = wbdaysdt(ind_wb2);
jpl2_days = jpl_days(ind_jpl);
jpldaysdt_lag2 = jpldaysdt_lag2(ind_jpl);

wb2_cs_tt = timetable(wb2daysdt,wb2_days);
wb2_cs_jpl_tt = retime(wb2_cs_tt,jpldaysdt_lag2,"linear");

res_wb2 = jpl2_days - wb2_cs_jpl_tt.wb2_days;
[~,locs,~,p] = findpeaks(-res_wb2);
wb2_peaks = wb2_cs_jpl_tt.wb2daysdt(locs);
wb2_peaks = wb2_peaks(p>2);
nexttile
plot(wb2_cs_jpl_tt.wb2daysdt,res_wb2)
ylabel('cm EWH')
xlim([datetime(2002,04,01) datetime(2020,12,31)])
ylim([-20 10])
