%% ERA 5 Evapotranspiration
close all
clear
clc 
%% Import
et = ncread("era5_e_tp_updated2024-04-05.nc",'e');

for j = 1:size(et,3)
    era5_et_25deg(:,:,j) = imresize(et(:,:,j),0.40,"box");
end
