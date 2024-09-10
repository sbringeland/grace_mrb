%% ERA 5 Precipitation
close all
clear
clc 
%% Import
tp = ncread("era5_e_tp_updated2024-04-05.nc",'tp');

for j = 1:size(tp,3)
    era5_tp_25deg(:,:,j) = imresize(tp(:,:,j),0.40,"box");
end
