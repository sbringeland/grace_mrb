%% EALCO Evapotranspiration and Evaporation in Canada
close all
clear
clc
%% Import water fraction values
fileID = fopen('WF_255_Background_5KM');
wf = transpose(fread(fileID,[1140,960],"int8"));
wf(wf==-1) = nan;
wf = wf/100;
%% Import EALCO evaporation
ealco_time = datetime(2002,1,15):calmonths:datetime(2016,12,15);
[y,~,~] = ymd(ealco_time');
ealco_dn_y = datenum(ealco_time)/365;
% Import Coordinates
coords = importdata("Canada_ET/coordinates_ealco.txt");
lon = coords.data(:,4);
lat = coords.data(:,5);
lon = reshape(lon,[1140,960])';
lat = reshape(lat,[1140,960])';
xind = 1:2:size(lon,1)-1;
yind = 1:2:size(lon,2)-1;
lon = lon(xind,yind);
lat = lat(xind,yind);
oceanflag = -32760;
nanflag = -32750;

% Import TIFF Files
dir_name = 'Canada_E0';names = dir(dir_name);
foldernamescell = struct2cell(names);
foldernamescell = foldernamescell(1,:)';
foldernamescell = foldernamescell(3:end);
e0 = zeros([size(lat,1) size(lat,2) size(ealco_time,1)]);

for i = 1:length(foldernamescell)
    filename = char(foldernamescell(i));
    fileID = fopen(sprintf('Canada_E0/%s',filename));
    t = fread(fileID,[570,480],"int16");
    t(t==oceanflag)=nan;
    t(t==nanflag)=nan;
    t = 0.1.*t;
    e0(:,:,i) = t';
end

%% Import EALCO Evapotranspiration
% Import Coordinates
coords = importdata("Canada_ET/coordinates_ealco.txt");
lon = coords.data(:,4);
lat = coords.data(:,5);
lon = reshape(lon,[1140,960])';
lat = reshape(lat,[1140,960])';

% Import TIFF Files
dir_name = 'Canada_ET/*.tif';names = dir(dir_name);
foldernamescell = struct2cell(names);
foldernamescell = foldernamescell(1,:)';
et = zeros([size(lat,1) size(lat,2) size(ealco_time,1)]);
nanflag = -32760;

for i = 1:length(foldernamescell)
    filename = char(foldernamescell(i));
    t = Tiff(sprintf('Canada_ET/%s',filename),"r");
    et(:,:,i) = read(t);
end

et(et==nanflag) = nan;
et = et*0.1; % mm

%% Apply water fraction to evapotranspiration
et_wf = zeros(size(et));
for i = 1:size(ealco_time,2)
    et_wf(:,:,i) = et(:,:,i).*(1-wf);
end

%% Apply water fraction to evapotranspiration
wf_10=blockproc(wf, [2 2], @(block_struct) mean(block_struct.data(:),'omitnan'));
e0_wf = zeros(size(e0));
for i = 1:size(ealco_time,2)
    e0_wf(:,:,i) = e0(:,:,i).*(wf_10);
end


