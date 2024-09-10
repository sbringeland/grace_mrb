function [spat_av_bsn] = spatavbasin_gen(bsn,grid,n_timesteps)
%% Basin Mask
canada = importdata('canada_extent_full_1.txt');
canada(canada==-9999) = nan;
canada(canada>=1) = 1;
%% Basin
default_nv = -9999; % Value assigned by ArcGIS to NV elements
u_bsn = unique(bsn);
bsn_basin_val = u_bsn(2);
bsn(bsn==default_nv)=0;
bsn(bsn==bsn_basin_val)=1;
bsn(bsn>1) = 1;

if size(grid,1) == 40
    bsn_can = bsn.*canada;
    bsn_can = bsn_can(:,~all(isnan(bsn_can)));
    bsn_can = bsn_can(~all(isnan(bsn_can),2),:);
    bsn_can(bsn_can==0)=nan;
else 
    bsn_can = bsn;
    bsn_can(bsn_can==0)=nan;
end
    
spat_av_bsn = ones(n_timesteps,1);
for i = 1:size(grid,3)
    bsn_clipped = bsn_can.*grid(:,:,i);
    bsn_clipped = bsn_clipped(:,~all(isnan(bsn_clipped)));
    bsn_clipped = bsn_clipped(~all(isnan(bsn_clipped),2),:);
    spat_av_bsn(i) = mean(bsn_clipped,[1 2],'omitnan');
end

