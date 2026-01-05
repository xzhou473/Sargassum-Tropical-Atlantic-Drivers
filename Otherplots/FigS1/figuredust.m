clear all
input_file = '/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/utils/newdata/duaod550_f64_res075_2003-2022_GASB_unmasked.nc';
%yeasm_before=load('/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/utils/fldmean_data/yseasm_1999-2011_mld_newdata.mat');
%yeasm_after=load('/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/utils/fldmean_data/yseasm_2011-2022_mld_newdata.mat');
yeasm_before=ncread('/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/utils/newdata/yseasm_2003-2011_dust.nc','duaod550');
yeasm_after=ncread('/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/utils/newdata/yseasm_2011-2022_dust.nc','duaod550');
change=yeasm_after-yeasm_before;
yseasm=yeasm_before;
mycm = table2array(readtable('/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/scripts/test_mycm.csv'));
lon = ncread(input_file, 'longitude');
lat = ncread(input_file, 'latitude');

data_var = 'duaod550'; % Change to the actual variable name in NetCDF
data = ncread(input_file, data_var);
aa=find(data<-30000);
data(aa)=NaN;

%%

TAB = xlsread('/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/utils/digitized_contour/Winter_Map_Jouanno_et_al_2023_DigitizedManual.xlsx');

lon_contour = TAB(:,1); 
lat_contour = TAB(:,2); 

%resultspath3 = '../outputs/simplified_timeser_changes/yseasm_changes_2011-2019_vs_2003-2011_dust.nc'
%ncid = netcdf.open(resultspath3)

%  time = double(netcdf.getVar(ncid,0));  
%  lon = double(netcdf.getVar(ncid,1));   
%  lat = double(netcdf.getVar(ncid,2));   
%netcdf.close(ncid);

[Xd Yd] = meshgrid(lon,lat);
Xd = Xd';
Yd = Yd';

% percentual changes of yseasm computed as 100*(yseasm_changes/yseam_before)
win_perc = 100*(change(:,:,1)./yseasm(:,:,1));
spr_perc = 100*(change(:,:,2)./yseasm(:,:,2));
summ_perc = 100*(change(:,:,3)./yseasm(:,:,3));
aut_perc = 100*(change(:,:,4)./yseasm(:,:,4));

all_prec(1,:,:)=win_perc;
all_prec(2,:,:)=spr_perc;
all_prec(3,:,:)=summ_perc;
all_prec(4,:,:)=aut_perc;


%%%%Xing replot this dust figure
ylimcc=[0 20]
xlimcc=[-89 0]
Boundary =shaperead('/home/xzhou473/MATLAB/R2023a/toolbox/map/mapdata/landareas.shp');
CCC0=[-100:2:100];
%CCC1=[-25:1:25];
CCC1=[-80:2:80];
nx=2; %figures you want put in x direction
ny=2; %figures you want put in y direction
gap=[0.05 0.05]; %the gap between each figures
mhorizon=[0.15 0.35]; %whole domain
mvertical=[0.40 0.20]; %whole domain
seasonstring={'Winter (DJF)','Spring (MAM)','Summer (JJA)','Fall (SON)'};
in1=inpolygon(Xd,Yd,lon_contour,lat_contour);
if true
figure(23)
clf
scale2=80;
scale=1; 
set(gcf,'position', [474   157     12*scale2 6*scale2])
set(gcf,'paperposition', [0   0   12*scale  6*scale])
[ha,axw,axh]=adjustsubplot(nx,ny,gap,mvertical,mhorizon)

for i=1:4
axes(ha(i))
contourf(Xd,Yd,squeeze(all_prec(i,:,:)),CCC0,'linecolor','none');hold on
%plot_prec=squeeze(all_prec(i,:,:));
%plot_prec(~in1)=NaN;
%contourf(X,Y,plot_prec,CCC0,'linecolor','none');hold on
plot(lon_contour,lat_contour,'k','LineWidth',2.0);hold on
geoshow(Boundary,'DefaultFaceColor',[0.8 0.8 0.8]);hold on
caxis([CCC1(1) CCC1(end)])
colormap(mycm)
set(gca,'fontsize',10)
ratio=0.3733;%(ylimcc(end)-ylimcc(1))/(xlimcc(end)-xlimcc(1))  %%y/x
          xlim([xlimcc(1) xlimcc(end)])
        ylim([ylimcc(1) ylimcc(end)])
       pbaspect([1 1*ratio 1])
title(seasonstring{i})
set(gca,'TickDir','out')
set(gca,'XMinorTick','on','YMinorTick','on')
if i<3
xticklabels('')
else
xlabel('Longitude (degree)')
end

if mod(i,2)==1
else
yticklabels('')
end

if i==1
ylabel('Latitude (degree)')
ylab=get(gca,'YLabel')
set(ylab, 'Position', get(ylab, 'Position') - [0 15 0]);
end

end

hbar = adjustcolorbar(0.17,0.28,1,22,1.0);
set(hbar,'TickDir','out')
ylabel(hbar,'Δ% Dust')
end
%outfile=['/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/paperfigure/figure_dust/figuredust.png'];print([outfile],'-dpng','-r300');


if false
win_perc = (change(:,:,1));
spr_perc = (change(:,:,2));
summ_perc = (change(:,:,3));
aut_perc = (change(:,:,4));

all_prec(1,:,:)=win_perc;
all_prec(2,:,:)=spr_perc;
all_prec(3,:,:)=summ_perc;
all_prec(4,:,:)=aut_perc;



ylimcc=[0 20]
xlimcc=[-89 0]
Boundary =shaperead('landareas.shp');
CCC0=[-0.02:1e-3:0.02];
%CCC1=[-25:1:25];
CCC1=[-0.02:1e-3:0.02];
nx=2; %figures you want put in x direction
ny=2; %figures you want put in y direction
gap=[0.05 0.05]; %the gap between each figures
mhorizon=[0.15 0.35]; %whole domain
mvertical=[0.40 0.20]; %whole domain
seasonstring={'Winter (a1)','Spring (a2)','Summer (a3)','Fall (a4)'};
in1=inpolygon(Xd,Yd,lon_contour,lat_contour);
figure(24)
clf
scale2=80;
scale=1; 
set(gcf,'position', [474   157     12*scale2 6*scale2])
set(gcf,'paperposition', [0   0   12*scale  6*scale])
[ha,axw,axh]=adjustsubplot(nx,ny,gap,mvertical,mhorizon)

for i=1:4
axes(ha(i))
contourf(Xd,Yd,squeeze(all_prec(i,:,:)),CCC0,'linecolor','none');hold on
%plot_prec=squeeze(all_prec(i,:,:));
%plot_prec(~in1)=NaN;
%contourf(X,Y,plot_prec,CCC0,'linecolor','none');hold on
plot(lon_contour,lat_contour,'k','LineWidth',2.0);hold on
geoshow(Boundary,'DefaultFaceColor',[0.8 0.8 0.8]);hold on
caxis([CCC1(1) CCC1(end)])
colormap(mycm)
set(gca,'fontsize',10)
ratio=0.3733;%(ylimcc(end)-ylimcc(1))/(xlimcc(end)-xlimcc(1))  %%y/x
          xlim([xlimcc(1) xlimcc(end)])
        ylim([ylimcc(1) ylimcc(end)])
       pbaspect([1 1*ratio 1])
title(seasonstring{i})
set(gca,'TickDir','out')
set(gca,'XMinorTick','on','YMinorTick','on')
if i<3
xticklabels('')
else
xlabel('Longitude (degree)')
end

if mod(i,2)==1
else
yticklabels('')
end

if i==1
ylabel('Latitude (degree)')
ylab=get(gca,'YLabel')
set(ylab, 'Position', get(ylab, 'Position') - [0 15 0]);
end

end

hbar = adjustcolorbar(0.17,0.28,1,22,1.0);
set(hbar,'TickDir','out')
ylabel(hbar,'Dust diffirence (m)')
%outfile=['/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/xingplot/deltaMLDabsolutevalue_20112022vs19992011.png'];print([outfile],'-dpng','-r300');
end

%%
CCC0=[-0.1:0.01:0.1];
%CCC1=[-25:1:25];
CCC1=[-0.1:0.01:0.1];
nx=1; %figures you want put in x direction
  ny=1; %figures you want put in y direction
gap=[0.05 0.05]; %the gap between each figures
mhorizon=[0.10 0.10]; %gap between margin and figures in x direction
mvertical=[0.15 0.15]; %gap between margin and figures in y direction
figure(13)
clf
[ha,axw,axh]=adjustsubplot(nx,ny,gap,mvertical,mhorizon)
axes(ha(1))
data1plot=data(:,:,108:end);
contourf(Xd,Yd,squeeze(nanmean(data1plot,3))-nanmean(data1plot(:)),CCC0,'linecolor','none');hold on
%plot_prec=squeeze(all_prec(i,:,:));
%plot_prec(~in1)=NaN;
%contourf(X,Y,plot_prec,CCC0,'linecolor','none');hold on
plot(lon_contour,lat_contour,'k','LineWidth',2.0);hold on
geoshow(Boundary,'DefaultFaceColor',[0.8 0.8 0.8]);hold on
caxis([CCC1(1) CCC1(end)])
colormap(rdbu)
set(gca,'fontsize',13)
ratio=0.3733;%(ylimcc(end)-ylimcc(1))/(xlimcc(end)-xlimcc(1))  %%y/x
          xlim([xlimcc(1) xlimcc(end)])
        ylim([ylimcc(1) ylimcc(end)])
       pbaspect([1 1*ratio 1])
set(gca,'TickDir','out')
set(gca,'XMinorTick','on','YMinorTick','on')
text(-84.0,4.0,['Mean = ' num2str(nanmean(data1plot(:))) ])
text(-84.0,2.0,['STD = 0.0182'])
cb = colorbar;
ylabel(cb,'DAOD anomaly');
%outfile=['/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/paperfigure/figure_nonlinearmodel/sigmaissue/DUST_sigma1.png'];print([outfile],'-dpng','-r300');

%% Select domain and calculate monthly averaged MLD
%%%%lon_contour, lat_contour is the GASB
%%%%lon_big, lat_big is lat 1 to 15

lon_big=[-89,-15,-15,-89,-89]
lat_big=[1,1,15,15,1];

inaa = inpolygon(Xd, Yd, lon_contour, lat_contour);
%inaa = inpolygon(Xd, Yd, lon_big, lat_big);
MLD = zeros(240, 1);

for i = 1:240
    tmp1 = squeeze(data(:,:,i));
    MLD(i) = mean(tmp1(inaa), 'omitnan');
end

MLD_climatology = zeros(12,1);
for month = 1:12
    MLD_climatology(month) = mean(MLD(month:12:end), 'omitnan');
end
MLD_deseasonalized = MLD - repmat(MLD_climatology, 20, 1);

[sigma_detrend,MLD_detrend] = detrend_and_std(MLD_deseasonalized(108:end))

sigma_MC = compute_uncertainty(data(:,:,108:end),inaa)

tsigma=[2003+1/12:1/12:2023];
figure(11)
clf
plot(tsigma,MLD_deseasonalized,'r','linewidth',2.0);hold on
plot(tsigma(108:end),MLD_detrend,'b','linewidth',2.0);hold on
set(gca,'fontsize',13)
legend('Deseasonalized data','Deseaonalized and detrent data')
ylabel('DAOD anomaly (m)')
xlim([2012 2023])
%ylim([-6 8])
text(2019, -0.035,['Mean = ' num2str(nanmean(MLD_detrend)) '']);
text(2019,-0.055, ['STD = 0.0192  '])
%outfile=['/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/paperfigure/figure_nonlinearmodel/sigmaissue/DUST_sigma2.png'];print([outfile],'-dpng','-r300');



;  % 240 = 12 * 20

MLD_annual_mean = NaN(5, 20);  % 2003–2022

for yr = 1:20
    idx_start = (yr - 1) * 12 + 1;
    idx_end   = yr * 12;

    % Annual mean
    MLD_annual_mean(1, yr) = mean(MLD_deseasonalized(idx_start:idx_end), 'omitnan');

    % DJF
    if yr == 1  % 2003: only Jan and Feb
        idxs = [idx_start, idx_start + 1];
    else        % Later years: Dec(prev) + Jan + Feb
        idxs = [idx_start - 1, idx_start, idx_start + 1];
    end
    MLD_annual_mean(2, yr) = mean(MLD_deseasonalized(idxs), 'omitnan');

    % MAM
    MLD_annual_mean(3, yr) = mean(MLD_deseasonalized(idx_start + 2 : idx_start + 4), 'omitnan');

    % JJA
    MLD_annual_mean(4, yr) = mean(MLD_deseasonalized(idx_start + 5 : idx_start + 7), 'omitnan');

    % SON
    MLD_annual_mean(5, yr) = mean(MLD_deseasonalized(idx_start + 8 : idx_start + 10), 'omitnan');
end




%%
years=[2003:1:2022];
titlestring={'Annual mean','Winter (DJF) mean (b1)','Spring (MAM) mean (b2)','Summer (JJA) mean (b3)','Fall (SON) mean (b4)'}
nx=2; %figures you want put in x direction
ny=2; %figures you want put in y direction
gap=[0.05 0.05]; %the gap between each figures
mhorizon=[0.10 0.10]; %whole domain
mvertical=[0.25 0.25]; %whole domain


figure(12)
clf
set_portrait
[ha,axw,axh]=adjustsubplot(nx,ny,gap,mvertical,mhorizon)
for i=2:5
axes(ha(i-1))
plot(years,MLD_annual_mean(i,:),'ko-','linewidth',2.0);hold on
%plot(years(1:12),ones(1,12).*mean(MLD_annual_mean(i,1:12),2),'-','color',[254,178,76]./255,'linewidth',2.0);hold on
%plot(years(13:24),ones(1,12).*mean(MLD_annual_mean(i,13:24),2),'-','color',[189,1,38]./255,'linewidth',2.0);hold on
xlim([1999 2022])
if mod(i,2)==0
ylabel('MLD anomaly (m)')
else
yticklabels('')
end
if i>3
xlabel('Years')
else
xticklabels('')
end
set(gca,'fontsize',9)
set(gca,'Tickdir','out')
title(titlestring{i})
%lgd=legend({['Mean MLD anomaly'], ...
%        ['1993-2010 mean: ' num2str(mean(MLD_annual_mean(i,1:12),2), '%.2f')], ...
%        ['2011-2022 mean: ' num2str(mean(MLD_annual_mean(i,13:24),2), '%.2f')]}, ...
%       'Location', 'northwest');
legend boxoff
%ylim([-1 1])
end


%outfile=['/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/paperfigure/figure_MLD/deltaMLDtimeseries_20112022vs19992011_newdata.png'];print([outfile],'-dpng','-r300');

%%
%%%%%VERSION 2
titlestringV1={'Winter (DJF)','Spring (MAM)','Summer (JJA)','Fall (SON)'};
titlestringV2={'Annual mean','Winter (DJF)','Spring (MAM)','Summer (JJA)','Fall (SON)'};
numstringV2={,'B0','B1','B2','B3','B4'};
numstringV1={'A1','A2','A3','A4'};
nx=4; %figures you want put in x direction
ny=2; %figures you want put in y direction
gap=[0.05 0.15]; %the gap between each figures
mhorizon=[0.10 0.05]; %whole domain
mvertical=[0.15 0.15]; %whole domain
CCC0=[-100:2:100];
%CCC1=[-25:1:25];
CCC1=[-80:2:80];
figure(19)
clf
scale2=80;
scale=1; 
set(gcf,'position', [474   157     12*scale2 12*scale2])
set(gcf,'paperposition', [0   0   12*scale  12*scale])
[ha,axw,axh]=adjustsubplot(nx,ny,gap,mvertical,mhorizon)

for i=2:5
axes(ha(2*i-2))
plot(years,MLD_annual_mean(i,:),'ko-','linewidth',2.0);hold on
plot(years(1:8),ones(1,8).*mean(MLD_annual_mean(i,1:8),2),'-','color',[254,178,76]./255,'linewidth',2.0);hold on
plot(years(9:20),ones(1,12).*mean(MLD_annual_mean(i,9:20),2),'-','color',[189,1,38]./255,'linewidth',2.0);hold on
xlim([2003 2022])


if i==2
ylabel('Dust Aerosol Optical Depth (DAOD) anomaly')
ylab=get(gca,'YLabel')
set(ylab, 'Position', get(ylab, 'Position') - [1 0.23 0]);
end


set(gca,'fontsize',10)
set(gca,'Tickdir','out')
title(titlestringV2{i})
lgd=legend({['Mean DAOD anomaly'], ...
        ['2003-2010 mean: ' num2str(mean(MLD_annual_mean(i,1:8),2), '%.1e')], ...
        ['2011-2022 mean: ' num2str(mean(MLD_annual_mean(i,9:20),2), '%.1e')]}, ...
       'Location', 'northwest');
legend boxoff
ylim([-0.04 0.08])
set(gca,'YMinorTick','on')
if i>4
xlabel('Years')
else
xticklabels('')
end
%text(2020,7.0,numstringV2{i})
end

% percentual changes of yseasm computed as 100*(yseasm_changes/yseam_before)
%win_perc = change(:,:,1);
%spr_perc = change(:,:,2);
%summ_perc = change(:,:,3);
%aut_perc = change(:,:,4);

%all_prec(1,:,:)=win_perc;
%all_prec(2,:,:)=spr_perc;
%all_prec(3,:,:)=summ_perc;
%all_prec(4,:,:)=aut_perc;


in1=inpolygon(Xd,Yd,lon_contour,lat_contour);

for i=1:4
axes(ha(2*i-1))
contourf(Xd,Yd,squeeze(all_prec(i,:,:)),CCC0,'linecolor','none');hold on
%plot_prec=squeeze(all_prec(i,:,:));
%plot_prec(~in1)=NaN;
%contourf(X,Y,plot_prec,CCC0,'linecolor','none');hold on
plot(lon_contour,lat_contour,'k','LineWidth',2.0);hold on
geoshow(Boundary,'DefaultFaceColor',[0.8 0.8 0.8]);hold on
caxis([CCC1(1) CCC1(end)])
colormap(mycm)
set(gca,'fontsize',10)
ratio=0.3733;%(ylimcc(end)-ylimcc(1))/(xlimcc(end)-xlimcc(1))  %%y/x
          xlim([xlimcc(1) xlimcc(end)])
        ylim([ylimcc(1) ylimcc(end)])
       pbaspect([1 1*ratio 1])
title(titlestringV1{i})
set(gca,'TickDir','out')
set(gca,'XMinorTick','on','YMinorTick','on')
if i<4
xticklabels('')
else
xlabel('Longitude (degree)')
end




if i==1
ylabel('Latitude (degree)')
ylab=get(gca,'YLabel')
set(ylab, 'Position', get(ylab, 'Position') - [0 40 0]);
end

%text(-8,19.0,numstringV1{i})

end

hbar = adjustcolorbar(0.47,0.15,2,0.3,6.5);
set(hbar,'TickDir','out')
ylabel(hbar,'Δ% Dust Aerosol Optical Depth')
outfile=['/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/paperfigure/figure_dust/figuredust.png'];print([outfile],'-dpng','-r300');
%outfile=['/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/utils/newdata/domain/Lat1to15/figuredust.png'];print([outfile],'-dpng','-r300');

%%
function sigma_MC = compute_uncertainty(data, inaa)
% COMPUTE_UNCERTAINTY  Spatial-variability based uncertainty for a
% domain-mean, deseasonalized predictor.
%
% Input:
%    data : 3D array (x, y, t)
%
% Output:
%    sigma_MC : scalar uncertainty (STD) for Monte-Carlo noise



%%%% === Extract dimensions ===
[~, ~, nt] = size(data);

%%%% === 1. Compute monthly climatology for each grid cell ===
% data_clim(x, y, 12)
data_clim = zeros(size(data,1), size(data,2), 12);

for m = 1:12
    idx = m:12:nt;  % all months m across years
    data_clim(:,:,m) = mean(data(:,:,idx), 3, 'omitnan');
end


%%%% === 2. Deseasonalize grid-by-grid ===
data_anom = zeros(size(data));

for t = 1:nt
    m = mod(t-1, 12) + 1;   % month index 1–12
    data_anom(:,:,t) = data(:,:,t) - data_clim(:,:,m);
end


%%%% === 3. Compute domain-mean anomaly time series ===
dom_anom = squeeze(mean(mean(data_anom,1,'omitnan'),2,'omitnan'));  % (t)
%display(dom_anom)

%%%% === 4. Compute spatial STD of anomalies each month ===
spatial_std = zeros(nt,1);

for t = 1:nt
    field = data_anom(:,:,t);
    dom_val = dom_anom(t);
    diff_field = field - dom_val;
    spatial_std(t) = std(diff_field(inaa), 'omitnan');
end


%%%% === 5. Estimate effective sample size ===
% Number of grid points
N = size(spatial_std);

% Ocean grid cells are correlated — assume 20:1 reduction (typical)
Neff = size(N,1) / 1%20;


%%%% === 6. Compute uncertainty of domain-mean anomaly ===
sigma_each_month = spatial_std ./ 1;%sqrt(Neff); %%%STD not SEM

sigma_MC = mean(sigma_each_month, 'omitnan');   % final uncertainty 


end

%%
function [sigma,data_detrended] = detrend_and_std(data)
% detrend_and_std
% ------------------------------------
% Input:
%   data : 1D time series (t x 1 or 1 x t)
%
% Output:
%   data_detrended : data with linear trend removed
%   sigma          : standard deviation of detrended data
%
% Example:
%   [x_dt, sig] = detrend_and_std(SST_dom);

    % Ensure column vector
    data = data(:);

    nt = length(data);
    t = (1:nt)';

    % ---- 1. Fit linear trend ----
    p = polyfit(t, data, 1);              % p(1)*t + p(2)

    % ---- 2. Remove trend ----
    trend = polyval(p, t);
    data_detrended = data - trend;

    % ---- 3. Compute STD of detrended data ----
    sigma = std(data_detrended, 'omitnan');

end