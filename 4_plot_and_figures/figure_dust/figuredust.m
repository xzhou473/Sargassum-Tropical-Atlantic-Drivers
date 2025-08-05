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
Boundary =shaperead('landareas.shp');
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
%% Select domain and calculate monthly averaged MLD
%%%%lon_contour, lat_contour is the GASB
%%%%lon_big, lat_big is lat 1 to 15

lon_big=[-89,-15,-15,-89,-89]
lat_big=[1,1,15,15,1];

%inaa = inpolygon(Xd, Yd, lon_contour, lat_contour);
inaa = inpolygon(Xd, Yd, lon_big, lat_big);
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
numstringV2={,'(b0)','(b1)','(b2)','(b3)','(b4)'};
numstringV1={'(a1)','(a2)','(a3)','(a4)'};
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
        ['1993-2010 mean: ' num2str(mean(MLD_annual_mean(i,1:8),2), '%.1e')], ...
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
%outfile=['/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/paperfigure/figure_dust/figuredust.png'];print([outfile],'-dpng','-r300');
outfile=['/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/utils/newdata/domain/Lat1to15/figuredust.png'];print([outfile],'-dpng','-r300');