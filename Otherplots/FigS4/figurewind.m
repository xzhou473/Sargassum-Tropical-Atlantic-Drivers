clear all

uw=ncread('/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/utils/newdata/Uwind_era5_1993-2022_f64.nc','u');
vw=ncread('/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/utils/newdata/Vwind_era5_1993-2022_f64.nc','v');
%uw=ncread('/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/utils/newdata/ERA5_winduv_1993_2022_xingcrosscheck.nc','u10');
%vw=ncread('/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/utils/newdata/ERA5_winduv_1993_2022_xingcrosscheck.nc','v10');
%gust=ncread('/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/utils/newdata/ERA5_windgust_1993_2022.nc','i10fg');

wind_u=uw(:,:,1:360);
wind_v=vw(:,:,1:360);
aa1=find(wind_u<-10000);
aa2=find(wind_v<-10000);
wind_u(aa1)=NaN;
wind_v(aa2)=NaN;
wind_speed=sqrt(wind_u.^2+wind_v.^2);


ncload('/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/utils/newdata/mlt_025deg_1993_01-2024_11_F64_masked.nc','time');
start_date = datetime(1993,1,1);
time_dates = start_date+seconds(time(1:360)+86401)-seconds(time(1));

input_file = '/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/utils/newdata/yseam_before_wind.nc';
lon = ncread(input_file, 'longitude');
lat = ncread(input_file, 'latitude');
seasons = {'DJF', 'MAM', 'JJA', 'SON'};
for s = 1:4
   switch seasons{s}
        case 'DJF'
            season_idx = find(month(time_dates) == 12 | ...
                              month(time_dates) == 1 | ...
                              month(time_dates) == 2);
        case 'MAM'
            season_idx = find(month(time_dates) == 3 | ...
                              month(time_dates) == 4 | ...
                              month(time_dates) == 5);
        case 'JJA'
            season_idx = find(month(time_dates) == 6 | ...
                              month(time_dates) == 7 | ...
                              month(time_dates) == 8);
        case 'SON'
            season_idx = find(month(time_dates) == 9 | ...
                              month(time_dates) == 10 | ...
                              month(time_dates) == 11);
   end
   season_allidx(:,s)=season_idx;
    %season_idx = season_idx(season_idx >= start_idx & season_idx <= end_idx);
end

aveindex1=[1 3:3:88 90];
aveindex2=[1:3:91];

for s=1:4
 if s==1
 aveindex=aveindex1;
 else
 aveindex=aveindex2;
 end
 for i=1:30
 wind_uall(i,s,:,:)=squeeze(mean(wind_u(:,:,season_allidx([aveindex(i):aveindex(i+1)-1],s)),3));
 wind_vall(i,s,:,:)=squeeze(mean(wind_v(:,:,season_allidx([aveindex(i):aveindex(i+1)-1],s)),3));
 wind_speedall(i,s,:,:)=squeeze(mean(wind_speed(:,:,season_allidx([aveindex(i):aveindex(i+1)-1],s)),3));
 end
end


mycm = table2array(readtable('/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/scripts/test_mycm.csv'));
TAB = xlsread('/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/utils/digitized_contour/Winter_Map_Jouanno_et_al_2023_DigitizedManual.xlsx');

%change=ncread('/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/utils/newdata/changes_yseasm_wind_after-before.nc','u');
%yseasm=ncread('/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/utils/newdata/yseam_before_wind.nc','u');

lon_contour = TAB(:,1); 
lat_contour = TAB(:,2); 

[Xd Yd] = meshgrid(lon,lat);
Xd = Xd';
Yd = Yd';
inaa = inpolygon(Xd,Yd,lon_contour,lat_contour);

%%


ylimcc=[0 20]
xlimcc=[-89 0]
Boundary =shaperead('landareas.shp');
CCC0=[-3:0.1:3];
%CCC1=[-25:1:25];
CCC1=[-1:0.1:1];
nx=2; %figures you want put in x direction
ny=2; %figures you want put in y direction
gap=[0.05 0.05]; %the gap between each figures
mhorizon=[0.15 0.35]; %whole domain
mvertical=[0.40 0.20]; %whole domain
seasonstring={'Winter (a)','Spring (b)','Summer (c)','Fall (d)'};
in1=inpolygon(Xd,Yd,lon_contour,lat_contour);
figure(23)
clf
scale2=80;
scale=1; 
set(gcf,'position', [474   157     12*scale2 6*scale2])
set(gcf,'paperposition', [0   0   12*scale  6*scale])
[ha,axw,axh]=adjustsubplot(nx,ny,gap,mvertical,mhorizon)

INT_INDEX=8;
mult=5;%0.05;%0.05;
scale=0;

for i=1:4

uplot=squeeze(mean(wind_uall(19:30,i,:,:),1))-squeeze(mean(wind_uall(7:18,i,:,:),1));
vplot=squeeze(mean(wind_vall(19:30,i,:,:),1))-squeeze(mean(wind_vall(7:18,i,:,:),1));
speedplot=squeeze(mean(wind_speedall(19:30,i,:,:),1))-squeeze(mean(wind_speedall(7:18,i,:,:),1));

axes(ha(i))
contourf(Xd,Yd,speedplot,CCC0,'linecolor','none');hold on
%plot_prec=squeeze(all_prec(i,:,:));
%plot_prec(~in1)=NaN;
%contourf(X,Y,plot_prec,CCC0,'linecolor','none');hold on
plot(lon_contour,lat_contour,'k','LineWidth',2.0);hold on
geoshow(Boundary,'DefaultFaceColor',[0.8 0.8 0.8]);hold on
quiver(Xd(1:INT_INDEX:end,1:INT_INDEX:end),Yd(1:INT_INDEX:end,1:INT_INDEX:end),uplot(1:INT_INDEX:end,1:INT_INDEX:end).*mult,vplot(1:INT_INDEX:end,1:INT_INDEX:end).*mult,scale,'k');hold on
quiver(-10,15,0,0.1.*mult,scale,'k'); hold on
%text(-15,13,'0.1m/s')


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
ylabel(hbar,'\delta% wind')

%outfile=['/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/xingplot/windtest/windgust_spatial.png'];print([outfile],'-dpng','-r300');

%%
data=wind_speed;
for i=1:360%size(data,3)
tmp1=squeeze(data(:,:,i));
datatmp1(i)=mean(tmp1(inaa), 'omitnan');
end


datatmp1_climatology = zeros(12,1);
for month = 1:12
    datatmp1_climatology(month) = mean(datatmp1(month:12:end), 'omitnan');
end
datatmp1_deseasonalized = datatmp1' - repmat(datatmp1_climatology, length(datatmp1)/12, 1);
clear MLD_annual_mean
for i = 1:24
    sidx = (i-1)*12 + 73;
    eidx = (i)*12 + 73-1;
   datatmp1_annual_mean(1,i) = mean(datatmp1_deseasonalized(sidx:eidx), 'omitnan'); % Annual mean
end

for i = 1:24
    sidx = (i-1)*12 + 72;
    eidx = (i-1)*12 + 74;
    datatmp1_annual_mean(2,i) = mean(datatmp1_deseasonalized(sidx:eidx), 'omitnan'); % Annual mean
end

for i = 1:24
    sidx = (i-1)*12 + 75;
    eidx = (i-1)*12 + 77;
    datatmp1_annual_mean(3,i) = mean(datatmp1_deseasonalized(sidx:eidx), 'omitnan'); % Annual mean
end

for i = 1:24
    sidx = (i-1)*12 + 78;
    eidx = (i-1)*12 + 80;
   datatmp1_annual_mean(4,i) = mean(datatmp1_deseasonalized(sidx:eidx), 'omitnan'); % Annual mean
end


for i = 1:24
    sidx = (i-1)*12 + 81;
    eidx = (i-1)*12 + 83;
   datatmp1_annual_mean(5,i) = mean(datatmp1_deseasonalized(sidx:eidx), 'omitnan'); % Annual mean
end

years=[1999:1:2022];
titlestring={'Annual mean','Winter(DJF) mean','Spring(MAM) mean','Summer(JJA) mean','Fall(SON) mean'}
figure(12)
clf
set_portrait
for i=2:5
subplot(3,2,i-1)
%plot(years,datatmp1_annual_mean(i,:),'bo-','linewidth',2.0);hold on
plot(years([1:1:12]),datatmp1_annual_mean(i,[1:1:12]),'bo-','linewidth',2.0);hold on
plot(years([13:24]),datatmp1_annual_mean(i,[13:24]),'bo-','linewidth',2.0);hold on
%plot(years(1:12),ones(1,12).*mean(datatmp1_annual_mean(i,1:12),2),'-','color',[254,178,76]./255,'linewidth',2.0);hold on
plot(years(1:12),ones(1,12).*mean(datatmp1_annual_mean(i,1:12),2),'-','color',[254,178,76]./255,'linewidth',2.0);hold on
plot(years(13:24),ones(1,12).*mean(datatmp1_annual_mean(i,13:24),2),'-','color',[189,1,38]./255,'linewidth',2.0);hold on
xlim([1999 2022])
ylabel('Wind speed anomaly (m/s)')
xlabel('Years')
set(gca,'fontsize',9)
set(gca,'Tickdir','out')
title(titlestring{i})
lgd=legend({['Mean Wind anomaly'], ...
        ['1993-2010 mean: ' num2str(mean(datatmp1_annual_mean(i,1:12),2), '%.2f')], ...
        ['2011-2022 mean: ' num2str(mean(datatmp1_annual_mean(i,13:24),2), '%.2f')]}, ...
       'Location', 'northwest');
legend boxoff
ylim([-1 1.5])
end

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
CCC0=[-3:0.1:3];
%CCC1=[-25:1:25];
CCC1=[-1:0.1:1];
figure(19)
clf
scale2=80;
scale=1; 
set(gcf,'position', [474   157     12*scale2 12*scale2])
set(gcf,'paperposition', [0   0   12*scale  12*scale])
[ha,axw,axh]=adjustsubplot(nx,ny,gap,mvertical,mhorizon)

for i=2:5
axes(ha(2*i-2))
plot(years,datatmp1_annual_mean(i,:),'ko-','linewidth',2.0);hold on
plot(years(1:12),ones(1,12).*mean(datatmp1_annual_mean(i,1:12),2),'-','color',[254,178,76]./255,'linewidth',2.0);hold on
plot(years(13:24),ones(1,12).*mean(datatmp1_annual_mean(i,13:24),2),'-','color',[189,1,38]./255,'linewidth',2.0);hold on
xlim([1999 2022])


if i==2
ylabel('Wind speed anomaly (m/s)')
ylab=get(gca,'YLabel')
set(ylab, 'Position', get(ylab, 'Position') - [1 5 0]);
end

if i>4
xlabel('Years')
else
xticklabels('')
end
set(gca,'fontsize',10)
set(gca,'Tickdir','out')
title(titlestringV2{i})
lgd=legend({['Wind speed anomaly (m/s)'], ...
        ['1993-2010 mean: ' num2str(mean(datatmp1_annual_mean(i,1:12),2), '%.2f')], ...
        ['2011-2022 mean: ' num2str(mean(datatmp1_annual_mean(i,13:24),2), '%.2f')]}, ...
       'Location', 'northwest');
legend boxoff
ylim([-1 1.5])
%text(2020,7.0,numstringV2{i})
end

% percentual changes of yseasm computed as 100*(yseasm_changes/yseam_before)
%win_perc = 100*(change(:,:,1)./yseasm(:,:,1));
%spr_perc = 100*(change(:,:,2)./yseasm(:,:,2));
%summ_perc = 100*(change(:,:,3)./yseasm(:,:,3));
%aut_perc = 100*(change(:,:,4)./yseasm(:,:,4));

%all_prec(1,:,:)=win_perc;
%all_prec(2,:,:)=spr_perc;
%all_prec(3,:,:)=summ_perc;
%all_prec(4,:,:)=aut_perc;


in1=inpolygon(Xd,Yd,lon_contour,lat_contour);

for i=1:4
uplot=squeeze(mean(wind_uall(19:30,i,:,:),1))-squeeze(mean(wind_uall(7:18,i,:,:),1));
vplot=squeeze(mean(wind_vall(19:30,i,:,:),1))-squeeze(mean(wind_vall(7:18,i,:,:),1));
speedplot=squeeze(mean(wind_speedall(19:30,i,:,:),1))-squeeze(mean(wind_speedall(7:18,i,:,:),1));


axes(ha(2*i-1))
contourf(Xd,Yd,speedplot,CCC0,'linecolor','none');hold on
%plot_prec=squeeze(all_prec(i,:,:));
%plot_prec(~in1)=NaN;
%contourf(X,Y,plot_prec,CCC0,'linecolor','none');hold on
plot(lon_contour,lat_contour,'k','LineWidth',2.0);hold on
geoshow(Boundary,'DefaultFaceColor',[0.8 0.8 0.8]);hold on
quiver(Xd(1:INT_INDEX:end,1:INT_INDEX:end),Yd(1:INT_INDEX:end,1:INT_INDEX:end),uplot(1:INT_INDEX:end,1:INT_INDEX:end).*mult,vplot(1:INT_INDEX:end,1:INT_INDEX:end).*mult,scale,'k');hold on
quiver(-10,15,0,0.1.*mult,scale,'k'); hold on
%text(-15,13,'0.1m/s')


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
ylabel(hbar,'Wind speed difference (m/s)')

%outfile=['/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/paperfigure/figure_wind/figurewind.png']; print([outfile],'-dpng','-r300');