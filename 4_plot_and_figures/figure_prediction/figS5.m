% ----------------------------------------
% Written by ljubanovi@yahoo.it
% ----------------------------------------

close all
clear all
clc

%addpath(genpath('../../../Matlab_Additional_Toolboxes/')) % if additional packages are needed

close all
clc 


%% NOTE: the MLT timeseries are obtained by the 1993-2022 data at original resolution (1/12deg).
%  When data are deseasonalized, the seasonal cycle was removed point-by-point in
%  the original high-res data PRIOR to compute the field average.
%  Deseasonalization and fldmeans are computed a priori using CDO.

% The fldmean are over the REDUCED lonlatbox = -89,-15,1,11 (or an ever different area, see script cdo: 1_mlt_preprocess.sh)

% Path to mlt timeseries: 
cd ('/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/scripts')
%resultspath1 = '../utils/fldmean_data/fldm_of_1993_2022_glorys_mlt_reducedGASB_masked_deseas_1-12deg.nc'     % fldm of deseasonalized mlt (see readme in that folder)
%resultspath1 = '/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/utils/newdata/fldm_of_1993_2022_glorys_mlt_reducedGASB_masked_deseas_1-12deg.nc' 
resultspath1 = '..//utils/newdata/domain/Lat1to15/fldm_of_1993_dec2024_glorys_mlt_masked_deseas_1-12deg_Lat_1-15.nc' 
%%%%% spatial field 1993_2022_glorys_mlt_TA_masked_025remapped_GASB_deseas.nc                                      

% Path to dust timeseries:
%resultspath3 = '../utils/fldmean_data/fldm_of_duaod550_f64_res075_2003-2022_deseas_reducedGASB.nc'                % fldmean of deseasonalized duaod550 (see readme in that folder)
resultspath3 = '/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/utils/newdata/fldm_of_duaod550_f64_res075_2003-2022_deseas_reducedGASB.nc' 
%%%% spatial field duaod550_f64_res075_2003-2022_GASB_deseas.nc

% -----------------
% MLT timeseries:
% -----------------

% Deseasonalized
ncid = netcdf.open(resultspath1)

  time = double(netcdf.getVar(ncid,0)); 
  mlt_flmd_deseas = double(netcdf.getVar(ncid,3));  
  netcdf.inqVar(ncid,3)
netcdf.close(ncid);

YEAR_INI = 1993                                     % Initial tear for mlt data (glorys rean + analysis)
datevector = datetime(YEAR_INI,1:length(time),15);  % Date vector for mlt data  



% -----------------
% Dust timeseries:
% -----------------

% Deseas
ncid = netcdf.open(resultspath3)
  duaod550_flmd_deseas = double(netcdf.getVar(ncid,3));  % This one
  netcdf.inqVar(ncid,3)
netcdf.close(ncid);


YEAR_INI_dust = 2003 
datevector_dust = datetime(YEAR_INI_dust,1:length(time),15);


%% PART TWO: Load Web-Digitized observations: MUST CITE DATA SOURCE AND SOFTWARE USED.
% ------------------------
%  ORIGINAL DATA SOURCE :
% ------------------------
%  These data have been obtained by manually extracting points from 
%  Fig. 1 (e) in Jouanno et al 2023 https://doi.org/10.1029/2023GL105545
%
% ------------------
%  DIGITIZER TOOL: 
% ------------------
%  To extract data from their figure, I used the online version of the free software 
%  WebPlot Digitizer  (https://automeris.io/WebPlotDigitizer.html). 
%  It is an open source computer vision assisted software that extracts data from plots. 
%  I didndt use the automatic extraction, but the manual extraction.
%  (Rohatgi, 2024, https://automeris.io/WebPlotDigitizer.html, WebPlotDigitizer version 4.7, 2024)
% 
%  The data are digitized starting from Jan 2011, i.e. all the 2010 points are not digitized 
%  because the Sargassum bloom conventionally started in year 2011. So that s why we have less points than Fig 1(e) in the original paper (they start from jan 2010 which for me is not important).

% ------------------------------------------------------------------------
% The data I want are stored in the second column of the .csv file: 
% ------------------------------------------------------------------------

TAB = csvread('/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/utils/digitized_sarg_biomass/WebDigitized_Jan2011-Dec2022_from_Jouanno-et-al-2023.csv');

obs_webdig = TAB(:,2); 

ncload('/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/paperfigure/figure_prediction/fldsum_biomass_MTon.nc','Mean_FC');
obs_webdig=[obs_webdig; Mean_FC];
% It's impossible to have negative concentration. Manual digitizing might
% introduce spurious negative values when concentrations are very small. 
% Hence, put to zero every negative point.

obs_webdig(obs_webdig<0) = 0;


% Create date vector for digitized observations:  Jan 2011 to Dec 2022
observation_dates = datetime(2011,1:length(obs_webdig),15)



%% TEMPORAL SUBSET OF MLT AND DUST TIMESERIES: 
% Select, in the glorys data vector, values on timepoints that overlap with observations: 
% i.e. subset the glorys dataset from Jan 2011 to Dec 2020.


%%%%xing comment, too ensure the exactly same, need calculate individually
%%%%for 144, and 168
% ====================
% 1. Focus on mlt 
% ====================

% =============================
% 1.1. Subset mlt timeseries: 
% =============================

a = 18*12+1             % index where I have Jan 2011 in the glorys date vector. Please always crosscheck.
b = length(datevector); % because Dec 2020 is already the last timestep for the glorys dataset loaded here. Pls always crosscheck. 
b1=360;
b2=length(datevector);


% Crosscheck: must be from Jan 2011 to Dec 2020: 
datevector(a)   % This must be Jan 2011
datevector(b)   % This must be Dec 2022

% Check on Jan 2011
if datevector(a) == datetime(2011,1,15)
   disp '' 
else 
   error 'datevector(a) must be Jan 2011 but it is not '
end 

% Check on Dec 2022
%if datevector(b) == datetime(2022,12,15)
%   disp '' 
%else 
%   error 'datevector(a) must be Dec 2022 but it is not '
%end 

% ========================================
% 1.2. Find mlt lag for regression model: 
% ========================================

% Find a lag within a +/- 12 months lag range that maximises the |c.c.| computed between: (Sarg Biomass - mean(SargBiomass)) and (smoothed MLT - mean(smoothed MLT)).
% Then I will use this lag in the regression equation below. 

 LAG = 6;    % Focus on +/- 6 months lag range only

% -----------------------------------------------------------
% 1.2.1 Apply a N-month moving mean to mlt to remove some noise  
% -----------------------------------------------------------

Window_for_MLT_timeser = 5;%5;     % N-month moving mean. 5

%%%%xingtestusing GASB wind instead
%load('/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/utils/wind/wind_regionserious.mat')
%MLT = (movmean(wind_deseasonalized(a:b),Window_for_MLT_timeser)); % Movmean mlt signal, less noise.
MLT1 = (movmean(squeeze(mlt_flmd_deseas(:,:,a:b1)),Window_for_MLT_timeser)); % Movmean mlt signal, less noise.
MLT2 = (movmean(squeeze(mlt_flmd_deseas(:,:,a:b2)),Window_for_MLT_timeser));

% --------------------------------------------------------------------------------------------------------------------------------------
% 1.2.2 Compute here the lagged correlation (xcorr) between :(Sarg Biomass - mean(SargBiomass)) and (smoothed MLT - mean(smoothed MLT)): 
% --------------------------------------------------------------------------------------------------------------------------------------

% Attenzione! Se usi xcorr devi togliere le medie, se i segnali non sono a media nulla, 
% per avere a zero-lag lo stesso risultato di corcoeff o corr (cioè il Pearson correlation coeff). 
MLT=MLT1;
[c,lags] = xcorr(obs_webdig(1:144)-nanmean(obs_webdig(1:144)),MLT-nanmean(MLT),'normalized')
figure 
s = stem(lags,c)
set(gca,'FontSize',18,'FontName','Times');
xlabel('months','Interpreter','latex');
ylabel('c.c','Interpreter','latex');
title ('xcorr(sarg,mlt)','FontSize',17,'FontName','Times');
xlim([-LAG LAG])

%print(gcf,'../Figures/XCorr_sarg_mld.png','-dpng','-r350')


% -------------------------------------------------------------------------------------------
%  1.2.3 Find which is the lag in the stem plot above, for which we have the max|c.c.|, 
% (limiting your search within +/- 12 months lag range)
% -------------------------------------------------------------------------------------------

% only consider the lag range of your choice: 
idx_from = find(lags == -LAG)
idx_to = find(lags == LAG)

lag_range = lags(idx_from:idx_to);
cc_range = c(idx_from:idx_to); 

% Find the max |c.c.| in this lag range and the lag at which it happens:

max_abs_cc = max(abs(cc_range))                 % the max abs(cc) in teh selected lag-range

idx_max_cc = find(abs(cc_range) == max_abs_cc); % where it happens, i.e. which index position in the lag and cc subset 

lag_at_max_abs_cc = lag_range(idx_max_cc)       % what is the corresponding lag --> To put in the following calculation

% ---------------------------------------------------------------------------------------------------------------
% 1.2.4 Then I shift the mlt timeseries of N months, such that the |c.c.| is max (i.e. lag as just computed)
% ---------------------------------------------------------------------------------------------------------------

% How do I know I have to shift mlt in the past? Because I used the xcorr(a,b)
% function and the resulting lag at which |cc| is max was positive. And for
% this function, xcorr(a,b) positive lag means b precedes a. 

% Do I shift in the past or in the future? It depends on the sign of the
% lag because of the way xcorr(a,b) works (lag > 0 means b precedes a.)

if lag_at_max_abs_cc > 0 
   MESI_mlt = lag_at_max_abs_cc             % if lag > 0 then shift in the past.
else 
   MESI_mlt = -lag_at_max_abs_cc            % else shift in the future (by simply changing that sign)
end 

MLT_shifted=circshift(MLT,MESI_mlt);

MLT_shifted1=MLT_shifted;

MLT=MLT2;
[c,lags] = xcorr(obs_webdig-nanmean(obs_webdig),MLT-nanmean(MLT),'normalized')
figure 
s = stem(lags,c)
set(gca,'FontSize',18,'FontName','Times');
xlabel('months','Interpreter','latex');
ylabel('c.c','Interpreter','latex');
title ('xcorr(sarg,mlt)','FontSize',17,'FontName','Times');
xlim([-LAG LAG])

%print(gcf,'../Figures/XCorr_sarg_mld.png','-dpng','-r350')


% -------------------------------------------------------------------------------------------
%  1.2.3 Find which is the lag in the stem plot above, for which we have the max|c.c.|, 
% (limiting your search within +/- 12 months lag range)
% -------------------------------------------------------------------------------------------

% only consider the lag range of your choice: 
idx_from = find(lags == -LAG)
idx_to = find(lags == LAG)

lag_range = lags(idx_from:idx_to);
cc_range = c(idx_from:idx_to); 

% Find the max |c.c.| in this lag range and the lag at which it happens:

max_abs_cc = max(abs(cc_range))                 % the max abs(cc) in teh selected lag-range

idx_max_cc = find(abs(cc_range) == max_abs_cc); % where it happens, i.e. which index position in the lag and cc subset 

lag_at_max_abs_cc = lag_range(idx_max_cc)       % what is the corresponding lag --> To put in the following calculation

% ---------------------------------------------------------------------------------------------------------------
% 1.2.4 Then I shift the mlt timeseries of N months, such that the |c.c.| is max (i.e. lag as just computed)
% ---------------------------------------------------------------------------------------------------------------

% How do I know I have to shift mlt in the past? Because I used the xcorr(a,b)
% function and the resulting lag at which |cc| is max was positive. And for
% this function, xcorr(a,b) positive lag means b precedes a. 

% Do I shift in the past or in the future? It depends on the sign of the
% lag because of the way xcorr(a,b) works (lag > 0 means b precedes a.)

if lag_at_max_abs_cc > 0 
   MESI_mlt = lag_at_max_abs_cc             % if lag > 0 then shift in the past.
else 
   MESI_mlt = -lag_at_max_abs_cc            % else shift in the future (by simply changing that sign)
end 

MLT_shifted=circshift(MLT,MESI_mlt);

MLT_shifted2=MLT_shifted;
MLT_shifted(1:144)=MLT_shifted1;

%MLT_shifted=1./(1+exp(-1*MLT_shifted));
% Initialize an array of zeros with size 144

%MLT_shifted=f_MLT;
%MLT_shifted=f_MLT;

%MLT_shifted_winter = circshift(MLT, 6)
%aa1=find(MLT_shifted_winter==0)
%MLT_shifted_winter(aa1)=9999;

%%                                             
% ==================
% 2. Focus on Dust 
% ==================

% Select in the dust data vector, values on timepoints that overlap with
% observations: i.e. subset the dust dataset from Jan 2011 to Dec 2022. 
a_dust = 97           % index where I have Jan 2011 in the dust date vetcor. 
b_dust = 240
% faccio  la riprova: deve venire da Jan 2011 a Dec 2020 
datevector_dust(a_dust)   % This must be Jan 2011
datevector_dust(b_dust)   % This must be Dec 2022


% Check on Jan 2011
if datevector_dust(a_dust) == datetime(2011,1,15)
   disp '' 
else 
    error 'datevector_dust(a) must be Jan 2011 but it is not '
end 

% Check on Dec 2022
if datevector_dust(b_dust) == datetime(2022,12,15)
   disp '' 
else 
    error 'datevector_dust(a) must be Dec 2022 but it is not '
end 


% -----------------------------------------------------------
% 2.1 Apply a X-month moving sum to dust
% -----------------------------------------------------------

DUST_MOV_WINDOW = 3;        % X-month for moving sum of dust

DUST = movsum(squeeze(duaod550_flmd_deseas(a_dust:b_dust)),[DUST_MOV_WINDOW 0]);  % X-month moving sum for the dust contribution

%DUST=(DUST-min(DUST))./(max(DUST)-min(DUST));
%DUST=(DUST-mean(DUST))./(std(DUST));

%%
% 3. Focus on SST 
%%%GLORY MISSING DATA FROM 2021-06 to 2022-06, temporally using previous year data here
%load('/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/scripts/SST_fornonlinearmodel.mat');
SSTmean=ncread('/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/utils/newdata/fldm_of_1993_2022_glorys_sst_reducedGASB_masked_deseas_1-12deg.nc','thetao')

SSTmean1=squeeze(SSTmean);

%SSTTIME1=zeros(360,1);
%SSTTIME1(1:341)=SSTTIME(1:341);
%SSTTIME1(342:353)=SSTTIME(330:341)+365.*86400;
%SSTTIME1(354:360)=SSTTIME(342:end);

%for tt=1:length(SSTTIME1)
%datevector_SST(tt)=greg2calendar(SSTTIME1(tt)+14.5.*86400,[1970 1 1]);
%end
datevector_SST=datevector ;


a_SST = 217;            % index where I have Jan 2011 in the dust date vetcor. 
b_SST = 360;
% faccio  la riprova: deve venire da Jan 2011 a Dec 2020 
datevector_SST(a_SST)   % This must be Jan 2011
datevector_SST(b_SST)   % This must be Dec 2022


% Check on Jan 2011
if datevector_SST(a_SST) == datetime(2011, 1, 15)
   disp '' 
else 
    error 'datevector_dust(a) must be Jan 2011 but it is not '
end 

% Check on Dec 2022
if datevector_SST(b_SST) == datetime(2022, 12, 15)
   disp '' 
else 
    error 'datevector_dust(a) must be Dec 2022 but it is not '
end 


% -----------------------------------------------------------
% 2.1 Apply a X-month moving sum to dust
% -----------------------------------------------------------

SST_MOV_WINDOW = 3;        % X-month for moving sum of dust

SST_nonlinear = movmean(squeeze(SSTmean1(a_SST:b_SST)),SST_MOV_WINDOW);  

%SST_nonlinear=(SST_nonlinear-min(SST_nonlinear))./(max(SST_nonlinear)-min(SST_nonlinear));

if false
CCCtest=[255,255,212
254,217,142
254,153,41
217,95,14
153,52,4]./255

figure(23)

for i=1:5
SST_MOV_WINDOW = i;        % X-month for moving sum of dust
SST_nonlinear = movmean(squeeze(SSTmean1(a_SST:b_SST)),SST_MOV_WINDOW);  
plot(SST_nonlinear,'linewidth',1.5,'color',CCCtest(i,:));hold on
end
plot(squeeze(SSTmean1(a_SST:b_SST)),'k','linewidth',1.0);hold on
set(gca,'fontsize',10)
xticks([1:60:144])
xticklabels([2011 2016 2021]);
set(gca,'TickDir','out')
ylabel('SST anomaly (\circ C)')
legend('Window=1', 'Window=2', 'Window=3', 'Window=4', 'Window=5','Raw data', 'Location', 'northwest', 'NumColumns', 2);
%outfile=['/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/xingplot/SSTmoveave.png'];print(outfile,'-dpng','-r300');
end

%% NON LINEAR REGRESSION MODEL: 

% -------------------------------------------
%  Try to Fit Non-linear regression model:
% -------------------------------------------

% Define table for variables and observations:
 
% Concentration of Sargassum in current month: i.e. at each month (timestep) the existing bloom (sargassum biomass) 

conc_current_year = obs_webdig; % THIS IS JUST BECAUSE I WANT TO CHANGE THE variable NAME... 
%conc_current_year = (conc_current_year-min(conc_current_year))./(max(conc_current_year)-min(conc_current_year));

% ====================================================================== 
%% Previous concentration of Sargassum, self-seeding effect:
% In each month, I use the average seasonal biomass 
% of the corresponding season in the previous year
% ======================================================================

% First, compute the seasonal averages:
clear x
%x(1:b_dust-a_dust+1) = NaN;
x(1:168) = NaN;
for k = 1:(length(x)/3)-1
    observation_dates(3*k:(3*k)+2);
    x(3*k:(3*k)+2) = mean(obs_webdig(3*k:(3*k)+2));
end
x(1:2) = mean(obs_webdig(1:2));
x(end) = mean(obs_webdig(end));
x = x';

% Then, use the values of 12 months earlier:
% =====================================
    CONC_LAG = 12;%12;                       % Number of months in the past for Sarg: I want to include the concentraction of "CONC_LAG" months before, in the current time step.
% =====================================

conc_previous = circshift(x,CONC_LAG); 

% I add a row of 0s as first indices because if I start in Jan 2011, the preceding year was 2010 
% and there was no-bloom. I need it because I want to have vectors of the same lenght starting
% from Jan 2011. So I have to put something before jan 2011 for the past concentration, i.e. no bloom = 0.  

if CONC_LAG>0
    conc_previous(1:CONC_LAG) = 0; 
else
    conc_previous(end+CONC_LAG+1:end) = 0;
end
remaining_conc = zeros(144, 1);
MLT_avg=zeros(1,144);
%MLT_winter_avg = nan(1, 144);
% Loop through each month (i from 1 to 144)
%for i = 1:length(MLT)
    % Determine the year and month for the current index
 %   year = floor((i-1)/12) + 2011;
 %   month = mod(i-1, 12) + 1;
    
    % Find indices for winter months based on the specified logic
 %   if month == 1  % January
 %       indices = [i-2, i-1, i];  % Nov, Dec of previous year and Jan of current year
 %   elseif month >= 2 && month <= 10  % February to October
 %       indices = [i-2, i-1, i, i+1];  % Nov, Dec of previous year and Jan, Feb of current year
 %   elseif month == 11  % November
 %       indices = [i-1, i, i-11, i-10];  % Dec of previous year and Jan, Feb, Nov of current year
 %   elseif month == 12  % December
 %       indices = [i, i-11, i-10, i-1];  % Jan, Feb, Nov, Dec of current year
 %   end
    
    % Remove indices that are out of bounds (less than 1 or greater than 144)
  %  indices = indices(indices >= 1 & indices <= 144);
    
    % Calculate the average for the selected winter months
   % MLT_winter_avg(i) = mean(MLT(indices));
%end
%%
clear Rall
Rall(1:30)=NaN;
%for toty=[6:1:15 17 18 20:1:30]
toty=12;

% CONC_LAG = toty;
%conc_previous = circshift(x,CONC_LAG);   
%if CONC_LAG>0
%    conc_previous(1:CONC_LAG) = 0; 
%else
%    conc_previous(end+CONC_LAG+1:end) = 0;
%end


% First, compute the seasonal averages:
clear MLT_seasonal1
MLT_seasonal1(1:168) = NaN;
%MLT_seasonal1(1:b_dust-a_dust+1) = NaN;
for k = 1:(length(MLT_seasonal1)/3)-1
    MLT_seasonal1(3*k:(3*k)+2) = mean(MLT2(3*k:(3*k)+2));
end
MLT_seasonal1(1:2) = mean(MLT2(1:2));
MLT_seasonal1(end) = mean(MLT2(end));
MLT_seasonal1 = MLT_seasonal1';
MLT_avg1 = circshift(MLT_seasonal1,toty); 


% First, compute the seasonal averages:
clear MLT_seasonal1
MLT_seasonal1(1:144) = NaN;
%MLT_seasonal1(1:b_dust-a_dust+1) = NaN;
for k = 1:(length(MLT_seasonal1)/3)-1
    MLT_seasonal1(3*k:(3*k)+2) = mean(MLT1(3*k:(3*k)+2));
end
MLT_seasonal1(1:2) = mean(MLT1(1:2));
MLT_seasonal1(end) = mean(MLT1(end));
MLT_seasonal1 = MLT_seasonal1';
MLT_avg2 = circshift(MLT_seasonal1,toty); 
MLT_avg=[MLT_avg2(1:144); MLT_avg1(145:end)];


for t = toty+1:144
    % Compute the 12-month average MLT anomaly and 12-month average conc
    %MLT_avg(t) = mean(MLT(t-toty:t));   
    remaining_conc(t) = mean(conc_current_year(t-toty:t));
end

NORMINDEX=144;%144;
NORMINDEX2=length(MLT_shifted)-NORMINDEX;
% MLT_avg=MLT_winter_avg;
 MLT_avgnorm1=(MLT_avg-min(MLT_avg(1:NORMINDEX)))./(max(MLT_avg(1:NORMINDEX))-min(MLT_avg(1:NORMINDEX)));
 %MLT_avgnorm1=ones(144,1);
 MLT_avgnorm0=(MLT_shifted-min(MLT_shifted(1:NORMINDEX)))./(max(MLT_shifted(1:NORMINDEX))-min(MLT_shifted(1:NORMINDEX)));

 MLT_avgnorm_test1=(MLT_avg-min(MLT_avg))./(max(MLT_avg)-min(MLT_avg));
 %MLT_avgnorm1=ones(144,1);
 MLT_avgnorm_test0=(MLT_shifted-min(MLT_shifted))./(max(MLT_shifted)-min(MLT_shifted));

 DUST_avgnorm=(DUST-min(DUST))./(max(DUST)-min(DUST));
 SST_avgnorm=(SST_nonlinear-min(SST_nonlinear))./(max(SST_nonlinear)-min(SST_nonlinear));
 if size(remaining_conc,2)==144
 remaining_conc=remaining_conc';
 end
alpha = 2.0;



conc_previous = [conc_previous zeros(1,length(MLT_shifted)-length(obs_webdig))]; % Se la plottassi vedrei esattamente cosa voglio. 
%conc_previous=remaining_conc';
%conc_previous=remaining_conc;
%conc_previous=(log10(conc_previous+0.01)-mean(conc_previous))./(std(conc_previous));
%conc_previous=(conc_previous-min(conc_previous))./(max(conc_previous)-min(conc_previous));
%conc_previous=remaining_conc;
%conc_previous =(1./(1+exp(-10.*MLT_shifted_winter))).*conc_previous0';
%conc_previous2 = (1./(1+exp(-MLT_shifted_winter'))).*conc_previous;
%%% Select Fit period: I use all the data to compute the regression: we don't want to forecats, we want to "explain".
%%%%%repeat last two year SST and DUST
remaining_conc=zeros(168,1);
DUST1=[DUST; DUST(end-23:end)];
SST_nonlinear1=[SST_nonlinear; SST_nonlinear(end-23:end)];

clear R1test
clear MSE1test
%for time_window=[1:1:6];%3;
counter=0;


clear yp_all
%for RRINDEX=[24:-time_window:0]%[27:-3:0]
for RRINDEX=[1 1 1 1 1 1 1 1 1].*NORMINDEX2;
counter=counter+1;
notemperature=4;
% I fit over all the data: 
Training_ini = 1;                        % starts from timestep 1
Training_end = length(observation_dates)-RRINDEX % ends at the last timestep
clear tbl

%%%%notemperature = 1 orignal
%%%%notemperature = 2 orignal + SST
%%%%notemperature = 3 Add MLD correction term
%%%%notemperature = 4 Add MLD correction term + SST
%%%%notemperature = 5 Add MLD correction term + SST + exclude MLD
%%%%notemperature = 6 Add MLD correction term + SST + exclude Clag
%%%%notemperature = 7 Add MLD correction term + SST + exclude Dust
%%%%notemperature == 8, all terms is linear
if notemperature == 1
tbl = table(MLT_shifted(Training_ini:Training_end),conc_previous(Training_ini:Training_end),DUST(Training_ini:Training_end),conc_current_year(Training_ini:Training_end)); % Orginal
elseif notemperature == 2
tbl =table(MLT_shifted(Training_ini:Training_end),conc_previous(Training_ini:Training_end),DUST(Training_ini:Training_end),SST_nonlinear(Training_ini:Training_end),conc_current_year(Training_ini:Training_end)); %Add SST
elseif notemperature == 3
tbl =table(MLT_avgnorm0(Training_ini:Training_end),conc_previous(Training_ini:Training_end),DUST(Training_ini:Training_end),SST_nonlinear(Training_ini:Training_end),MLT_avgnorm1(Training_ini:Training_end),remaining_conc(Training_ini:Training_end),conc_current_year(Training_ini:Training_end)); 
elseif notemperature == 4
tbl =table(MLT_shifted(Training_ini:Training_end),conc_previous(Training_ini:Training_end),DUST1(Training_ini:Training_end),SST_nonlinear1(Training_ini:Training_end),MLT_avgnorm1(Training_ini:Training_end),remaining_conc(Training_ini:Training_end),conc_current_year(Training_ini:Training_end)); 
%tbl =table(MLT_avgnorm0(Training_ini:Training_end),conc_previous(Training_ini:Training_end),DUST1(Training_ini:Training_end),SST_nonlinear1(Training_ini:Training_end),MLT_avgnorm1(Training_ini:Training_end),remaining_conc(Training_ini:Training_end),conc_current_year(Training_ini:Training_end)); 
elseif notemperature == 5
tbl =table(MLT_avgnorm0(Training_ini:Training_end),conc_previous(Training_ini:Training_end),DUST(Training_ini:Training_end),SST_nonlinear(Training_ini:Training_end),MLT_avgnorm1(Training_ini:Training_end),remaining_conc(Training_ini:Training_end),conc_current_year(Training_ini:Training_end)); 
elseif notemperature ==6
tbl =table(MLT_avgnorm0(Training_ini:Training_end),conc_previous(Training_ini:Training_end),DUST(Training_ini:Training_end),SST_nonlinear(Training_ini:Training_end),MLT_avgnorm1(Training_ini:Training_end),remaining_conc(Training_ini:Training_end),conc_current_year(Training_ini:Training_end)); 
elseif notemperature == 7
tbl =table(MLT_avgnorm0(Training_ini:Training_end),conc_previous(Training_ini:Training_end),DUST(Training_ini:Training_end),SST_nonlinear(Training_ini:Training_end),MLT_avgnorm1(Training_ini:Training_end),remaining_conc(Training_ini:Training_end),conc_current_year(Training_ini:Training_end)); 
elseif notemperature == 8
tbl =table(MLT_avgnorm0(Training_ini:Training_end),conc_previous(Training_ini:Training_end),DUST(Training_ini:Training_end),SST_nonlinear(Training_ini:Training_end),MLT_avgnorm1(Training_ini:Training_end),remaining_conc(Training_ini:Training_end),conc_current_year(Training_ini:Training_end)); 
else
end

%---------------------------------------------------------------------------------------------
% x1 = mlt, x2 = prev conc , x3 = dust component 
%-----------------------------------------------
% With a model function del tipo: 
% ----------------------------------------------

% Function of this Type X1:
if notemperature == 1
modelfun = @(b,x)(b(1)*(exp(x(:,1)).^abs(b(2))) + ((x(:,2).^abs(b(3)))) + b(4)*((x(:,3)))).*sign(heaviside((b(1)*exp((x(:,1))).^abs(b(2)) + ((x(:,2).^abs(b(3)))) + b(4)*((x(:,3)))))); %%%ORGINAL
elseif notemperature == 2
modelfun = @(b,x)(b(1)*(exp(x(:,1)).^abs(b(2))) + ((x(:,2).^abs(b(3)))) + b(4)*((x(:,3))) + b(5)*((x(:,4)))).*sign(heaviside((b(1)*exp((x(:,1))).^abs(b(2)) + ((x(:,2).^abs(b(3)))) + b(4)*((x(:,3))) + b(5)*((x(:,4)))))); %%%Add SST
%modelfun = @(b,x)(((x(:,2).^abs(b(1)))) + b(2)*((x(:,3))) + b(3)*((x(:,4)))).*sign(heaviside(((x(:,2).^abs(b(1)))) + b(2)*((x(:,3))) + b(3)*((x(:,4))))); %%%Add SST
elseif notemperature == 3
modelfun = @(b,x)(b(1).*(x(:,1)).^abs(b(2)) + b(3)*exp(1-x(:,5)).*x(:,2) + b(4).*((x(:,3)))).*sign(heaviside(b(1).*(x(:,1)).^abs(b(2)) + b(3)*exp(1-x(:,5)).*x(:,2) + b(4).*((x(:,3))))); %%%Add SST
elseif notemperature == 4
%modelfun = @(b,x)(b(1)*(exp(x(:,1)).^abs(b(2))) + ((x(:,2).^abs(b(3)))) + b(4)*((x(:,3))) + b(5)*((x(:,4)))+b(6)*(x(:,2).*exp(b(7).*(x(:,5)-b(8))))).*sign(heaviside((b(1)*exp((x(:,1))).^abs(b(2)) + ((x(:,2).^abs(b(3)))) + b(4)*((x(:,3))) +b(5)*((x(:,4)))+b(6)*(x(:,2).*exp(b(7).*(x(:,5)-b(8))))))); %%%Add SST
%modelfun = @(b,x)(b(1).*(x(:,1)).^abs(b(5)) + b(2)*exp(1-x(:,5)).*x(:,2) + b(3).*((x(:,3))) + b(4)*((x(:,4)))).*sign(heaviside(b(1).*(x(:,1)).^abs(b(5)) + b(2)*exp(1-x(:,5)).*x(:,2) + b(3).*((x(:,3))) +b(4)*((x(:,4))))); %%%Add SST
%modelfun = @(b,x)(b(1).*(x(:,1)).^abs(3) + b(2)*(1-x(:,5).^3).*x(:,2) + b(3).*((x(:,3))) + b(4)*((x(:,4)))).*sign(heaviside(b(1).*(x(:,1)).^abs(3) +b(2)*(1-x(:,5).^3).*x(:,2) + b(3).*((x(:,3))) +b(4)*((x(:,4))))); %%%Add SST
%modelfun = @(b,x)(-abs(b(6))+b(1).*(x(:,1)).^abs(b(5)) + b(2)*(1-x(:,5).^3).*x(:,2) + b(3).*((x(:,3))) + b(4)*((x(:,4)))).*sign(heaviside(-abs(b(6))+b(1).*(x(:,1)).^abs(b(5)) +b(2)*exp(1-x(:,5)).*x(:,2) + b(3).*((x(:,3))) +b(4)*((x(:,4)))));

modelfun = @(b,x) ( ...
    (x(:,1) >= 0) .* ( ...
        (b(1) .* (x(:,1)).^abs(b(5)) + b(2) .* (1 - x(:,5).^3) .* x(:,2) + b(3) .* x(:,3) + b(4) .* x(:,4)) .* ...
        sign(heaviside(b(1) .* (x(:,1)).^abs(b(5)) + b(2) .* (1 - x(:,5).^3) .* x(:,2) + b(3) .* x(:,3) + b(4) .* x(:,4))) ...
    ) + ...
    (x(:,1) < 0) .* ( ...
        (-b(1) .* abs(x(:,1)).^abs(b(5)) + b(2) .* (1 - x(:,5).^3) .* x(:,2) + b(3) .* x(:,3) + b(4) .* x(:,4)) .* ...
        sign(heaviside(-b(1) .* abs(x(:,1)).^abs(b(5)) + b(2) .* (1 - x(:,5).^3) .* x(:,2) + b(3) .* x(:,3) + b(4) .* x(:,4))) ...
    ) ...
);
elseif notemperature ==5 %%%exclude MLD
modelfun = @(b,x)(b(1).*x(:,2) + b(2).*((x(:,3))) + b(3)*((x(:,4)))).*sign(heaviside(b(1).*x(:,2) + b(2).*((x(:,3))) +b(3)*((x(:,4))))); %%%Add SST
%modelfun = @(b,x)(b(1)*(exp(x(:,1)).^abs(b(2))) + ((x(:,2).^abs(b(3)))) + b(4)*((x(:,3))) + b(5)*((x(:,4)))+b(6)*(x(:,6).*exp(b(7).*(x(:,5)-b(8))))).*sign(heaviside((b(1)*exp((x(:,1))).^abs(b(2)) + ((x(:,2).^abs(b(3)))) + b(4)*((x(:,3))) +b(5)*((x(:,4)))+b(6)*(x(:,6).*exp(b(7).*(x(:,5)-b(8))))))); 
elseif notemperature == 6
modelfun = @(b,x)(b(1).*(x(:,1)).^abs(b(5)) + b(2)*exp(1-x(:,5)) + b(3).*((x(:,3))) + b(4)*((x(:,4)))).*sign(heaviside(b(1).*(x(:,1)).^abs(b(5)) + b(2)*exp(1-x(:,5)) + b(3).*((x(:,3))) +b(4)*((x(:,4))))); %%%Add SST
elseif notemperature == 7
modelfun = @(b,x)(b(1).*(x(:,1)).^abs(b(2)) + b(3)*exp(1-x(:,5)).*x(:,2)  + b(4)*((x(:,4)))).*sign(heaviside(b(1).*(x(:,1)).^abs(b(2)) + b(3)*exp(1-x(:,5)).*x(:,2) +b(4)*((x(:,4))))); %%%Add SST
elseif notemperature == 8
modelfun = @(b,x)(b(1).*(x(:,1)) + b(2).*x(:,2) + b(3).*((x(:,3))) + b(4).*((x(:,4)))).*sign(heaviside(b(1).*(x(:,1)) + b(2).*x(:,2) + b(3).*((x(:,3))) +b(4).*((x(:,4))))); %%%Add SST
else
end

% Pseudo Random initial values chosen between at i-th iteration: 
% pseudorandom values are drawn from the standard uniform distribution on the open interval(0,1).
% (uniform distributions are essentially probability distributions with equally likely outcomes)
if notemperature == 1
     beta0 = [rand(1,1) rand(1,1) rand(1,1) rand(1,1)];
elseif notemperature == 2
     beta0 = [rand(1,1) rand(1,1) rand(1,1) rand(1,1) rand(1,1)];
      %  beta0 = [rand(1,1) rand(1,1) rand(1,1)];
elseif notemperature == 3
        beta0 = 5.0.*[0.1 0.1 0.1 0.1];
elseif notemperature == 4     
     %beta0 = [rand(1,1) rand(1,1) rand(1,1) rand(1,1) rand(1,1)];
    % beta0 = [rand(1,1) rand(1,1) rand(1,1) rand(1,1) rand(1,1) rand(1,1) rand(1,1) rand(1,1)];
     % beta0 = [0.9 0.9 0.9 0.9 0.9 0.9 0.9];
     beta0 = 5.0.*[0.1 0.1 0.1 0.1 0.1];
    % beta0 = 5.0.*[0.1 0.1 0.1 0.1];
elseif notemperature == 5
   beta0 = 5.0.*[0.1 0.1 0.1];
 % beta0 = [rand(1,1) rand(1,1) rand(1,1) rand(1,1)];
elseif notemperature == 6
   beta0 = 5.0.*[0.1 0.1 0.1 0.1 0.1];
elseif notemperature == 7
   beta0 = 5.0.*[0.1 0.1 0.1 0.1];
elseif notemperature == 8
   beta0 = 5.0.*[0.1 0.1 0.1 0.1];
end
   
   
% Fit non linear model:    
mdl_nl = fitnlm(tbl,modelfun,beta0); 

% Table to array to extract the coefficients value 
Bi = table2array(mdl_nl.Coefficients);
b1 = Bi(1,1);
b2 = Bi(2,1);
b3 = Bi(3,1);
if notemperature ~=5
b4 = Bi(4,1);
end

if notemperature ~= 1
if notemperature == 2
    b5 = Bi(5,1);
    %clear b4
elseif notemperature == 3
   %b4 = Bi(4,1);
   % b6 = Bi(5,1);  
   % b7 = Bi(6,1);
   % b8 = Bi(7,1);
elseif notemperature == 4
    b5 = Bi(5,1);
  %  b6 = Bi(6,1);  
  %  b7 = Bi(7,1);
  %  b8 = Bi(8,1);
elseif notemperature == 5
 %clear b4
elseif notemperature == 6
    b4 = Bi(4,1);
    b5 = Bi(5,1);
elseif notemperature == 7
    b4 = Bi(4,1);
    %b5 = Bi(5,1);
    %b6 = Bi(6,1);  
    %b7 = Bi(7,1);
elseif notemperature == 8
     b4 = Bi(4,1);
else
end
end


%%%
% ==========================================================================
%  Final Estimator: with Model Function of the chosen Type above: 
% ==========================================================================
%b3=0;
%b4=0;

% Since I used the Type X1 above, it is: 
if notemperature == 1
yp1 = (b1*exp(MLT_shifted(1:end)).^abs(b2) + (conc_previous.^abs(b3)) + b4*(DUST(1:end))).*sign(heaviside(b1*exp(MLT_shifted(1:end)).^abs(b2) + (conc_previous.^abs(b3)) + b4*(DUST(1:end))));
elseif notemperature == 2
yp1 = (b1*exp(MLT_shifted(1:end)).^abs(b2) + (conc_previous.^abs(b3)) + b4*(DUST(1:end))+b5*(SST_nonlinear(1:end))).*sign(heaviside(b1*exp(MLT_shifted(1:end)).^abs(b2) + (conc_previous.^abs(b3)) + b4*(DUST(1:end))+b5*(SST_nonlinear(1:end))));
%yp1 = ((conc_previous.^abs(b3)) + b4*(DUST(1:end))+b5*(SST_nonlinear(1:end))).*sign(heaviside((conc_previous.^abs(b3)) + b4*(DUST(1:end))+b5*(SST_nonlinear(1:end))));
elseif notemperature == 3
yp1 = (b1.*(MLT_avgnorm0(1:end)).^abs(b2) + b3*exp(1-MLT_avgnorm1(1:end)).*conc_previous + b4.*(DUST(1:end))).*sign(heaviside(b1.*(MLT_avgnorm0(1:end)).^abs(b2) + b3*exp(1-MLT_avgnorm1(1:end)).*conc_previous + b4.*(DUST(1:end))));
elseif notemperature== 4
%yp1 = (b1*exp(MLT_shifted(1:end)).^abs(b2) + (conc_previous.^abs(b3)) + b4*(DUST(1:end))+b5*(SST_nonlinear(1:end))+b6*(conc_previous(1:end).*exp(b7.*(MLT_avgnorm(1:end)-b8)))).*sign(heaviside(b1*exp(MLT_shifted(1:end)).^abs(b2) + (conc_previous.^abs(b3)) + b4*(DUST(1:end))+b5*(SST_nonlinear(1:end))+b6*(conc_previous(1:end).*exp(b7.*(MLT_avgnorm(1:end)-b8)))));
%yp1 = (b1.*(MLT_avgnorm0(1:end)).^abs(b5) + b2*exp(1-MLT_avgnorm1(1:end)).*conc_previous + b3.*(DUST(1:end))+b4*(SST_nonlinear(1:end))).*sign(heaviside(b1.*(MLT_avgnorm0(1:end)).^abs(b5) + b2*exp(1-MLT_avgnorm1(1:end)).*conc_previous + b3.*(DUST(1:end))+b4*(SST_nonlinear(1:end))));
%ypold1 = (b1.*(MLT_shifted(1:end)).^abs(3) + b2*(1-MLT_avgnorm1(1:end).^3).*conc_previous(1:end) + b3.*(DUST1(1:end))+b4*(SST_nonlinear1(1:end))).*sign(heaviside(b1.*(MLT_shifted(1:end)).^abs(3) + b2*(1-MLT_avgnorm1(1:end).^3).*conc_previous(1:end) + b3.*(DUST1(1:end))+b4*(SST_nonlinear1(1:end))));
%yp1 = (b1.*(MLT_avgnorm0(1:end)).^abs(b5) + b2*(1-MLT_avgnorm1(1:end).^3).*conc_previous)*sign(heaviside(b1.*(MLT_avgnorm0(1:end)).^abs(b5) + b2*exp(1-MLT_avgnorm1(1:end)).*conc_previous));
%yp1 = (b1.*(MLT_shifted(1:end)).^abs(3) + b2*(1- MLT_avgnorm_test1(1:end).^3).*conc_previous).*sign(heaviside(b1.*(MLT_shifted(1:end)).^abs(3) + b2*(1- MLT_avgnorm_test1(1:end).^3).*conc_previous));

for yptest1=1:length(MLT_shifted)
if MLT_shifted(yptest1) >=0
ypold1(yptest1) = (b1.*(MLT_shifted(yptest1)).^abs(b5) + b2*(1-MLT_avgnorm1(yptest1).^3).*conc_previous(yptest1) + b3.*(DUST1(yptest1))+b4*(SST_nonlinear1(yptest1))).*sign(heaviside(b1.*(MLT_shifted(yptest1)).^abs(b5) + b2*(1-MLT_avgnorm1(yptest1).^3).*conc_previous(yptest1) + b3.*(DUST1(yptest1))+b4*(SST_nonlinear1(yptest1))));
yp1(yptest1) = (b1.*(MLT_shifted(yptest1)).^abs(b5) + b2*(1- MLT_avgnorm_test1(yptest1).^3).*conc_previous(yptest1)).*sign(heaviside(b1.*(MLT_shifted(yptest1)).^abs(b5) + b2*(1- MLT_avgnorm_test1(yptest1).^3).*conc_previous(yptest1)));
else
ypold1(yptest1) = (-b1.*(abs(MLT_shifted(yptest1))).^abs(b5) + b2*(1-MLT_avgnorm1(yptest1).^3).*conc_previous(yptest1) + b3.*(DUST1(yptest1))+b4*(SST_nonlinear1(yptest1))).*sign(heaviside(-b1.*(abs(MLT_shifted(yptest1))).^abs(b5) + b2*(1-MLT_avgnorm1(yptest1).^3).*conc_previous(yptest1) + b3.*(DUST1(yptest1))+b4*(SST_nonlinear1(yptest1))));
yp1(yptest1) = (-b1.*(abs(MLT_shifted(yptest1))).^abs(b5) + b2*(1- MLT_avgnorm_test1(yptest1).^3).*conc_previous(yptest1)).*sign(heaviside(-b1.*(abs(MLT_shifted(yptest1))).^abs(b5) + b2*(1- MLT_avgnorm_test1(yptest1).^3).*conc_previous(yptest1)));
end 
end

elseif notemperature == 5
yp1 = (b1.*conc_previous + b2.*(DUST(1:end))+b3*(SST_nonlinear(1:end))).*sign(heaviside(b1.*conc_previous + b2.*(DUST(1:end))+b3*(SST_nonlinear(1:end))));
elseif notemperature == 6
yp1 = (b1.*(MLT_avgnorm0(1:end)).^abs(b5) + b2*exp(1-MLT_avgnorm1(1:end)) + b3.*(DUST(1:end))+b4*(SST_nonlinear(1:end))).*sign(heaviside(b1.*(MLT_avgnorm0(1:end)).^abs(b5) + b2*exp(1-MLT_avgnorm1(1:end)) + b3.*(DUST(1:end))+b4*(SST_nonlinear(1:end))));
elseif notemperature == 7
yp1 = (b1.*(MLT_avgnorm0(1:end)).^abs(b2) + b3*exp(1-MLT_avgnorm1(1:end)).*conc_previous + b4*(SST_nonlinear(1:end))).*sign(heaviside(b1.*(MLT_avgnorm0(1:end)).^abs(b2) + b3*exp(1-MLT_avgnorm1(1:end)).*conc_previous + b4*( SST_nonlinear(1:end))));
elseif notemperature == 8
yp1 = (b1.*(MLT_avgnorm0(1:end)) + b2.*conc_previous + b3.*(DUST(1:end))+b4*(SST_nonlinear(1:end))).*sign(heaviside(b1.*(MLT_avgnorm0(1:end)) + b2.*conc_previous + b3.*(DUST(1:end))+b4*(SST_nonlinear(1:end))));
else
end

yp_all(counter,:)=yp1;

end
%%% Plot the resulting estimation 
%Rall(toty)=corr(obs_webdig,yp1);
%end
%%

yp2=zeros(168,1);
%yp2(1:144)=yp_all(1,1:144);
yp2=yp_all(1,1:168)';

%for i=1:size(yp_all,1)-1
%yp2(144+time_window*(i-1)+1:144+time_window*i)=yp_all(i,144+time_window*(i-1)+1:144+time_window*i)'
%end


%%%%here is b nochange
%%%%control xing
%yp1 = ( b2.*conc_previous +b3.*(DUST(1:end))+b4*(SST_nonlinear(1:end))).*sign(heaviside(b2.*conc_previous+ b3.*(DUST(1:end))+b4*(SST_nonlinear(1:end)))); %exclude MLD
%%%exp()
%yp1 = (b1.*(MLT_avgnorm0(1:end)).^abs(b5) + b2*exp(1-MLT_avgnorm1(1:end)) + b3.*(DUST(1:end))+b4*(SST_nonlinear(1:end))).*sign(heaviside(b1.*(MLT_avgnorm0(1:end)).^abs(b5) + b2*exp(1-MLT_avgnorm1(1:end)) + b3.*(DUST(1:end))+b4*(SST_nonlinear(1:end)))); %exclude Clag
%yp1 = (b1.*(MLT_avgnorm0(1:end)).^abs(b5) + b2*exp(1-MLT_avgnorm1(1:end)).*conc_previous +b4*(SST_nonlinear(1:end))).*sign(heaviside(b1.*(MLT_avgnorm0(1:end)).^abs(b5) + b2*exp(1-MLT_avgnorm1(1:end)).*conc_previous +b4*(SST_nonlinear(1:end))));
%yp1 = (b1.*(MLT_avgnorm0(1:end)).^abs(b5) + b2*exp(1-MLT_avgnorm1(1:end)).*conc_previous + b3.*(DUST(1:end))).*sign(heaviside(b1.*(MLT_avgnorm0(1:end)).^abs(b5) + b2*exp(1-MLT_avgnorm1(1:end)).*conc_previous + b3.*(DUST(1:end))));
%%%MLD^3
%yp1 = (b1.*(MLT_avgnorm0(1:end)).^abs(b5) + b2*(1-MLT_avgnorm1(1:end).^3) + b3.*(DUST(1:end))+b4*(SST_nonlinear(1:end))).*sign(heaviside(b1.*(MLT_avgnorm0(1:end)).^abs(b5) + b2*exp(1-MLT_avgnorm1(1:end)) + b3.*(DUST(1:end))+b4*(SST_nonlinear(1:end)))); %exclude Clag
%yp1 = (b1.*(MLT_avgnorm0(1:end)).^abs(b5) + b2*(1-MLT_avgnorm1(1:end).^3).*conc_previous +b4*(SST_nonlinear(1:end))).*sign(heaviside(b1.*(MLT_avgnorm0(1:end)).^abs(b5) + b2*exp(1-MLT_avgnorm1(1:end)).*conc_previous +b4*(SST_nonlinear(1:end))));
%yp2 = (b1.*(MLT_avgnorm0(1:end)).^abs(b5) + b2*(1-MLT_avgnorm1(1:end).^3).*conc_previous + b3.*(DUST(1:end))).*sign(heaviside(b1.*(MLT_avgnorm0(1:end)).^abs(b5) + b2*exp(1-MLT_avgnorm1(1:end)).*conc_previous + b3.*(DUST(1:end))));
%yp2 = (b1.*(MLT_avgnorm0(1:end)).^abs(b5) + b2*(1-MLT_avgnorm1(1:end).^3).*conc_previous).*sign(heaviside(b1.*(MLT_avgnorm0(1:end)).^abs(b5) + b2*exp(1-MLT_avgnorm1(1:end)).*conc_previous));
%%%% control Lyuba
%yp1 = ((conc_previous.^abs(b3)) + b4*(DUST(1:end))).*sign(heaviside((conc_previous.^abs(b3)) + b4*(DUST(1:end)))); %exclude MLD
%yp1 = (b1*exp(MLT_shifted(1:end)).^abs(b2) + b4*(DUST(1:end))).*sign(heaviside(b1*exp(MLT_shifted(1:end)).^abs(b2)+ b4*(DUST(1:end)))); %exclude Clag
%yp1 = (b1*exp(MLT_shifted(1:end)).^abs(b2) + (conc_previous.^abs(b3)) ).*sign(heaviside(b1*exp(MLT_shifted(1:end)).^abs(b2) + (conc_previous.^abs(b3)))); %exclude DUST

yp3=[ypold1(1:NORMINDEX);  yp2(NORMINDEX+1:end)];
figure(23)
clf
scale2=80;
scale=1; 
set(gcf,'position', [474   157     12*scale2 6*scale2])
set(gcf,'paperposition', [0   0   12*scale  6*scale])

set(gcf,'units','centimeters','paperunits','centimeters');
%setfigsize(gcf,25,10) % 30 15 
box on 
%xlabel('','Interpreter','latex');
%ylabel('Million tons','Interpreter','latex');
%title('Monthly mean Sargassum biomass','Interpreter','latex');
grid on 
box on 
set(gca,'GridLineStyle',':')

%set(gca,'FontSize',15,'FontName','Times');


% Observations:
%if false
%h3 = bar(datevector(a:b),obs_webdig)
%h3 = bar(datevector(a:b),conc_current_year)
%h3.FaceColor = [0. 0.55 0.9] % azzurrognolo
%h3.FaceColor = [0.9 0.75 0.]  % giallognolo
%h3.FaceColor = [65,174,118]./256;  % giallognolo
%h3.FaceColor = [.3 .3 .3];
%h3.FaceAlpha = 0.7;%1.0;%0.7;
%h3.FaceAlpha = 1.0;%1.0;%0.7;
%h3.EdgeColor = 'none';

%else
%end

hold on
a2=area(datevector(a:b-NORMINDEX2),yp3(1:NORMINDEX));
a2.FaceColor=[65,174,118]./256;
a2.FaceAlpha = 0.7;
a2.EdgeColor = 'none';
h2=plot(datevector(a:b-NORMINDEX2),yp3(1:NORMINDEX),'-o');hold on
h2.Color=[65,174,118]./256;
h2.LineWidth = 1.5;
h2.MarkerSize = 4;
h2.MarkerFaceColor = [65,174,118]./256;
%h2 = bar(datevector(a:b-NORMINDEX2),yp3(1:NORMINDEX)) % non serve più normalizzare.
%h2.FaceColor = [.5 .5 .5] % if the other is azzurrognolo
%h2.FaceColor = [.3 .3 .3];
%h2.FaceColor = [65,174,118]./256; 
%h2.FaceColor = [5,112,176]./256;
%h2.BarWidth = 0.6
%h2.FaceAlpha = 1.0;
%h2.FaceAlpha = 0.7;
%h2.EdgeColor = 'none';
%hold on 
%title('')

hold on
a4=area(datevector(end-NORMINDEX2:end),yp3(end-NORMINDEX2:end));
a4.FaceColor=[34,94,168]./256;
a4.FaceAlpha = 0.7;
a4.EdgeColor = 'none';
h4=plot(datevector(end-NORMINDEX2:end),yp3(end-NORMINDEX2:end),'-o');hold on
h4.Color=[34,94,168]./256;
h4.LineWidth = 1.5;
h4.MarkerSize = 4;
h4.MarkerFaceColor = [34,94,168]./256;
%h4 = bar(datevector(end-NORMINDEX2:end),yp2(end-NORMINDEX2:end)) 

%h2.FaceColor = [.5 .5 .5] % if the other is azzurrognolo
%h2.FaceColor = [.3 .3 .3];
%h4.FaceColor = [34,94,168]./256; 
%h2.FaceColor = [5,112,176]./256;
%h4.BarWidth = 0.6
%h2.FaceAlpha = 1.0;
%h4.FaceAlpha = 0.7;
%h4.EdgeColor = 'none';
hold on 
h3 = plot(datevector(a:b), obs_webdig, '-o');hold on  % Use '-o' for line with circle markers
% Set line color (gray)
h3.Color = [.3 .3 .3];

% Optional: adjust line width and marker style
h3.LineWidth = 1.5;
h3.MarkerSize = 4;
h3.MarkerFaceColor = [.3 .3 .3];
%title('')


legend([h3 h2 h4],'obs (digitized)','model (trained)','model (predicted)','Location','northwest')

ylim([0 25])
set(gca,'FontSize',15); 
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
set(gca,'TickDir','out')
xlabel('Years','Fontsize',15)
ylabel('Sargassum concentrations [Mton]')
%%%ALL R VALUE


R1=corr(obs_webdig(1:NORMINDEX), yp3(1:NORMINDEX))
MSE1= mean((yp3(1:NORMINDEX) - obs_webdig(1:NORMINDEX)).^2)
%text(2025,23,['R=' sprintf('%.2f', corr(obs_webdig(1:NORMINDEX), ypold1(1:NORMINDEX)))],'Fontsize',15,'color',[65,174,118]./256)
%text(2025, 21, ['MSE = ' sprintf('%.2f', mean((ypold1(1:NORMINDEX) - obs_webdig(1:NORMINDEX)).^2))], 'FontSize', 15,'color',[65,174,118]./256)
text(2025,23,['R=' sprintf('%.2f', corr(obs_webdig(1:end), yp3(1:end)))],'Fontsize',15,'color','k')
text(2025, 21, ['MSE = ' sprintf('%.2f', mean((yp3(1:end) - obs_webdig(1:end)).^2))], 'FontSize', 15,'color','k')


text(2025,19,['R=' sprintf('%.2f', corr(obs_webdig(NORMINDEX+1:end), yp3(NORMINDEX+1:end)))],'Fontsize',15,'color',[34,94,168]./256)
text(2025, 17, ['MSE = ' sprintf('%.2f', mean((yp3(NORMINDEX+1:end) - obs_webdig(NORMINDEX+1:end)).^2))], 'FontSize', 15,'color',[34,94,168]./256)
%%%2024 only
%text(2025,15,['R=' sprintf('%.2f', corr(obs_webdig(end-12:end), yp2(end-12:end)))],'Fontsize',15,'color','k')
%text(2025, 13, ['MSE = ' sprintf('%.2f', mean((yp2(end-12:end) - obs_webdig(end-12:end)).^2))], 'FontSize', 15,'color','k')
%title('(c)')


%R1test(time_window,1)=corr(obs_webdig(end-24:end), yp2(end-24:end));
%R1test(time_window,2)=corr(obs_webdig(end-24:end-13), yp2(end-24:end-13));
%R1test(time_window,3)=corr(obs_webdig(end-12:end), yp2(end-12:end));


%MSE1test(time_window,1)=mean((yp2(end-24:end) - obs_webdig(end-24:end)).^2);
%MSE1test(time_window,2)=mean((yp2(end-24:end-13) - obs_webdig(end-24:end-13)).^2);
%MSE1test(time_window,3)=mean((yp2(end-12:end) - obs_webdig(end-12:end)).^2);

%end

% ===================================================
% PUT THE RIGHT NAME FOR YOUR SAVED FIGURE HERE: 

%print(gcf,'../Figures/Fig_regression_model_reducedGASB_giallino.png','-dpng','-r350')

 % ===================================================
%outfile=['/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/xingplot/statistical_anlyz_linear.png'];print([outfile],'-dpng','-r300');
%outfile=['/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/xingplot/factortest/exclude_dust.png'];print([outfile],'-dpng','-r300');
%outfile=['/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/xingplot/statistical_anlyz_control.png'];print([outfile],'-dpng','-r300');
%outfile=['/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/paperfigure/figure_nonlinearmodel/Fig3_MLDrag3.png'];print([outfile],'-dpng','-r300');

%outfile=['/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/paperfigure/figure_prediction/FigS5_2022to2024_newequation.png'];print([outfile],'-dpng','-r300');
outfile=['/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/utils/newdata/domain/Lat1to15/Regressiontest/FigS5_4_20232024_version2.png'];print([outfile],'-dpng','-r300');
%% 
figure(43)
clf
set_portrait
subplot(2,1,1)
plot(datevector(a:b), MLT, 'k.-', 'LineWidth', 2.0); hold on;
ylabel('MLD anomaly (m)')
set(gca,'fontsize',12)
subplot(2,1,2)
plot(datevector(a:b), b1.*(MLT_avgnorm_test0(1:end)).^abs(b5), 'k.-', 'LineWidth', 2.0); hold on;
plot(datevector(a:b), b2*(1-MLT_avgnorm_test1(1:end).^3).*conc_previous, 'b.-', 'LineWidth', 2.0); hold on;
plot(datevector(a:b), -abs(b6)+b1.*(MLT_avgnorm_test0(1:end)).^abs(b5), 'r.-', 'LineWidth', 2.0); hold on;
plot(datevector(a:b),  b3.*(DUST1(1:end)), 'm.-', 'LineWidth', 2.0); hold on;
plot(datevector(a:b), b4*(SST_nonlinear1(1:end)), 'y.-', 'LineWidth', 2.0); hold on;
ylabel('Terms contribution')
legend('MLD term','MLD and Conc Previous term','location','northwest')
set(gca,'fontsize',12)
%outfile=['/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/paperfigure/figure_prediction/FigS5_explain.png'];print([outfile],'-dpng','-r300');
outfile=['/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/utils/newdata/domain/Lat1to15/Regressiontest/FigS5_1_20232024.png'];print([outfile],'-dpng','-r300');
%outfile=['/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/utils/newdata/domain/Lat1to15/Regressiontest/FigS5_1_20212024.png'];print([outfile],'-dpng','-r300');

%% NOT SURE IF USEFUL: BUT WE COULD VERIFY THE STABILITY OF THE DETERMINED COEFFICIENTS: i.e. by repeating the fitting N times, all the bi should be more or less the same 
clc

% Select how many iterations: N
%close all
N = 50; 

clear b1_i b2_i b3_i b4_i 
clear initial_values 
% Initialize dummy vectors for coefficients 
b1_i(1:N) = -999;
b2_i(1:N) = -999;
b3_i(1:N) = -999;
b4_i(1:N) = -999;
initial_values(1:N,1:4) = NaN;

% Create a loop where you evaluate all the bi coefficietns from scratch at
% each cycle, and you store each of them in a vector. 

% Essentially, you just re-run the non linear regression to obtain the coefficients: 
% (no need to re-define the model function nor the table, because it's the one defined above)

for i = 1:N 
   
   
    % Pseudo Random initial values chosen between at i-th iteration: 
    % pseudorandom values drawn from the standard uniform distribution on the open interval(0,1).
    % (uniform distributions are essentially probability distributions with equally likely outcomes)
     
     beta0_1 = rand(1,1);    
     beta0_2 = rand(1,1); 
     beta0_3 = rand(1,1);  
     beta0_4 = rand(1,1);
      
     beta0_i = [beta0_1 beta0_2 beta0_3 beta0_4];
     
    % Store initial vaklues here just for a lter check 
    initial_values(i,1:4) = beta0_i;
    
    % Fit non linear model, at i-th iteration: (model function and table are the same already defined above.)   
    mdl_nl_i = fitnlm(tbl,modelfun,beta0_i); 

    % Table to array to extract the coefficients value, at i-th iteration:  
    Bi_i = table2array(mdl_nl_i.Coefficients);
    
    b1_i(i) = Bi_i(1,1);
    b2_i(i) = Bi_i(2,1);
    b3_i(i) = Bi_i(3,1);
    b4_i(i) = Bi_i(4,1);
    
    clear beta0_i   
    clear mdl_nl_i
    clear Bi_i
    i
end


% Finally plot the resulting coefficients againts the iteration 
figure 
subplot(2,2,1)
plot(b1_i,'k.')
ylabel 'b1'
xlabel 'iter'

subplot(2,2,2)
plot(abs(b2_i),'k.')
ylabel '|b2|'
xlabel 'iter'

subplot(2,2,3)
plot(abs(b3_i),'k.')
ylabel '|b3|'
xlabel 'iter'

subplot(2,2,4)
plot(b4_i,'k.')
ylabel 'b4'
xlabel 'iter'