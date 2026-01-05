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
resultspath1 = '../utils/newdata/domain/Lat1to15/fldm_of_1993_2022_glorys_mlt_reducedGASB_masked_deseas_1-12deg_Lat_1-15.nc' 
                                   
% Path to dust timeseries:

resultspath3 = '../utils/newdata/domain/Lat1to15/fldm_of_duaod550_f64_res075_2003-2022_deseas_reducedGASB_Lat_1-15.nc' 

SSTmean=ncread('/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/utils/newdata/domain/Lat1to15/fldm_of_1993_2022_glorys_sst_reducedGASB_masked_deseas_1-12deg_Lat_1-15.nc','thetao');
SSTmean1=squeeze(SSTmean);

%%
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

TAB = csvread('../utils/digitized_sarg_biomass/WebDigitized_Jan2011-Dec2022_from_Jouanno-et-al-2023.csv');

obs_webdig = TAB(:,2); 

% It's impossible to have negative concentration. Manual digitizing might
% introduce spurious negative values when concentrations are very small. 
% Hence, put to zero every negative point.

obs_webdig(obs_webdig<0) = 0;


% Create date vector for digitized observations:  Jan 2011 to Dec 2022
observation_dates = datetime(2011,1:length(obs_webdig),15)



%% TEMPORAL SUBSET OF MLT AND DUST TIMESERIES: 
% Select, in the glorys data vector, values on timepoints that overlap with observations: 
% i.e. subset the glorys dataset from Jan 2011 to Dec 2020.

% ====================
% 1. Focus on mlt 
% ====================

% =============================
% 1.1. Subset mlt timeseries: 
% =============================

a = 18*12+1             % index where I have Jan 2011 in the glorys date vector. Please always crosscheck.
b = length(datevector); % because Dec 2020 is already the last timestep for the glorys dataset loaded here. Pls always crosscheck. 

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
if datevector(b) == datetime(2022,12,15)
   disp '' 
else 
   error 'datevector(a) must be Dec 2022 but it is not '
end 

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


MLT = (movmean(squeeze(mlt_flmd_deseas(a:b)),Window_for_MLT_timeser)); % Movmean mlt signal, less noise.
% --------------------------------------------------------------------------------------------------------------------------------------
% 1.2.2 Compute here the lagged correlation (xcorr) between :(Sarg Biomass - mean(SargBiomass)) and (smoothed MLT - mean(smoothed MLT)): 
% --------------------------------------------------------------------------------------------------------------------------------------

% Attenzione! Se usi xcorr devi togliere le medie, se i segnali non sono a media nulla, 
% per avere a zero-lag lo stesso risultato di corcoeff o corr (cioè il Pearson correlation coeff). 

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


%%
% 3. Focus on SST 


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


%% NON LINEAR REGRESSION MODEL: 

% -------------------------------------------
%  Try to Fit Non-linear regression model:
% -------------------------------------------

% Define table for variables and observations:
 
% Concentration of Sargassum in current month: i.e. at each month (timestep) the existing bloom (sargassum biomass) 

conc_current_year = obs_webdig; % THIS IS JUST BECAUSE I WANT TO CHANGE THE variable NAME... 
%conc_current_year = (conc_current_year-min(conc_current_year))./(max(conc_current_year)-min(conc_current_year));

% ====================================================================== 
%%% Previous concentration of Sargassum, self-seeding effect:
% In each month, I use the average seasonal biomass 
% of the corresponding season in the previous year
% ======================================================================

% First, compute the seasonal averages:
clear x
x(1:b_dust-a_dust+1) = NaN;

weight=[ones(6,1);zeros(6,1)];%1-exp(-0.1.*[12:-1:0])';
lagmonth=12;

%MLT_avgnorm2=MLT_avgnorm1;
%MLT_avgnorm1=zeros(144,1);

%%%%Sargassyn biomass movmean for 3 months
x=movmean(obs_webdig,3);

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

%%
clear Rall
Rall(1:30)=NaN;
for toty=[12] %%% redo 12 months lag here, no diffirence from previous 
%toty=12;

CONC_LAG = toty;
conc_previous = circshift(x,CONC_LAG);   
if CONC_LAG>0
    conc_previous(1:CONC_LAG) = 0; 
else
    conc_previous(end+CONC_LAG+1:end) = 0;
end


% First, compute the seasonal averages:
clear MLT_seasonal1
MLT_seasonal1(1:b_dust-a_dust+1) = NaN;
for k = 1:(length(MLT_seasonal1)/3)-1
    MLT_seasonal1(3*k:(3*k)+2) = mean(MLT(3*k:(3*k)+2));
end
MLT_seasonal1(1:2) = mean(MLT(1:2));
MLT_seasonal1(end) = mean(MLT(end));
MLT_seasonal1 = MLT_seasonal1';
MLT_avg = circshift(MLT_seasonal1,toty); 


for t = toty+1:144
    % Compute the 12-month average MLT anomaly and 12-month average conc
    %MLT_avg(t) = mean(MLT(t-toty:t));   
    remaining_conc(t) = mean(conc_current_year(t-toty:t));
end

%%%doing normalizarion
 MLT_avgnorm1=(MLT_avg-min(MLT_avg))./(max(MLT_avg)-min(MLT_avg));
 MLT_avgnorm0=(MLT_shifted-min(MLT_shifted))./(max(MLT_shifted)-min(MLT_shifted));
 DUST_avgnorm=(DUST-min(DUST))./(max(DUST)-min(DUST));
 SST_avgnorm=(SST_nonlinear-min(SST_nonlinear))./(max(SST_nonlinear)-min(SST_nonlinear));
 if size(remaining_conc,2)==144
 remaining_conc=remaining_conc';
 end
alpha = 2.0;



conc_previous = [conc_previous zeros(1,length(MLT_shifted)-length(obs_webdig))]; % Se la plottassi vedrei esattamente cosa voglio. 

%%% Select Fit period: I use all the data to compute the regression: we don't want to forecats, we want to "explain".




% I fit over all the data: 
Training_ini = 1;                        % starts from timestep 1
Training_end = length(observation_dates); % ends at the last timestep
clear tbl

%%%%notemperature = 1 orignal
%%%%notemperature = 2 orignal + SST
%%%%notemperature = 3 Add MLD correction term
%%%%notemperature = 4 Add MLD correction term + SST
%%%%notemperature = 5 Add MLD correction term + SST + exclude MLD
%%%%notemperature = 6 Add MLD correction term + SST + exclude Clag
%%%%notemperature = 7 Add MLD correction term + SST + exclude Dust
%%%%notemperature == 8, all terms is linear

tbl =table(MLT_shifted(Training_ini:Training_end),conc_previous(Training_ini:Training_end),DUST(Training_ini:Training_end),SST_nonlinear(Training_ini:Training_end),MLT_avgnorm1(Training_ini:Training_end),remaining_conc(Training_ini:Training_end),conc_current_year(Training_ini:Training_end)); 
%tbl =table(MLT_avgnorm0(Training_ini:Training_end),conc_previous(Training_ini:Training_end),DUST(Training_ini:Training_end),SST_nonlinear(Training_ini:Training_end),MLT_avgnorm1(Training_ini:Training_end),remaining_conc(Training_ini:Training_end),conc_current_year(Training_ini:Training_end)); 


%---------------------------------------------------------------------------------------------
% x1 = mlt, x2 = prev conc , x3 = dust component 
%-----------------------------------------------
% With a model function del tipo: 
% ----------------------------------------------

% Function of this Type X1:

%modelfun = @(b,x)(b(1)*(exp(x(:,1)).^abs(b(2))) + ((x(:,2).^abs(b(3)))) + b(4)*((x(:,3))) + b(5)*((x(:,4)))+b(6)*(x(:,2).*exp(b(7).*(x(:,5)-b(8))))).*sign(heaviside((b(1)*exp((x(:,1))).^abs(b(2)) + ((x(:,2).^abs(b(3)))) + b(4)*((x(:,3))) +b(5)*((x(:,4)))+b(6)*(x(:,2).*exp(b(7).*(x(:,5)-b(8))))))); %%%Add SST
%modelfun = @(b,x)(b(1).*(x(:,1)).^abs(b(5)) + b(2)*exp(1-x(:,5)).*x(:,2) + b(3).*((x(:,3))) + b(4)*((x(:,4)))).*sign(heaviside(b(1).*(x(:,1)).^abs(b(5)) + b(2)*exp(1-x(:,5)).*x(:,2) + b(3).*((x(:,3))) +b(4)*((x(:,4))))); %%%Add SST
%modelfun = @(b,x)(b(1).*(x(:,1)).^abs(b(5)) + b(2)*(1-x(:,5).^3).*x(:,2) + b(3).*((x(:,3))) + b(4)*((x(:,4)))).*sign(heaviside(b(1).*(x(:,1)).^abs(b(5)) + b(2)*exp(1-x(:,5)).*x(:,2) + b(3).*((x(:,3))) +b(4)*((x(:,4))))); %%%Add SST
%modelfun = @(b,x)(b(1).*(x(:,1)).^abs(b(5)) + b(2)*(1-x(:,5).^3).*x(:,2) + b(3).*((x(:,3))) + b(4)*((x(:,4)))).*sign(heaviside(b(1).*(x(:,1)).^abs(b(5)) + b(2)*(1-x(:,5).^3).*x(:,2) + b(3).*((x(:,3))) +b(4)*((x(:,4))))); %%%Add SST
%modelfun = @(b,x)(b(1).*(x(:,1)).^abs(b(5)) + b(2)*(1-x(:,5).^3).*x(:,2) + b(3).*((x(:,3))) + b(4)*((x(:,4)))).*sign(heaviside(b(1).*(x(:,1)).^abs(b(5)) + b(2)*(1-x(:,5).^3).*x(:,2) + b(3).*((x(:,3))) +b(4)*((x(:,4))))); %%%Add SST
modelfun = @(b,x) ( ...
    (x(:,1) >= 0) .* ( ...
        (b(1) .* (x(:,1)).^abs(b(5)) + b(2) .* (1-x(:,5).^3) .* x(:,2) + b(3) .* x(:,3) + b(4) .* x(:,4)) .* ...
        sign(heaviside(b(1) .* (x(:,1)).^abs(b(5)) + b(2) .* (1-x(:,5).^3) .* x(:,2) + b(3) .* x(:,3) + b(4) .* x(:,4))) ...
    ) + ...
    (x(:,1) < 0) .* ( ...
        (-b(1) .* abs(x(:,1)).^abs(b(5)) + b(2) .* (1-x(:,5).^3) .* x(:,2) + b(3) .* x(:,3) + b(4) .* x(:,4)) .* ...
        sign(heaviside(-b(1) .* abs(x(:,1)).^abs(b(5)) + b(2) .* (1-x(:,5).^3) .* x(:,2) + b(3) .* x(:,3) + b(4) .* x(:,4))) ...
    ) ...
);



% Pseudo Random initial values chosen between at i-th iteration: 
% pseudorandom values are drawn from the standard uniform distribution on the open interval(0,1).
% (uniform distributions are essentially probability distributions with equally likely outcomes)

     % beta0 = [0.9 0.9 0.9 0.9 0.9 0.9 0.9];
   %  beta0 = 5.0.*[0.1 0.1 0.1 0.1 0.1 0.1];
     beta0 = 5.0.*[0.1 0.1 0.1 0.1 0.1];

   
   
% Fit non linear model:    
mdl_nl = fitnlm(tbl,modelfun,beta0); 

% Table to array to extract the coefficients value 
Bi = table2array(mdl_nl.Coefficients);
b1 = Bi(1,1);
b2 = Bi(2,1);
b3 = Bi(3,1);
b4 = Bi(4,1);
b5 = Bi(5,1);




%%%
% ==========================================================================
%  Final Estimator: with Model Function of the chosen Type above: 
% ==========================================================================



%yp1 = (b1*exp(MLT_shifted(1:end)).^abs(b2) + (conc_previous.^abs(b3)) + b4*(DUST(1:end))+b5*(SST_nonlinear(1:end))+b6*(conc_previous(1:end).*exp(b7.*(MLT_avgnorm(1:end)-b8)))).*sign(heaviside(b1*exp(MLT_shifted(1:end)).^abs(b2) + (conc_previous.^abs(b3)) + b4*(DUST(1:end))+b5*(SST_nonlinear(1:end))+b6*(conc_previous(1:end).*exp(b7.*(MLT_avgnorm(1:end)-b8)))));
%yp1 = (b1.*(MLT_avgnorm0(1:end)).^abs(b5) + b2*exp(1-MLT_avgnorm1(1:end)).*conc_previous + b3.*(DUST(1:end))+b4*(SST_nonlinear(1:end))).*sign(heaviside(b1.*(MLT_avgnorm0(1:end)).^abs(b5) + b2*exp(1-MLT_avgnorm1(1:end)).*conc_previous + b3.*(DUST(1:end))+b4*(SST_nonlinear(1:end))));
%yp1 = (b1.*(MLT_avgnorm0(1:end)).^abs(b5) + b2*(1-MLT_avgnorm1(1:end).^3).*conc_previous + b3.*(DUST(1:end))+b4*(SST_nonlinear(1:end))).*sign(heaviside(b1.*(MLT_avgnorm0(1:end)).^abs(b5) + b2*exp(1-MLT_avgnorm1(1:end)).*conc_previous + b3.*(DUST(1:end))+b4*(SST_nonlinear(1:end))));
%yp1 = (b1.*(MLT_shifted(1:end)).^abs(b5) + b2*(1-MLT_avgnorm1(1:end).^3).*conc_previous + b3.*(DUST(1:end))+b4*(SST_nonlinear(1:end))).*sign(heaviside(b1.*(MLT_shifted(1:end)).^abs(b5) + b2*(1-MLT_avgnorm1(1:end).^3).*conc_previous + b3.*(DUST(1:end))+b4*(SST_nonlinear(1:end))));
%yp1 = (b1.*(MLT_shifted(1:end)).^abs(round(b5)) + b2*(1-MLT_avgnorm1(1:end).^3).*conc_previous + b3.*(DUST(1:end))+b4*(SST_nonlinear(1:end))).*sign(heaviside(b1.*(MLT_shifted(1:end)).^abs(round(b5)) + b2*(1-MLT_avgnorm1(1:end).^3).*conc_previous + b3.*(DUST(1:end))+b4*(SST_nonlinear(1:end))));
for yptest1=1:length(MLT_shifted)
if MLT_shifted(yptest1) >=0
%yp1(yptest1) = (b1.*(MLT_shifted(yptest1)).^abs(b5) + b2*(1-MLT_avgnorm1(yptest1).^3).*conc_previous(yptest1) + b3.*(DUST(yptest1))+b4*(SST_nonlinear(yptest1))).*sign(heaviside(b1.*(MLT_shifted(yptest1)).^abs(b5) + b2*(1-MLT_avgnorm1(yptest1).^3).*conc_previous(yptest1) + b3.*(DUST(yptest1))+b4*(SST_nonlinear(yptest1))));
yp1(yptest1) = (b1.*(MLT_shifted(yptest1)).^abs(b5) + b2*(1-MLT_avgnorm1(yptest1).^3).*conc_previous(yptest1) + b3.*(DUST(yptest1))+b4*(SST_nonlinear(yptest1))).*sign(heaviside(b1.*(MLT_shifted(yptest1)).^abs(b5) + b2*(1-MLT_avgnorm1(yptest1).^3).*conc_previous(yptest1) + b3.*(DUST(yptest1))+b4*(SST_nonlinear(yptest1))));
else
%yp1(yptest1) = (-b1.*(abs(MLT_shifted(yptest1))).^abs(b5) + b2*(1-MLT_avgnorm1(yptest1).^3).*conc_previous(yptest1) + b3.*(DUST(yptest1))+b4*(SST_nonlinear(yptest1))).*sign(heaviside(-b1.*(abs(MLT_shifted(yptest1))).^abs(b5) + b2*(1-MLT_avgnorm1(yptest1).^3).*conc_previous(yptest1) + b3.*(DUST(yptest1))+b4*(SST_nonlinear(yptest1))));
yp1(yptest1) = (-b1.*(abs(MLT_shifted(yptest1))).^abs(b5) + b2*(1-MLT_avgnorm1(yptest1).^3).*conc_previous(yptest1) + b3.*(DUST(yptest1))+b4*(SST_nonlinear(yptest1))).*sign(heaviside(-b1.*(abs(MLT_shifted(yptest1))).^abs(b5) + b2*(1-MLT_avgnorm1(yptest1).^3).*conc_previous(yptest1) + b3.*(DUST(yptest1))+b4*(SST_nonlinear(yptest1))));
end
end
%yp1 = (b1.*(MLT_avgnorm0(1:end)).^abs(b5) + b2*(1-MLT_avgnorm1(1:end).^3).*conc_previous + b3.*(DUST(1:end))+b4*(SST_nonlinear(1:end))).*sign(heaviside(b1.*(MLT_avgnorm0(1:end)).^abs(b5) + b2*(1-MLT_avgnorm1(1:end).^3).*conc_previous + b3.*(DUST(1:end))+b4*(SST_nonlinear(1:end))));



end
%%
%%%%validation

for yptest1=1:length(MLT_shifted)
if MLT_shifted(yptest1) >=0
%%%% default %%%%
yp1(yptest1) = (b1.*(MLT_shifted(yptest1)).^abs(b5) + b2*(1-MLT_avgnorm1(yptest1).^3).*conc_previous(yptest1) + b3.*(DUST(yptest1))+b4*(SST_nonlinear(yptest1))).*sign(heaviside(b1.*(MLT_shifted(yptest1)).^abs(b5) + b2*(1-MLT_avgnorm1(yptest1).^3).*conc_previous(yptest1) + b3.*(DUST(yptest1))+b4*(SST_nonlinear(yptest1)))); 
%%%% no mlt %%%%
%yp1(yptest1) = (b2.*conc_previous(yptest1) +b3.*(DUST(yptest1))+b4*(SST_nonlinear(yptest1))).*sign(heaviside(b2.*conc_previous(yptest1)+ b3.*(DUST(yptest1))+b4*(SST_nonlinear(yptest1))));
%%%% no self seeding %%%%
%yp1(yptest1) = (b1.*(MLT_shifted(yptest1)).^abs(b5) + b2*(1-MLT_avgnorm1(yptest1).^3) + b3.*(DUST(yptest1))+b4*(SST_nonlinear(yptest1))).*sign(heaviside(b1.*(MLT_shifted(yptest1)).^abs(b5) + b2*(1-MLT_avgnorm1(yptest1).^3) + b3.*(DUST(yptest1))+b4*(SST_nonlinear(yptest1))));
%%%% no dust %%%%
%yp1(yptest1) = (b1.*(MLT_shifted(yptest1)).^abs(b5) + b2*(1-MLT_avgnorm1(yptest1).^3).*conc_previous(yptest1)+b4*(SST_nonlinear(yptest1))).*sign(heaviside(b1.*(MLT_shifted(yptest1)).^abs(b5) + b2*(1-MLT_avgnorm1(yptest1).^3).*conc_previous(yptest1)+b4*(SST_nonlinear(yptest1))));
%%%% no SST %%%%
%yp1(yptest1) = (b1.*(MLT_shifted(yptest1)).^abs(b5) + b2*(1-MLT_avgnorm1(yptest1).^3).*conc_previous(yptest1) + b3.*(DUST(yptest1))).*sign(heaviside(b1.*(MLT_shifted(yptest1)).^abs(b5) + b2*(1-MLT_avgnorm1(yptest1).^3).*conc_previous(yptest1) + b3.*(DUST(yptest1))));
else
%%%% default %%%%
yp1(yptest1) = (-b1.*(abs(MLT_shifted(yptest1))).^abs(b5) + b2*(1-MLT_avgnorm1(yptest1).^3).*conc_previous(yptest1) + b3.*(DUST(yptest1))+b4*(SST_nonlinear(yptest1))).*sign(heaviside(-b1.*(abs(MLT_shifted(yptest1))).^abs(b5) + b2*(1-MLT_avgnorm1(yptest1).^3).*conc_previous(yptest1) + b3.*(DUST(yptest1))+b4*(SST_nonlinear(yptest1))));
%%%% no mlt %%%%
%yp1(yptest1) = (b2.*conc_previous(yptest1) + b3.*(DUST(yptest1))+b4*(SST_nonlinear(yptest1))).*sign(heaviside(b2.*conc_previous(yptest1) + b3.*(DUST(yptest1))+b4*(SST_nonlinear(yptest1))));
%%%% no self seeding %%%%
%yp1(yptest1) = (-b1.*(abs(MLT_shifted(yptest1))).^abs(b5) + b2*(1-MLT_avgnorm1(yptest1).^3) + b3.*(DUST(yptest1))+b4*(SST_nonlinear(yptest1))).*sign(heaviside(-b1.*(abs(MLT_shifted(yptest1))).^abs(b5) + b2*(1-MLT_avgnorm1(yptest1).^3) + b3.*(DUST(yptest1))+b4*(SST_nonlinear(yptest1))));
%%%% no dust %%%%
%yp1(yptest1) = (-b1.*(abs(MLT_shifted(yptest1))).^abs(b5) + b2*(1-MLT_avgnorm1(yptest1).^3).*conc_previous(yptest1)+b4*(SST_nonlinear(yptest1))).*sign(heaviside(-b1.*(abs(MLT_shifted(yptest1))).^abs(b5) + b2*(1-MLT_avgnorm1(yptest1).^3).*conc_previous(yptest1)+b4*(SST_nonlinear(yptest1))));
%%%% no SST %%%%
%yp1(yptest1) = (-b1.*(abs(MLT_shifted(yptest1))).^abs(b5) + b2*(1-MLT_avgnorm1(yptest1).^3).*conc_previous(yptest1) + b3.*(DUST(yptest1))).*sign(heaviside(-b1.*(abs(MLT_shifted(yptest1))).^abs(b5) + b2*(1-MLT_avgnorm1(yptest1).^3).*conc_previous(yptest1) + b3.*(DUST(yptest1))));
end
end

%%

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

RRINDEX=36;
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

a2=area(datevector(a:b),yp1);
a2.FaceColor=[65,174,118]./256;
a2.FaceAlpha = 0.7;
a2.EdgeColor = 'none';
h2=plot(datevector(a:b),yp1,'-o');hold on
h2.Color=[65,174,118]./256;
h2.LineWidth = 1.5;
h2.MarkerSize = 4;
h2.MarkerFaceColor = [65,174,118]./256;



%h2 = bar(datevector(a:b),yp1) % non serve più normalizzare.
%h2.FaceColor = [.5 .5 .5] % if the other is azzurrognolo
%h2.FaceColor = [.3 .3 .3];
%h2.FaceColor = [65,174,118]./256; 
%h2.FaceColor = [5,112,176]./256;
%h2.BarWidth = 0.6
%h2.FaceAlpha = 1.0;
%h2.FaceAlpha = 0.7;
%h2.EdgeColor = 'none';
hold on 
%title('(e)')
h3 = plot(datevector(a:b), obs_webdig, '-o');hold on  % Use '-o' for line with circle markers
% Set line color (gray)
h3.Color = [.3 .3 .3];

% Optional: adjust line width and marker style
h3.LineWidth = 1.5;
h3.MarkerSize = 4;
h3.MarkerFaceColor = [.3 .3 .3];

legend([h3 h2],'Observations','Model','Location','northwest')

ylim([0 25])
set(gca,'FontSize',15); 
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
set(gca,'TickDir','out')
xlabel('Years','Fontsize',15)
ylabel('Sargassum concentrations [Mton]')
%%%ALL R VALUE
text(2025,23,['R=' sprintf('%.2f', corr(obs_webdig, yp1'))],'Fontsize',15)
text(2025, 21, ['MSE = ' sprintf('%.2f', mean((yp1' - obs_webdig).^2))], 'FontSize', 15)
%title( 'E', 'FontSize', 15)


%text(2025,19,['R=' sprintf('%.2f', corr(obs_webdig(end-RRINDEX:end), yp1(end-RRINDEX:end)))],'Fontsize',15,'color','b')
%text(2025, 17, ['MSE = ' sprintf('%.2f', mean((yp1(end-RRINDEX:end) - obs_webdig(end-RRINDEX:end)).^2))], 'FontSize', 15,'color','b')
%title('(c)')

% ===================================================
% PUT THE RIGHT NAME FOR YOUR SAVED FIGURE HERE: 

%print(gcd cf,'../Figures/Fig_regression_model_reducedGASB_giallino.png','-dpng','-r350')

 % ===================================================
%outfile=['/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/xingplot/statistical_anlyz_linear.png'];print([outfile],'-dpng','-r300');
%outfile=['/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/xingplot/factortest/exclude_dust.png'];print([outfile],'-dpng','-r300');
%outfile=['/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/xingplot/statistical_anlyz_control.png'];print([outfile],'-dpng','-r300');
%outfile=['/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/paperfigure/figure_nonlinearmodel/Fig3A.png'];print([outfile],'-dpng','-r300');
%outfile=['/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/utils/newdata/domain/Lat1to15/Regressiontest/Fig3_4_version2.png'];print([outfile],'-dpng','-r300');
%outfile=['/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/utils/newdata/domain/Lat1to11/Fig3_SM.png'];print([outfile],'-dpng','-r300');
%outfile=['/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/utils/newdata/domain/GASB/Fig3_GASB.png'];print([outfile],'-dpng','-r300');
%outfile=['/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/paperfigure/FigS4/FigS4e.png'];print([outfile],'-dpng','-r300');
%outfile=['/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/paperfigure/figure_nonlinearmodel/optimization/integrationcase_yp.png'];print([outfile],'-dpng','-r300');
%%
%%%%calculate each term contribution
for yptest1=1:length(MLT_shifted)
if MLT_shifted(yptest1) >=0
%yp1(yptest1) = (b1.*(MLT_shifted(yptest1)).^abs(b5) + b2*(1-MLT_avgnorm1(yptest1).^3).*conc_previous(yptest1) + b3.*(DUST(yptest1))+b4*(SST_nonlinear(yptest1))).*sign(heaviside(b1.*(MLT_shifted(yptest1)).^abs(b5) + b2*(1-MLT_avgnorm1(yptest1).^3).*conc_previous(yptest1) + b3.*(DUST(yptest1))+b4*(SST_nonlinear(yptest1))));
test1(yptest1)=b1.*(MLT_shifted(yptest1)).^abs(b5);%.*sign(heaviside(b1.*(MLT_shifted(yptest1)).^abs(b5)));
test2(yptest1)=b2*(1-MLT_avgnorm1(yptest1).^3).*conc_previous(yptest1);%.*sign(b2*(1-MLT_avgnorm1(yptest1).^3).*conc_previous(yptest1));
else
%yp1(yptest1) = (-b1.*(abs(MLT_shifted(yptest1))).^abs(b5) + b2*(1-MLT_avgnorm1(yptest1).^3).*conc_previous(yptest1) + b3.*(DUST(yptest1))+b4*(SST_nonlinear(yptest1))).*sign(heaviside(-b1.*(abs(MLT_shifted(yptest1))).^abs(b5) + b2*(1-MLT_avgnorm1(yptest1).^3).*conc_previous(yptest1) + b3.*(DUST(yptest1))+b4*(SST_nonlinear(yptest1))));
test1(yptest1)=-b1.*(abs(MLT_shifted(yptest1))).^abs(b5);%.*sign(heaviside(-b1.*(abs(MLT_shifted(yptest1))).^abs(b5)));
test2(yptest1)=b2*(1-MLT_avgnorm1(yptest1).^3).*conc_previous(yptest1);%.*sign(b2*(1-MLT_avgnorm1(yptest1).^3).*conc_previous(yptest1));

end
end

figure(11)
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
plot(datevector(a:b),test1,'o-','color',[0.01 0.01 0.5],'linewidth',1.5,'MarkerSize',4);hold on
plot(datevector(a:b),test2,'o-','color',[ 0.7461 0.5039 0.1758],'linewidth',1.5,'MarkerSize',4);hold on
set(gca,'FontSize',15); 
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
set(gca,'TickDir','out')
xlabel('Years','Fontsize',15)
ylabel('Sargassum concentrations [Mton]')
ylim([-1 15])
legend('MLD term','Sargassumsphere term','Location','northwest')
%title( 'B', 'FontSize', 15)
%outfile=['/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/paperfigure/figure_nonlinearmodel/Fig3B.png'];print([outfile],'-dpng','-r300');

