
clear all
close all
clc

%addpath('/Users/ljubanovi/Documents/Matlab_Additional_Toolboxes/github_repo_cmocean/')
%addpath('/Users/ljubanovi/Documents/Matlab_Additional_Toolboxes/github_repo_cmocean/')
%addpath(genpath('/Users/lnovi3/Documents/additional_matlab_packages_local'));

%addpath('~/GT/CoralTri/scripts')  
%addpath('~/Matlab_Additional_Toolboxes/EOF_v1.1/EOF/')
%addpath('~/Matlab_Additional_Toolboxes/plotpsd/plotpsd')
%addpath('~/Matlab_Additional_Toolboxes/ColorBrewer_colormaps/github_repo') % Please cite it if you use it! (see .txt file in folder )
%addpath('/Users/ljubanovi/Documents/Matlab_Additional_Toolboxes/github_repo_cmocean/')
%addpath('/Users/ljubanovi/Documents/Matlab_Additional_Toolboxes/rgbmap_v2') % Cite uf you this to do custom colormaps
%addpath('/Users/ljubanovi/Documents/Matlab_Additional_Toolboxes')
%addpath('/Users/lnovi3/Documents/additional_matlab_packages_local') % mac
%addpath('/Users/ljubanovi/Documents/Matlab_Additional_Toolboxes/xkcd_rgb_v1/XKDC_RGB/');
%addpath(genpath('/Users/ljubanovi/Documents/Matlab_Additional_Toolboxes'))



%% Part 1: compute the N necessary to sustain sargassum biomass:

% To do this, I digitized Figure 1d of Bach et al 2021, using
% WebPlotDogitizer online version4. 
% To extract the data from Fig 1d of Bach el at 2021,
% I have used the manual extraction option of the software WebPlotDigitizer:
% (Rohatgi, 2024, https://automeris.io/WebPlotDigitizer.html, WebPlotDigitizer version 4.7, 2024), 
% available at https://automeris.io/WebPlotDigitizer.html; 
% -------------------------------------------------------------------------------

Ctab = csvread('../utils/Fig1_d_Bach_et_al_2021_MillionTon.csv');

C = Ctab(:,2);
C = unique(C,'stable');

x = Ctab(:,1);
x = unique(x,'stable');


x1=1:96; % months from Jan 2011 to Dec 2018 because their figure spans this period.


C_2011_2018_Mton = interp1(x,C,x1);

datevector_C = datetime(2011,1:length(x1),1) % Create date-timevectore form Jan 2011 to Dec 2018 as Fig 1d in Bach et al 2021 



%% Plot necessary C and N 

figure
set(gcf,'units','centimeters','paperunits','centimeters');
setfigsize(gcf,25,10) % 30 15 
box on 
xlabel('','Interpreter','latex');
ylabel('C [Mton]','Interpreter','latex');
title('POC Bach et al. 2021','Interpreter','latex');
grid on 
box on 
set(gca,'GridLineStyle',':')
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')

set(gca,'FontSize',18,'FontName','Times');
hold on 
h2 = bar(datevector_C,C_2011_2018_Mton) 

h2.FaceColor = [.5 .5 .5]
h2.BarWidth = 0.5
h2.FaceAlpha = 1

print(gcf,'../Figures/POC_estimated_from_Bach_et_al_2021.png','-dpng','-r350')

%% ===================================================

% Compute the corresponding necessary N, using the percentages of C and N
% of Fig 2a (2010s) in Lapointe et al 2021 https://doi.org/10.1038/s41467-021-23135-7

perc_C  = 0.3   % circa il 30 per cento della massa di sarg, è fatta da C
perc_N = 0.012  % circa il 1.2 per cento della massa di sarg, è fatta da N

massa_sarg = C_2011_2018_Mton./perc_C;

N_2011_2018_Mton = massa_sarg.*perc_N;


figure 
set(gcf,'units','centimeters','paperunits','centimeters');
setfigsize(gcf,25,10) % 30 15 
box on 
xlabel('','Interpreter','latex');
ylabel('N [Mton]','Interpreter','latex');

title('Necessary N (estimated from POC - Bach et al. 2021)','Interpreter','latex');

grid on 
box on 
set(gca,'GridLineStyle',':')
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')

set(gca,'FontSize',18,'FontName','Times');
hold on 
h2 = bar(datevector_C,N_2011_2018_Mton) 
h2.FaceColor = [.5 .5 .5]
h2.BarWidth = 0.5
h2.FaceAlpha = 1

print(gcf,'../Figures/Necessary_N_allyears.png','-dpng','-r350')



%% Part 2: compute the N that we have using MLT and NO3:  

% load the monthly MLT after bloom (Deeper, monthly Z2). It's a 2D*Time. 

 ncid = netcdf.open('../utils/nutriflux_new_bloomcontour/MLT_2011-2022_f64_025_masked_bloomcontour.nc')  
  lon = double(netcdf.getVar(ncid,1));   %1 Always cross-check if it's 1
  lat = double(netcdf.getVar(ncid,2));   %2 Always cross-check if it's 2
  z2 = double(netcdf.getVar(ncid,3));    %3 always check

netcdf.close(ncid);

% load the ymonmean MLT pre-bloom (Shallower, monthly climatology Z1) - i soliti 12 mesi ripetuti
% ogni anno per 12 anni, vedi sh script. 
ncid = netcdf.open('../utils/nutriflux_new_bloomcontour/ymonmean_mlt_1999-2010_f64_025_masked_repeat_bloomcontour.nc') 
  lon = double(netcdf.getVar(ncid,2));   %1 Always cross-check if it's 1
  lat = double(netcdf.getVar(ncid,3));   %2 Always cross-check if it's 2
  z1 = double(netcdf.getVar(ncid,4));    %3 always check
netcdf.close(ncid);

% Load monthly climatology Nutrients (NO3) before bloom - i soliti 12 mesi ripetuti
% ogni anno per 12 anni, vedi sh script. 
ncid = netcdf.open('../utils/nutriflux_new_bloomcontour/ymonmean_NUTRI_1999-2010_repeat_bloomcontour.nc') % *** LN: QUESTI SONO DA RICONTROLLARE, E CAPIRE SE DEVONO ESSERE DI PERIODO CONSISTENTE CON DATI MLT OPPURE NO...
  lon = double(netcdf.getVar(ncid,2));    % 1 Always cross-check if it's 1
  lat = double(netcdf.getVar(ncid,3));    % 2 Always cross-check if it's 2
  depth_nutri = (netcdf.getVar(ncid,4));  % 3 always check
  no3 = double(netcdf.getVar(ncid,5));    % 4 always check
netcdf.close(ncid);

% Load grid cell area in m2 computed with cdo (cdo gridarea Z2_yseasm_2020_2022_glorys_mlt_GASB_masked_025remapped.nc gridcell_area_025.nc)

ncid = netcdf.open('../utils/nutriflux_new_bloomcontour/gridcell_area_new_025.nc')
  lon = double(netcdf.getVar(ncid,0));    % Always cross-check if it's 1
  lat = double(netcdf.getVar(ncid,1));    % Always cross-check if it's 2
  cell_area_m2 = (netcdf.getVar(ncid,2)); % always check
netcdf.close(ncid);


[X Y]=meshgrid(lon,lat);
X=X';
Y=Y';


%% Remove land masks
land1 = max(max(max(max(max(no3)))))        % Land no3
land2 = min(min(min(z1)))                   % Land mlt


no3(no3 == land1) = NaN;

z1(z1 == land2) = NaN;
z2(z2 == land2) = NaN;




%% Per ogni mese nel post-bloom calcola l'integrale verticale dei Nutrienti prebloom, 
%% compresi tra z2 di quel mese e z1 (climatology) del mese corrispondente 

% 1. In ogni mese trova la depth di no3 piu vicina a z1 e a z2, in ogni
% punto. 

% 2. Poi prendi, in ogni punto (e in ogni istante), solo i no3 PRE-BLOOM compresi tra quelle due depths, 
% e metti a NaN tutti gli altri sulla verticale. Qui sei ancora in mmol\m3 

% 3. in ogni punto (i.e. in ogni grid point) moltiplica per l'area del grid point (m2) 
% e per lo spessore del grid point (m). Adesso sei in mmol. 

% 4. Next, somma quel valore in mmol in verticale in ogni punto tra z1 e
% z2. E questo è l'integrale di volume che vuoi (una mappa 2D, punto a punto, e sei sempre in mmol perchè sommi mmol).

% 5. Infine, vai in mol poi trova le moli di N da quelle di NO3, e poi vai in grammi. 
% 5.1 - SE USI LA BOX: definisci una box orizzontale che comprenda l'area di
% interesse e per ogni mese fai la media spaziale su quell'area, e poi vai
% in Mton. Quindi hai delle bars mensili di N in MTon (da confrontare col N
% necessario credo).
% 5.2 - SE NON USI LA BOX MA USI IL BLOOM CONTOUR: same as point 5.1, ma la
% media spaziale la fai solo sul bloomcontour. 

%% 

% 1. In each grid poin of the nutrients matrices, find the depths closest to the value of z1 and z2 in that grid point 

dim1 = length(no3(:,1,1,1))  % lon
dim2 = length(no3(1,:,1,1))  % lat
dim3 = length(no3(1,1,:,1))  % depth 
dim4 = length(no3(1,1,1,:))  % time longer (144)

datevector_MLT = datetime(2011,1:dim4,1); % POST-BLOOM: Create date-timevectore form Jan 2011 for the entire length of MLT data we have (post bloom here).
datevector_N = datetime(1999,1:dim4,1);   % Nutri pre-bloom start in Jan 1999. Create date-timevectore for the entire length of MLT data we have, because we will have to use it together with the mlt postbloom

%%

% Initialize dummy arrays 
no3_sum_nearest(1:dim1,1:dim2,1:dim4) = NaN;   % dummy, to be filled with vertical integral at each month. Dimensions x,y,time (no more vertical dim)

no3_at_Z1_nearest(1:dim1,1:dim2,1:dim4) = NaN; % dummy, to be filled with no3 at Z1, at each month
no3_at_Z2_nearest(1:dim1,1:dim2,1:dim4) = NaN; % dummy, to be filled with no3 at Z2, at each month

depth_z1(1:dim1,1:dim2,1:dim4) = NaN;      
depth_z2(1:dim1,1:dim2,1:dim4) = NaN;



%% 

for time = 1:dim4  % time dimension (the longer one, 12 years)
    time
    for i = 1:dim1;  % lon dimension

        for j = 1:dim2;  % lat dimension

            for k = 1:dim3;  % depth dimensiomn
           
        
                z1_current_point = z1(i,j,time);        % Depth Z1 a cui vorrei i nutri in questo grid point i,j e in questo istante. It is a number, i.e. a point value.
            
                z2_current_point = z2(i,j,time);        % Depth Z2 a cui vorrei i nutri in questo grid point i,j e in questo istante. It is a number, i.e. a point value. 
            
               
                if  isnan(z1_current_point) | isnan(z2_current_point)    % If at least one is NaN, leave it NaN everywhere because it's a land point
              
                    no3_nearest(i,j,time) = NaN;
       
                else 
               
                % First, find the value in depth_nutri that is nearest to the depth a cui li vorrei: 
            
                mydiff_z1 = depth_nutri -  z1_current_point;
                mydiff_z2 = depth_nutri -  z2_current_point;
                
                % Trovo l'indice a cui si trova la depth piu' vicina al
                % current z1 e z2 (quindi l'indice dove ho la differenza più piccola in valore assouto tra la depth nutri e quella che vorrei)
 
                % Prima devo capire che segno va bene per z1 (uno dei due sarà empty): 
                idx1_z1 =  find(depth_nutri == (z1_current_point - min(abs(mydiff_z1))));
                idx2_z1 =  find(depth_nutri == (z1_current_point + min(abs(mydiff_z1)))); 
          

                % Prima devo capire che segno va bene per z2 (uno dei due sarà empty): 
                idx1_z2 =  find(depth_nutri == (z2_current_point - min(abs(mydiff_z2))));
                idx2_z2 =  find(depth_nutri == (z2_current_point + min(abs(mydiff_z2)))); 
           
         
                
                  % --------------------------------------------------
                    % Z1:
                    idx_z1 = unique([idx1_z1 idx2_z1]);                    % DEPTH IDX AT WHICH I HAVE NEAREST DEPTH to Z1(scelgo l'unico non empty - l'altro e' empty sicuramente, ma non lo so a priori). 
           
                    nearest_depth_current_z1 = depth_nutri(idx_z1);        % DEPTH NEAREST to Z1 in current grid point 
                    
                    depth_z1(i,j,time) = nearest_depth_current_z1;         % This is the Z1 matrix (mlt prima)
                    
                    % Extract Nutrient ONLY at this depth Z1 (nearest depth), in mmol/m3: 
                    no3_at_Z1_nearest(i,j,time) = no3(i,j,idx_z1,time); 
                    
                    % Z2:
                    idx_z2 = unique([idx1_z2 idx2_z2]);               % DEPTH IDX AT WHICH I HAVE NEAREST DEPTH to Z2(scelgo l'unico non empty - l'altro e' empty sicuramente, ma non lo so a priori). 
            
                    nearest_depth_current_z2 = depth_nutri(idx_z2);   % DEPTH NEAREST to Z2
                    
                    depth_z2(i,j,time) = nearest_depth_current_z2;         % This is the Z2 matrix (mlt dopo)
                    
                    % Extract Nutrient ONLY at this depth Z2 (nearest depth), in mmol/m3: 
                    no3_at_Z2_nearest(i,j,time) = no3(i,j,idx_z2,time); 
                  % --------------------------------------------------
       
              
                
              end  % close the if-statement loop 
            
          end  % close the loop on k 
    
       end % close the loop on j 
    
    end  % close the loop on i
end

disp ' loop 1 ended'

%%
   
% in ogni strato verticale per ogni istante in ogni punto voglio dA*dZ*dNutri e poi sommo tutto in verticale tra quelle due depths appena trovate punto a punto (e divento in mmol visto che dA*dZ è m2*m)

dim_vert = length(depth_nutri);                 % Total length of vertical dimension of the nutrient array. Depth nutri is a vector, extracted from the netcdf file. It contains the vertical levels at which the nutients were stored in the data, as they came from the provider.

dVolume(1:dim1,1:dim2,1:dim_vert,1:dim4) = NaN; % Initialize
dNO3vol(1:dim1,1:dim2,1:dim_vert,1:dim4) = NaN; % Initialize

dZ(1:dim3) = NaN;                        % Initialize. dim3 is the total lenght of the vertical dimension of the nutrient array (should be 29)
dZ(1) = depth_nutri(1);                  % Initialize, start from the shallowest depth of the nutrient array. So far I haven't excluded any depths yet.

for time = 1:dim4
    for k = 2:dim_vert   % This consider all the depths in the nutri files, not only those between z1 and z2. Then I will select those ones later, befor summing in the vertical. 
         for i = 1:dim1
             for j = 1:dim2 
            
                 dZ(k) = depth_nutri(k)-depth_nutri(k-1);                    % dZ in m, con dZ = Zdeeper - Zshallower. 
            
                 dVolume(i,j,k,time) = cell_area_m2(i,j).*dZ(k);             % dVol in m3
            
                 dNO3vol(i,j,k,time) = no3(i,j,k,time)*dVolume(i,j,k,time);  % [mmol]: perche' dNO3 at kth depth is in in mmol/m3, so I multiply dNO3*dVol and I have dNO3vol in mmol (in each grid point, i.e. it's not yet a vertical sum)
            
             end 
         end
   
    end
end

% Now compute the vertical integral (i.e. sum) of NO3 between Z1 and Z2 (and convert to [mol]). 
% To this aim, we only consider the points where the MLT got deeper, and
% put 0 elsewhere. And of course we only consider depths between Z1 and Z2.

NO3vol_mol(1:dim1,1:dim2,1:dim4) = NaN; % Initialize


DELTA_Z = depth_z2-depth_z1;            % Il mlt ha fatto deepening solo dove questa roba e' > 0. Per cui considero in horizontal solo questi punti, gli altri metto NaN. 

mask_deepening = DELTA_Z;
mask_deepening(DELTA_Z <= 0 ) = NaN;
mask_deepening(DELTA_Z > 0 ) = 1;
%%

dNO3vol_deepening(1:dim1,1:dim2,1:dim3,1:dim4) = NaN;   % Initialize

for i = 1:dim3 
    i
    dNO3vol_deepening(:,:,i,:) = squeeze(dNO3vol(:,:,i,:)).*mask_deepening;  % Ha NaN dove MLT non ha fatto deepening
end


%% Devo sommare gli strati verticali di dNO3vol_deepening, tra z1 e z2 in ogni grid point. 

%
% I need to select only the depths between z1 and z2. 
%
for time =1:dim4
    for i = 1:dim1
        for j = 1:dim2 
       
            current_z1 = depth_z1(i,j,time);  % shallower
            current_z2 = depth_z2(i,j,time);  % deeper
          
            idx_z1 = find(depth_nutri == current_z1);
            idx_z2 = find(depth_nutri == current_z2);
         
            if (current_z1 ~= NaN) & (current_z2~= NaN) % sommo in vertical solo dove non ho NaN

               % Considero solo i puti dove il MLT ha fatto deepening quindi
               % sommo in verticale (cioè lungo la dimensione 3) sul dNO3vol_deepening tra z1 e z2
    
                 NO3vol_mol(i,j,time) = 10^-3*nansum(dNO3vol_deepening(i,j,idx_z1:idx_z2,time),3);    % [mol] perche' moltiplico per 10^-3
            else 
                 NO3vol_mol(i,j,time) = NaN; % se c'è nan, non sommo, resta nan.
            end
        end 
    end
end

% And put NaN over land points

landmask = (squeeze(depth_z2(:,:,1)));
landmask = landmask - landmask;
landmask = landmask + 1;   % Mask che ha NaN sul land e 1 elsewhere, e ha dimensioni 2D dim1*dim2

%NO3vol_mol(isnan(depth_z2)) = NaN;  
NO3vol_mol = NO3vol_mol.*landmask; % [mol]

% Quante moli di N ci sono in questo NO3? NO3 e' espresso in mol 
% Siccome in NO3 c'è un solo atomo di N, ogni mole di NO3 contiene una mole di N


N_Vol_in_mol = (NO3vol_mol);%.*FACTOR;   %



%% This block if you Select a sub-box 

% If you select a rectangular box:

  lon1 = -89; % -89
  lon2 = -15; % -15
  lat1 = 1; % 1
  lat2 = 15; % 15 reduced box to avoid Sargassum Sea

Xbox = X;
Ybox = Y;

Xbox(Xbox<lon1) = NaN;
Xbox(Xbox>lon2) = NaN;

Ybox(Ybox<lat1) = NaN;
Ybox(Ybox>lat2) = NaN;

xx = Xbox;
yy = Ybox;
xx(~isnan(Xbox)) = 1;
yy(~isnan(Ybox)) = 1;

%MYBOX = xx.*yy; % Uncomment this to use the box instead

%% Otherwise uncomment this line below to use the bloom contour directly instead of the box
MYBOX = 1 % This will deactivate the box selection as now we have the bloomcontours in the MLT and nans elsewhere

% Put a NaN to the nutrient flux matrices if a grid point is outside the box
% N 
N_Vol_in_mol_box = N_Vol_in_mol.*MYBOX;


% Also in the mlt because I want to see the time series here too
depth_z2_box = depth_z2.*MYBOX;
depth_z1_box = depth_z1.*MYBOX;


% and compute the spatial average of the MLT (i.e. of these z1 and z2) over
% the box or the bloom contour (depending which one you selected above):
depth_z2_box_r = reshape(depth_z2_box,[dim1*dim2,dim4]);
depth_z1_box_r = reshape(depth_z1_box,[dim1*dim2,dim4]);

depth_z2_box_mean = nanmean(depth_z2_box_r,1);
depth_z1_box_mean = nanmean(depth_z1_box_r,1);

clear depth_z1_box_r depth_z2_box_r

% Verify what happens at NO3 (not N) at Z2 and Z1 over the box in 2012... is there a peak or sth? 

no3_at_Z2_nearest_box = no3_at_Z2_nearest.*MYBOX;
no3_at_Z1_nearest_box = no3_at_Z1_nearest.*MYBOX;

no3_at_Z2_nearest_r = reshape(no3_at_Z2_nearest_box,[dim1*dim2,dim4]);
no3_at_Z1_nearest_r = reshape(no3_at_Z1_nearest_box,[dim1*dim2,dim4]);

no3_at_Z2_nearest_mean = nanmean(no3_at_Z2_nearest_r,1);
no3_at_Z1_nearest_mean = nanmean(no3_at_Z1_nearest_r,1);


% For each time step, compute the spatially-averaged nutrient flux over the box (i.e. mean over all the grid points) 

N_Vol_in_Mton_box_spaceavg(1:dim4) = NaN; % Initialize 
N_Vol_in_Mton_box_spacesum(1:dim4) = NaN; % Initialize 

for time = 1:dim4 ;
    time

    % reshape and remove NaNs entries

    currentN = squeeze(N_Vol_in_mol_box(:,:,time));

    r_n = reshape(currentN,[1,dim1*dim2]);

    r_n(isnan(r_n)) = [];

    SpaceMean_N_box = nanmean(r_n) % spatial mean over the box, of integrated N  in mol
    SpaceSum_N_box = nansum(r_n)  % spatial sum over the box, of integrated N  in mol


    % Converto in grammi di N

    % Massa molare di N = 14 gr/mol 
    SpaceMean_N_box_gr = 14 * SpaceMean_N_box      % Ora in Grammi
    SpaceSum_N_box_gr = 14 * SpaceSum_N_box        % Ora in Grammi
    
    N_Vol_in_Mton_box_spaceavg(time) = SpaceMean_N_box_gr*(10^-12);  % Million Ton

    N_Vol_in_Mton_box_spacesum(time) = SpaceSum_N_box_gr*(10^-12);  % Million Ton
end

disp ' last loop ended'




%% Prepare stuff for the final n figure: 
%% TO OVERPLOT SARGASSUM N OBTAINED BY OBSERVED BIOMASS: 
% Load Web-Digitized observations: MUST CITE DATA SOURCE AND SOFTWARE USED.
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

obs_webdig = TAB(:,2); % Observed Sarg Biomass in Million Tons

% It's impossible to have negative concentration. Manual digitizing might
% introduce spurious negative values when concentrations are very small. 
% Hence, put to zero every negative point.

obs_webdig(obs_webdig<0) = 0;


% Create date vector for digitized observations:  Jan 2011 to Dec 2022
observation_dates = datetime(2011,1:length(obs_webdig),1) % Centered at the 1st of the month to be consistent with the other timeseries of this script, which is fine as they are always monthly data so no matter the day.




%% FINAL FIGURE: Plot the N that we estimated as a time series - FIGURE FOR THE MANUSCRIPT: 
% ==================================================================================
figure
set(gcf,'units','centimeters','paperunits','centimeters');
setfigsize(gcf,25,10) % 30 15 
box on 
xlabel('','Interpreter','latex');
ylabel('N [Mton]','Interpreter','latex');
%title('estimated N from mlt','Interpreter','latex');
grid on 
box on 
set(gca,'GridLineStyle',':')
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')

set(gca,'FontSize',18,'FontName','Times');
hold on 
%h2 = bar(datevector_MLT,N_Vol_in_Mton_box_spaceavg) % SPACE AVERAGE 
h2 = bar(datevector_MLT,N_Vol_in_Mton_box_spacesum)  % SPACE SUM 

h2.FaceColor = [0.9 0.75 0.]%[.4 .4 .4]
h2.BarWidth = 0.75
h2.FaceAlpha = 0.7 %1

% Overplot: now take the time series of sargassum and transform those into N 
% (1.2% of sargassum is made of N) and add that on top as a line though 
% (so have the bars and then plot in a connected way the Sarg concentrations multiplied by 0.012)

hold on 

plot(observation_dates,obs_webdig*perc_N,'o-','Color','k','LineWidth',0.9,'MarkerSize',2.5,'MarkerFaceColor','k');

ylim([0 0.26])
print(gcf,'../Figures/NutriFlux_allyears_nolag.png','-dpng','-r350')

% ==============================================================================






%% Save the final variable obtained: 
%save('N_Vol_in_Mton_box_spacesum_bloomcontour.mat','N_Vol_in_Mton_box_spacesum')

%% Save the final workspace: 
save Nutrients_workspace.mat

