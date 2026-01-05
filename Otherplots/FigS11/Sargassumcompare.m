%%
clear all
load('/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/srg_rev_for_Xing/srg_rev_for_Xing/digitized_from_wang_et_al_2019_for_revision.mat')

TAB = csvread('../utils/digitized_sarg_biomass/WebDigitized_Jan2011-Dec2022_from_Jouanno-et-al-2023.csv');

obs_webdig = TAB(:,2); 

% It's impossible to have negative concentration. Manual digitizing might
% introduce spurious negative values when concentrations are very small. 
% Hence, put to zero every negative point.

obs_webdig(obs_webdig<0) = 0;


% Create date vector for digitized observations:  Jan 2011 to Dec 2022
observation_dates = datetime(2011,1:length(obs_webdig),15)

%%
figure(22)
clf
plot(datevector_wang,digitized_wang_2019,'-o','Color',[0.85 0.75 0.45],'LineWidth',1.5,'MarkerSize',4,'MarkerFaceColor',[0.85 0.75 0.45]);hold on
plot(observation_dates,obs_webdig,'-o','Color',[0.3 0.3 0.3],'LineWidth',1.5,'MarkerSize',4,'MarkerFaceColor',[0.3 0.3 0.3]);hold on
legend('Adapted from Wang et al., (2019)','Adapted from Jouanno et al., (2023)','Location','Northwest')
set(gca,'fontsize',13)
set(gca,'Tickdir','out')
ylim([0 25])
ylabel('Sargassum concentrations [Mton]')
outfile=['/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/srg_rev_for_Xing/srg_rev_for_Xing/Sargassumcompare.png'];print([outfile],'-dpng','-r300');