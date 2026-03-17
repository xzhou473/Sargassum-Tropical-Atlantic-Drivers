%%
clear all
load('/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/paperfigure/figure_nutrient/nutrientplotdata.mat');

%load('/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/paperfigure/figure_nutrient/nutrients_ljuba/nutrients_codes_and_output/N_Vol_in_Mton_box_spacesum_Lat_1-15_lag0.mat');
load('/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/paperfigure/figure_nutrient/nutrients_ljuba/nutrients_codes_and_output/N_Vol_in_Mton_box_spacesum_bloomcontour_lag0.mat');
%%
figure(12)
clf
scale2=80;
scale=1; 
set(gcf,'position', [474   157     12*scale2 6*scale2])
set(gcf,'paperposition', [0   0   12*scale  6*scale])

% Set percentage scaling
perc_N = 1; % 100%

% Define colors
%bar_color0
bar_color = [31,120,180] / 256;
line_color =  [65, 174, 118] / 256; % black

% Left Y-axis for bar plot
yyaxis left
a2=area(datevector_MLT, circshift(N_Vol_in_Mton_box_spacesum, 3));hold on
a2.FaceColor=bar_color;
a2.FaceAlpha = 0.7;
a2.EdgeColor = 'none';
h2=plot(datevector_MLT, circshift(N_Vol_in_Mton_box_spacesum, 3),'-o');hold on
h2.Color=bar_color;
h2.LineWidth = 1.5;
h2.MarkerSize = 4;
h2.MarkerFaceColor = bar_color;


%h2 = bar(datevector_MLT, N_Vol_in_Mton_box_spacesum, 0.6); % 0.6 is BarWidth
%h2.FaceColor = bar_color;
%h2.EdgeColor = 'none';
%h2.FaceAlpha = 0.7;
ylabel('Additional N input (Mton)')

% Set the left axis color to match bar color
ax = gca;
ax.YColor = bar_color;

hold on
ylim([0 0.08])
% Right Y-axis for line plot
yyaxis right
plot(observation_dates, obs_webdig * perc_N, 'o-', ...
    'Color', line_color, 'LineWidth', 1.5, ...
    'MarkerSize', 2.5, 'MarkerFaceColor', line_color);
ylabel('Sargassum concentrations (Mton)')

% Set the right axis color to match line color
ax.YAxis(2).Color = line_color;

% X-axis label
set(gca,'fontsize',10)
set(gca,'Tickdir','out')
xlabel('Years')
%outfile=['/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/paperfigure/figure_nutrient/Fig4_GASB_lag3months.png'];print([outfile],'-dpng','-r300');