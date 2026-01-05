%%
%%%%
%%%%this is a monte carlo simulation for adding error bar to Fig.3, with
%%%%given uncertainty to MLD, Dust and SST
%%%%running this program require the input of NLregression model

%%
%%%calculated use sigma_MC
std_MLD=1.2567;%0.89;%1.8141;%4.8868;
std_Dust=0.0;%0.0192;%0.0182;
std_SST=0.04;%0.2308;%0.3043;

Nsamples=1000;

yp_MC=zeros(Nsamples,length(MLT_shifted))

for i=1:Nsamples
noise_MLT=std_MLD*randn(1,length(MLT))'+MLT;
noise_DUST=std_Dust*randn(1,length(DUST))'+DUST;
noise_SST=std_SST*randn(1,length(SST_nonlinear))'+SST_nonlinear;


noise_MLT_shifted=circshift(noise_MLT,MESI_mlt);
clear MLT_seasonal1
noise_MLT_seasonal1(1:b_dust-a_dust+1) = NaN;
for k = 1:(length(noise_MLT_seasonal1)/3)-1
    noise_MLT_seasonal1(3*k:(3*k)+2) = mean(noise_MLT(3*k:(3*k)+2));
end
noise_MLT_seasonal1(1:2) = mean(noise_MLT(1:2));
noise_MLT_seasonal1(end) = mean(noise_MLT(end));
noise_MLT_seasonal1 = noise_MLT_seasonal1';
noise_MLT_avg = circshift(noise_MLT_seasonal1,toty); 
noise_MLT_avgnorm1=(noise_MLT_avg-min(noise_MLT_avg))./(max(noise_MLT_avg)-min(noise_MLT_avg));


for yptest1=1:length(MLT_shifted)
if noise_MLT_shifted(yptest1) >=0
yp_MC(i,yptest1) = (b1.*(noise_MLT_shifted(yptest1)).^abs(b5) + b2*(1-noise_MLT_avgnorm1(yptest1).^3).*conc_previous(yptest1) + b3.*(noise_DUST(yptest1))+b4*(noise_SST(yptest1))).*sign(heaviside(b1.*(noise_MLT_shifted(yptest1)).^abs(b5) + b2*(1-noise_MLT_avgnorm1(yptest1).^3).*conc_previous(yptest1) + b3.*(noise_DUST(yptest1))+b4*(noise_SST(yptest1))));
else
yp_MC(i,yptest1) = (-b1.*(abs(noise_MLT_shifted(yptest1))).^abs(b5) + b2*(1-noise_MLT_avgnorm1(yptest1).^3).*conc_previous(yptest1) + b3.*(noise_DUST(yptest1))+b4*(noise_SST(yptest1))).*sign(heaviside(-b1.*(abs(noise_MLT_shifted(yptest1))).^abs(b5) + b2*(1-noise_MLT_avgnorm1(yptest1).^3).*conc_previous(yptest1) + b3.*(noise_DUST(yptest1))+b4*(noise_SST(yptest1))));
end
end



end



%%

for yptest1=1:length(MLT_shifted)
if MLT_shifted(yptest1) >=0
yp1(yptest1) = (b1.*(MLT_shifted(yptest1)).^abs(b5) + b2*(1-MLT_avgnorm1(yptest1).^3).*conc_previous(yptest1) + b3.*(DUST(yptest1))+b4*(SST_nonlinear(yptest1))).*sign(heaviside(b1.*(MLT_shifted(yptest1)).^abs(b5) + b2*(1-MLT_avgnorm1(yptest1).^3).*conc_previous(yptest1) + b3.*(DUST(yptest1))+b4*(SST_nonlinear(yptest1))));
else
yp1(yptest1) = (-b1.*(abs(MLT_shifted(yptest1))).^abs(b5) + b2*(1-MLT_avgnorm1(yptest1).^3).*conc_previous(yptest1) + b3.*(DUST(yptest1))+b4*(SST_nonlinear(yptest1))).*sign(heaviside(-b1.*(abs(MLT_shifted(yptest1))).^abs(b5) + b2*(1-MLT_avgnorm1(yptest1).^3).*conc_previous(yptest1) + b3.*(DUST(yptest1))+b4*(SST_nonlinear(yptest1))));
end
end

figure(23)
clf
scale2=80;
scale=1; 
set(gcf,'position', [474   157     12*scale2 6*scale2])
set(gcf,'paperposition', [0   0   12*scale  6*scale])

set(gcf,'units','centimeters','paperunits','centimeters');
box on 
grid on 
box on 
set(gca,'GridLineStyle',':')

RRINDEX=36;

hold on

a2=area(datevector(a:b),mean(yp_MC,1));
a2.FaceColor=[65,174,118]./256;
a2.FaceAlpha = 0.7;
a2.EdgeColor = 'none';
h2=plot(datevector(a:b),mean(yp_MC,1),'-o');hold on
h2.Color=[65,174,118]./256;
h2.LineWidth = 1.5;
h2.MarkerSize = 4;
h2.MarkerFaceColor = [65,174,118]./256;
errorbar(datevector(a:b), mean(yp_MC,1), std(yp_MC,1), 'LineStyle', 'none', 'Color', [65,174,118]./256, 'LineWidth', 2.0);


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

text(2025,23,['R=' sprintf('%.2f', corr(obs_webdig, mean(yp_MC,1)'))],'Fontsize',15)
text(2025, 21, ['MSE = ' sprintf('%.2f', mean((mean(yp_MC,1)' - obs_webdig).^2))], 'FontSize', 15)
%outfile=['/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/paperfigure/figure_nonlinearmodel/Fig3_MC.png'];print([outfile],'-dpng','-r300');

%%
%%%%%residual test
res=obs_webdig-yp1'

%%%%autocorrelation short memory, not long memory
figure(11)
clf
autocorr(res,'NumLags',24)
ylabel('ACF')
xlabel('Lag')
set(gca,'fontsize',13)
xlabel('Month')
%outfile=['/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/paperfigure/figure_nonlinearmodel/Fig3_ACF_exp.png'];print([outfile],'-dpng','-r300');
%%
%%%%
[pval,bp] = archtest(res), %%%%bp test show heteroscedasticity, residual change influenced by dependent variable changes

%res_index=(obs_webdig>5);
%figure(12)
%clf
%plot(yp1(res_index)',res(res_index),'r.');hold on
%plot(yp1(~res_index)',res(~res_index),'b.');hold on
%xlabel('Fitted values')
%ylabel('Residuals')
%yline(0,'k--')
%%
%%%%whether this is normal distribution
qqplot(res)
set(gca,'fontsize',13)
title('B')
%outfile=['/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/paperfigure/figure_nonlinearmodel/Fig3_QQplot.png'];print([outfile],'-dpng','-r300');
outfile=['/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/paperfigure/figure_nonlinearmodel/Fig3_QQplot_b6constant.png'];print([outfile],'-dpng','-r300');
sk=skewness(res) %more large positive residuals than normal distribution
ku=kurtosis(res) %underpredicts extreme event, residuals become more postiive during extreme event.

