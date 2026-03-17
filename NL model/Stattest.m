%%
%%%%
%%%%this script include: 
%%%%	1.monte carlo simulation for adding error bar to Fig.3, with
%%%%	given uncertainty to MLD, Dust and SST
%%%%    2. ACF and QQ plot shown in Fig.S11
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
for i=1:1000
res_MC(i,:)=obs_webdig-yp_MC(i,:)';
end

res_mean = mean(res_MC,1);
res_std  = std(res_MC,0,1);
res_min  = min(res_MC,[],1);
res_max  = max(res_MC,[],1);
n = length(res_mean);
res_sorted = sort(res_mean);
p = ((1:n)-0.5)/n;
q_theory = norminv(p, mean(res_mean), std(res_mean));

%%%std
if true
res_upper_std = res_mean + res_std;
res_lower_std = res_mean - res_std;
res_upper_std_sorted = sort(res_upper_std);
res_lower_std_sorted = sort(res_lower_std);
figure(11)
clf
hold on

fill([q_theory fliplr(q_theory)], ...
     [res_lower_std_sorted fliplr(res_upper_std_sorted)], ...
     [0.8 0.8 0.8], ...
     'EdgeColor','none', ...
     'FaceAlpha',0.5);
% 1:1 reference line
plot(q_theory, q_theory, 'r-', 'LineWidth',1.5)
plot(q_theory, res_sorted, 'kx', 'MarkerSize',6, 'MarkerFaceColor','b')




set(gca,'FontSize',13)
xlabel('Standard Normal Quantiles')
ylabel('Quality of Input Sample')

% title('QQ Plot with Monte Carlo Uncertainty')

hold off
else
%%%min-max
figure(12)
clf
hold on
res_upper_mm_sorted = sort(res_max);
res_lower_mm_sorted = sort(res_min);
fill([q_theory fliplr(q_theory)], ...
     [res_lower_mm_sorted fliplr(res_upper_mm_sorted)], ...
     [0.8 0.8 0.8], ...
     'EdgeColor','none', ...
     'FaceAlpha',0.5);
plot(q_theory, res_sorted, 'bx', 'MarkerSize',4, 'MarkerFaceColor','k')
plot(q_theory, q_theory, 'r-', 'LineWidth',1.5)

set(gca,'FontSize',13)
xlabel('Standard Normal Quantiles')
ylabel('Quality of Input Sample')
%title('QQ Plot with Monte Carlo Uncertainty (min-max)')
title('')
hold off

end

%outfile=['/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/paperfigure/figure_nonlinearmodel/Fig3_QQplot_MCtest.png'];print([outfile],'-dpng','-r300');
%outfile=['/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/paperfigure/figure_nonlinearmodel/Fig3_QQplot_b6constant.png'];print([outfile],'-dpng','-r300');


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
xlabel('Lag (months)')
set(gca,'fontsize',13)
title('')    
%outfile=['/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/paperfigure/figure_nonlinearmodel/Fig3_ACF.png'];print([outfile],'-dpng','-r300');

%%
figure(12)
autocorr(obs_webdig,'NumLags',24)
ylabel('ACF')
xlabel('Lag')
set(gca,'fontsize',13)
xlabel('Month')
%outfile=['/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/paperfigure/figure_nonlinearmodel/Fig3_ACF_observation.png'];print([outfile],'-dpng','-r300');

%% =========================
%  Ljung–Box Q-test (custmoize range)
%  =========================
n = numel(res);

% sample autocorrelations up to lag 12
maxLag = 24;
acf = autocorr(res,'NumLags',maxLag);   % returns lag 0..12
rho = acf(2:end);                       % lags 1..12

% choose lag range to test
lags = [4:1:16];

% Ljung–Box-type Q over selected lags only
Q = n*(n+2) * sum( (rho(lags).^2) ./ (n - lags)' );
df = numel(lags);
p = 1 - chi2cdf(Q, df);

fprintf('Custom LBQ over lags %d–%d: Q=%.3f, df=%d, p=%.4g\n', ...
        lags(1), lags(end), Q, df, p);



%%
%%%%
 [h,pValue,stat,cValue] = archtest(res), %%%%bp test show heteroscedasticity, residual change influenced by dependent variable changes

%res_index=(obs_webdig>5);
%figure(12)
%clf
%plot(yp1(res_index)',res(res_index),'r.');hold on
%plot(yp1(~res_index)',res(~res_index),'b.');hold on
%xlabel('Fitted values')
%ylabel('Residuals')
%yline(0,'k--')

%%%Spread–Location plot
res_std=res./std(res);
spread=sqrt(abs(res_std));
figure(23)
clf
scatter(yp1, spread, 25, 'filled');hold on
p = polyfit(yp1, spread, 1);
plot(yp1, polyval(p,yp1),'LineWidth',2)
xlabel('Fitted values')
ylabel('Sqrt(|Standardized residuals|)')
title('Scale-Location Plot')
grid on
set(gca,'fontsize',13)

outfile=['/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/paperfigure/figure_nonlinearmodel/Fig3_Scale-Location.png'];print([outfile],'-dpng','-r300');


%%
%%%%whether this is normal distribution
qqplot(res)
set(gca,'fontsize',13)
title('B')
%outfile=['/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/paperfigure/figure_nonlinearmodel/Fig3_QQplot.png'];print([outfile],'-dpng','-r300');
%outfile=['/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/paperfigure/figure_nonlinearmodel/Fig3_QQplot_b6constant.png'];print([outfile],'-dpng','-r300');
sk=skewness(res) %more large positive residuals than normal distribution
ku=kurtosis(res) %underpredicts extreme event, residuals become more postiive during extreme event.

