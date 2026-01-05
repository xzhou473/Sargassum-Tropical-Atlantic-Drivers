%%
clear all
%%%%% deal with river discharge data
file = '/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/riverdata/Obidos_discharge.csv';

opts = detectImportOptions(file);
opts.Delimiter = ',';

% First tell MATLAB that column 1 is datetime
opts.VariableTypes{1} = 'datetime';

% Now set the input format for that datetime column
opts = setvaropts(opts, 1, 'InputFormat', 'yyyy/MM/dd HH:mm');

data_A = readtable(file, opts);
file = '/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/riverdata/Ciudad_Bolivar_discharge.csv';

opts = detectImportOptions(file);
opts.Delimiter = ',';

% First tell MATLAB that column 1 is datetime
opts.VariableTypes{1} = 'datetime';

% Now set the input format for that datetime column
opts = setvaropts(opts, 1, 'InputFormat', 'yyyy/MM/dd HH:mm');

data_B = readtable(file, opts);


data_A=data_A([9372:end],:); %%1993-2023
data_B=data_B([24314:end],:); %%1993-2018


%%
data_A1 = process_discharge_data(data_A);
data_B1 = process_discharge_data(data_B);

[data_A1_month(:,1), data_A1_month(:,2),data_A1_month(:,3)] = daily2monthly_datenum(data_A1(:,1),data_A1(:,2))
[data_B1_month(:,1), data_B1_month(:,2),data_B1_month(:,3)] = daily2monthly_datenum(data_B1(:,1),data_B1(:,2))





%%
%%%%%deal with NO3 data
file = '/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/riverdata/Ciudad_Bolivar_NO3.csv';

opts = detectImportOptions(file);
opts.Delimiter = ',';

% First tell MATLAB that column 1 is datetime
opts.VariableTypes{1} = 'datetime';

% Now set the input format for that datetime column
opts = setvaropts(opts, 1, 'InputFormat', 'yyyy/MM/dd');

NO3_B = readtable(file, opts);

file = '/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/riverdata/Obidos_NO3.csv';

opts = detectImportOptions(file);
opts.Delimiter = ',';

% First tell MATLAB that column 1 is datetime
opts.VariableTypes{1} = 'datetime';

% Now set the input format for that datetime column
opts = setvaropts(opts, 1, 'InputFormat', 'yyyy/MM/dd');

NO3_A = readtable(file, opts);

%%
NO3_A1 = process_NO3_monthly(NO3_A);
NO3_B1 = process_NO3_monthly(NO3_B);


%%
MWNO3 = 62.01/1000;   % µM → mg/L → g/m^3
xtickyears = 1990:5:2025;
xt = datenum(xtickyears,1,1);



%%%% ============================================================
% 1. DISCHARGE — dual y-axis
%%%% ============================================================
figure(10)
clf
scale2=80;
scale=1; 
set(gcf,'position', [474   157     12*scale2 6*scale2])
set(gcf,'paperposition', [0   0   12*scale  6*scale])
yyaxis left
plot(data_A1_month(1:end,1), data_A1_month(1:end,2), 'color',[26,26,26]./256, 'LineWidth', 2.0); hold on
ylabel('Amazon discharge (m^3/s)')
set(gca,'YColor',[26,26,26]./256)

yyaxis right
plot(data_B1_month(1:end,1), data_B1_month(1:end,2), 'color',[178,24,43]./256, 'LineWidth', 2.0); hold on
ylabel('Orinoco discharge (m^3/s)')
set(gca,'YColor',[178,24,43]./256)

set(gca,'FontSize',13)
set(gca,'XTick', xt, 'XTickLabel', string(xtickyears))
xlabel('Year')

hl1 = legend('Amazon river (Obidos)', ...
             'Orinoco river (Ciudad Bolivar)', ...
             'Location','North');
outfile=['/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/riverdata/riverinformation1.png']; print([outfile],'-dpng','-r300');
% optional tiny legend box at top
%set(hl1,'Position',[0.50 0.96 0.00001 0.00001]);
%%
%%%% ============================================================
% 2. NITRATE CONCENTRATION — dual y-axis
%%%% ============================================================
figure(11)
clf
scale2=80;
scale=1; 
set(gcf,'position', [474   157     12*scale2 6*scale2])
set(gcf,'paperposition', [0   0   12*scale  6*scale])
yyaxis left
plot(NO3_A1(1:end,1), NO3_A1(1:end,2)*MWNO3, '-o','color',[26,26,26]./256, 'LineWidth', 2.0); hold on
ylabel('Amazon NO_3 (g/m^3)')
set(gca,'YColor',[26,26,26]./256)

yyaxis right
plot(NO3_B1(1:end,1), NO3_B1(1:end,2)*MWNO3, '-o','color',[178,24,43]./256, 'LineWidth', 2.0); hold on
ylabel('Orinoco NO_3 (g/m^3)')
set(gca,'YColor',[178,24,43]./256)
set(gca,'FontSize',13)
set(gca,'XTick', xt, 'XTickLabel', string(xtickyears))
xlabel('Year')
hl1 = legend('Amazon river (Obidos)', ...
             'Orinoco river (Ciudad Bolivar)', ...
             'Location','Northwest');
% optional tiny legend box at top
%set(hl1,'Position',[0.50 0.96 0.00001 0.00001]);
outfile=['/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/riverdata/riverinformatio2.png']; print([outfile],'-dpng','-r300');
%%
figure(12)
clf
scale2=80;
scale=1; 
set(gcf,'position', [474   157     12*scale2 6*scale2])
set(gcf,'paperposition', [0   0   12*scale  6*scale])
%%%% ============================================================
% 3. NITRATE FLUX — dual y-axis
%%%% ============================================================
flux_A = NO3_A1(:,2).*data_A1_month(121:end,2).*1e-3.*MWNO3;        % kg/s
flux_B = NO3_B1(1:end-14,2).*data_B1_month(121:end,2).*1e-3.*MWNO3; % kg/s

yyaxis left
plot(NO3_A1(:,1), flux_A, 'r','color',[26,26,26]./256, 'LineWidth', 2.0); hold on
ylabel('Amazon NO_3 flux (kg/s)')
set(gca,'YColor',[26,26,26]./256)

yyaxis right
plot(NO3_B1(1:end-14,1), flux_B, 'color',[178,24,43]./256, 'LineWidth', 2.0); hold on
ylabel('Orinoco NO_3 flux (kg/s)')
set(gca,'YColor',[178,24,43]./256)

set(gca,'FontSize',13)
set(gca,'XTick', xt, 'XTickLabel', string(xtickyears))
xlabel('Year')
hl1 = legend('Amazon river (Obidos)', ...
             'Orinoco river (Ciudad Bolivar)', ...
             'Location','Northwest');

outfile=['/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/riverdata/riverinformation3.png']; print([outfile],'-dpng','-r300');

%%
%%%%%%read sargassum data here, and do a correlation
TAB = csvread('../utils/digitized_sarg_biomass/WebDigitized_Jan2011-Dec2022_from_Jouanno-et-al-2023.csv');

obs_webdig = TAB(:,2); 

% It's impossible to have negative concentration. Manual digitizing might
% introduce spurious negative values when concentrations are very small. 
% Hence, put to zero every negative point.

obs_webdig(obs_webdig<0) = 0;


% Create date vector for digitized observations:  Jan 2011 to Dec 2022
observation_dates = data_A1_month([217:1:217+143],1)


  [anom_A1, clim] = deseasonalize_monthly(data_A1_month(:,1), data_A1_month(:,2));
  [anom_B1, clim] = deseasonalize_monthly(data_B1_month(:,1), data_B1_month(:,2));
  [anom_all, clim] = deseasonalize_monthly(data_B1_month(:,1), data_B1_month(:,2)+data_A1_month(1:end-64,2));

corr(data_A1_month([217:1:217+143],2),obs_webdig)
corr(data_B1_month([217:1:end],2),obs_webdig(1:end-52))
corr(data_B1_month([217:1:end],2)+data_A1_month([217:1:217+143-52],2),obs_webdig(1:end-52))

corr(anom_A1([217:1:217+143]),obs_webdig)
corr(anom_B1([217:1:end]),obs_webdig(1:end-52))
corr(anom_all([217:1:end]),obs_webdig(1:end-52))

corr(flux_A([97:end-12]),obs_webdig, 'Rows','complete')
corr(flux_B([97:end]),obs_webdig(1:end-52), 'Rows','complete')
corr(flux_B([97:end])+flux_A([97:end-12-52]),obs_webdig(1:end-52), 'Rows','complete')

flux_A_filled = fillmissing(flux_A, 'linear')

%%
function data_A1 = process_discharge_data(data_A)
% PROCESS_DISCHARGE_DATA
% Converts irregular m3/s discharge data to continuous daily m3/day time series.
%
% INPUT:
%   data_A: table
%       column 1 = datetime
%       column 3 = discharge (m3/s)
%
% OUTPUT:
%   data_A1: Nx2 matrix
%       column 1 = datenum (daily)
%       column 2 = discharge (m3/day), with interpolation where appropriate

    % -----------------------------
    % 1. Extract time and discharge
    % -----------------------------
    time = data_A{:,1};
    Q0 = data_A{:,3};  % m3/s
    for i=1:length(Q0)
    Q1=Q0{i};;
    Q(i,:)=str2double(Q1);
    end

    % Round to day
    day_only = dateshift(time, "start", "day");

    % -------------------------------------------------
    % 2. Average multiple observations within each day
    % -------------------------------------------------
    [unique_days, ~, idx] = unique(day_only);
    daily_mean_Q = accumarray(idx, Q, [], @mean);   % m3/s

    % -------------------------------------------------
    % 3. Convert discharge from m3/s → m3/day
    % -------------------------------------------------
    daily_m3day = daily_mean_Q; %* 86400;   % 86400 sec/day

    % -------------------------------------------------
    % 4. Build continuous daily time vector
    % -------------------------------------------------
    full_days = (unique_days(1):unique_days(end))';
    N = numel(full_days);
    full_days_num = datenum(full_days);

    full_Q = NaN(N,1);

    % Fill known values
    [~, loc] = ismember(unique_days, full_days);
    full_Q(loc) = daily_m3day;

    % -------------------------------------------------
    % 5. Identify missing segments
    % -------------------------------------------------
    isn_missing = isnan(full_Q);

    % Find contiguous NaN blocks
    d = diff([false; isn_missing; false]);
    start_missing = find(d == 1);
    end_missing   = find(d == -1) - 1;

    % -------------------------------------------------
    % 6. Interpolate short gaps (≤ 10 days), but ONLY 
    %    if both boundaries exist and are valid
    % -------------------------------------------------
    maxGap = 3;   % threshold for interpolation

    for i = 1:length(start_missing)
        s = start_missing(i);
        e = end_missing(i);
        gap_length = e - s + 1;

        if gap_length <= maxGap
            left_idx  = s - 1;
            right_idx = e + 1;

            % Check boundary validity
            if left_idx >= 1 && right_idx <= N && ...
               ~isnan(full_Q(left_idx)) && ~isnan(full_Q(right_idx))

                full_Q(s:e) = interp1( ...
                    [full_days_num(left_idx); full_days_num(right_idx)], ...
                    [full_Q(left_idx);        full_Q(right_idx)], ...
                    full_days_num(s:e) ...
                );
            else
                % Missing boundary → keep NaN
            end
        else
            % Gap > 10 days → keep NaN
        end
    end

    % -------------------------------------------------
    % 7. Final output: datenum + discharge
    % -------------------------------------------------
    data_A1 = [full_days_num, full_Q];

end

%%
function NO3_monthly = process_NO3_monthly(data_NO3)
% PROCESS_NO3_MONTHLY
% Input: table with:
%   column 1 = datetime
%   column 2 = NO3 concentration (may be char/cell/string)
%
% Output: NO3_monthly (Nx3 matrix)
%   Column 1: datenum of month start (YYYY-MM-01)
%   Column 2: monthly mean NO3 (mg/L)
%   Column 3: number of samples in that month

    % -------------------------------
    % 1. Extract datetime
    % -------------------------------
    time = data_NO3{:,1};

    % -------------------------------
    % 2. Extract NO3 column and convert to numeric
    % -------------------------------
    NO3 = data_NO3{:,2}; 
 % -------------------------------
    % 3. Convert to Year-Month only
    % -------------------------------
    ym = dateshift(time, "start", "month");

    % -------------------------------
    % 4. Group by month
    % -------------------------------
    [unique_months, ~, idx] = unique(ym);

    monthly_mean = accumarray(idx, NO3, [], @mean);
    monthly_count = accumarray(idx, NO3, [], @(x) sum(~isnan(x)));

    % -------------------------------
    % 5. Build continuous month vector
    % -------------------------------
    full_months = (unique_months(1):calmonths(1):unique_months(end))';
    N = numel(full_months);

    full_mean = NaN(N,1);
    full_count = zeros(N,1);

    [~, loc] = ismember(unique_months, full_months);
    full_mean(loc) = monthly_mean;
    full_count(loc) = monthly_count;

    % -------------------------------
    % 6. Convert month -> datenum
    % -------------------------------
    month_num = datenum(full_months);

    % -------------------------------
    % 7. Final output
    % -------------------------------
    NO3_monthly = [month_num, full_mean, full_count];
end
%%
function [monthly_dn, monthly_Q, monthly_N] = daily2monthly_datenum(date_dn, Q)
% DAILY2MONTHLY_DATENUM
% Convert daily discharge (Q) into monthly averages
% date_dn: vector of datenum dates
% Q:       vector of daily discharge
%
% Outputs:
% monthly_dn : datenum of first day of each month
% monthly_Q  : monthly averaged discharge
% monthly_N  : number of daily points in each month

    % Convert datenum → datetime
    t = datetime(date_dn, 'ConvertFrom', 'datenum');

    % Extract year and month
    [y, m] = ymd(t);

    % Group by year-month
    G = findgroups(y, m);

    % Monthly mean
    monthly_Q = splitapply(@mean, Q, G);

    % Number of daily values per month
    monthly_N = splitapply(@numel, Q, G);

    % Monthly datenum (first day of month)
    monthly_time = datetime(y, m, 1);
    monthly_time = monthly_time(splitapply(@(x) x(1), (1:length(t))', G));

    % Convert back to datenum
    monthly_dn = datenum(monthly_time);
end

%%
function [anom, clim] = deseasonalize_monthly(time_dn, data)
% DESEASONALIZE_MONTHLY
% Remove the monthly climatology from a monthly time series.
%
% INPUTS:
%   time_dn : datenum vector (monthly timestamps)
%   data    : monthly data (same length)
%
% OUTPUTS:
%   anom : deseasonalized anomalies
%   clim : monthly climatology (12x1)
%
% USAGE:
%   [anom, clim] = deseasonalize_monthly(t, Q);
%

    % ensure column
    time_dn = time_dn(:);
    data = data(:);

    % convert to datetime
    t = time_dn;%datetime(time_dn, 'ConvertFrom', 'datenum');

    % extract month number (1–12)
    m = month(t);

    % -------------------------------------------------------------
    % 1. Compute monthly climatology (mean for each month 1–12)
    % -------------------------------------------------------------
    clim = nan(12,1);
    for mm = 1:12
        idx = (m == mm);
        clim(mm) = mean(data(idx), 'omitnan');
    end

    % -------------------------------------------------------------
    % 2. Generate the deseasonalized anomaly time series
    % -------------------------------------------------------------
    anom = data - clim(m);

end
