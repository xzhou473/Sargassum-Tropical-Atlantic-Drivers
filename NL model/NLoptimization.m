%% ===================== MAIN SCRIPT =====================

% --- Basic setup ---
lagmonth = 20;
b_dust   = b_dust;   
a_dust   = a_dust;  
    MLT_seasonal1(1:b_dust-a_dust+1) = NaN;
    for k = 1:(length(MLT_seasonal1)/3)-1
        idx = 3*k:(3*k)+2;
        MLT_seasonal1(idx) = mean(MLT(idx));
    end
    MLT_seasonal1(1:2) = mean(MLT(1:2));
    MLT_seasonal1(end) = mean(MLT(end));
    MLT_seasonal1 = MLT_seasonal1';

    MLT_avg = circshift(MLT_seasonal1, 0);
    MLT_avgnorm=(MLT_avg-min(MLT_avg))./(max(MLT_avg)-min(MLT_avg));

%%

% obs_webdig, MLT, MLT_shifted, DUST, SST_nonlinear,
% MLT_avgnorm2, conc_current_year, observation_dates
% should already be in workspace.
for run=1:50
% Initial weight guess: normalized random
try
weight0 = rand(lagmonth,1);
weight0 = weight0 ./ sum(weight0);

% Objective function: NEGATIVE correlation (since fminsearch MINIMIZES)
objfun = @(w)( -compute_R_for_weight(w, lagmonth, ...
                                     obs_webdig, ...
                                     MLT, MLT_shifted, ...
                                     DUST, SST_nonlinear, ...
                                     MLT_avgnorm, ...
                                     conc_current_year, ...
                                     observation_dates, ...
                                     a_dust, b_dust) );

% Optimization options
opts = optimset('MaxIter',800, ...
                'Display','iter', ...
                'TolX',1e-5, ...
                'TolFun',1e-5);

% Run optimization
[weight_opt, fval_opt] = fminsearch(objfun, weight0, opts);



% Best correlation
%%%%%

R_best = compute_R_for_weight(weight_opt, lagmonth, ...
                                  obs_webdig, ...
                                  MLT, MLT_shifted, ...
                                  DUST, SST_nonlinear, ...
                                  MLT_avgnorm, ...
                                  conc_current_year, ...
                                  observation_dates, ...
                                  a_dust, b_dust)

    weight_final = max(weight_opt(:), 0);             
    %if sum(weight) == 0
    %    weight = ones(size(weight));
    %end
    weight_final(end-2:end)=0;
    weight_final2 = weight_final ./ sum(weight_final);
weight_final_all(run,:)=weight_final;
Rall(run)=R_best;
catch ME
     fprintf('Run %d FAILED: %s\n', run, ME.message);
     weight_final_all(run,:) = weight_final;
     Rall(run) = R_best;   
     continue
end
end
%%    
for i=1:50
     weight_final2_all(i,:)=weight_final_all(i,:) ./ sum(weight_final_all(i,:));
end

figure(13)
clf
w_mean = mean(weight_final2_all, 1, 'omitnan');
w_min  = prctile(weight_final2_all, 25, 1);
w_max  = prctile(weight_final2_all, 75, 1);
lower = w_mean - w_min;   % distance from mean to min
upper = w_max - w_mean;   % distance from mean to max

[R_test,yp1] = compute_R_for_weight2(w_mean, lagmonth, ...
                                  obs_webdig, ...
                                  MLT, MLT_shifted, ...
                                  DUST, SST_nonlinear, ...
                                  MLT_avgnorm, ...
                                  conc_current_year, ...
                                  observation_dates, ...
                                  a_dust, b_dust);

errorbar([-19:1:0], w_mean, lower,upper,'-o','Color',[0.1 0.1 0.1],'LineWidth',1.5,'MarkerSize',4,'MarkerFaceColor',[0.1 0.1 0.1]); 
set(gca,'fontsize',13)
set(gca,'Tickdir','out')
xlabel('Lag Month');
ylabel('Weight');
ylim([0 0.5])
text(-8, 0.4, ...
     ['R = ' num2str(R_test, '%.2f') ' ± ' num2str(std(Rall), '%.2f')], ...
     'fontsize', 13);

%outfile=['/media/xzhou473/Seagate Backup Plus Drive/Sargassum/Sarg_for_Xing/Sarg_for_Xing/paperfigure/figure_nonlinearmodel/optimization/weight1.png'];print([outfile],'-dpng','-r300');
%%
weight_test=zeros(20,1);
weight_test(8:13)=1;

[R_test2,yp1] = compute_R_for_weight2(weight_test, lagmonth, ...
                                  obs_webdig, ...
                                  MLT, MLT_shifted, ...
                                  DUST, SST_nonlinear, ...
                                  MLT_avgnorm, ...
                                  conc_current_year, ...
                                  observation_dates, ...
                                  a_dust, b_dust)

%%
function R = compute_R_for_weight(weight, lagmonth, ...
                                  obs_webdig, ...
                                  MLT, MLT_shifted, ...
                                  DUST, SST_nonlinear, ...
                                  MLT_avgnorm2, ...
                                  conc_current_year, ...
                                  observation_dates, ...
                                  a_dust, b_dust)

    % ---------- 1. Ensure weight is "nice" ----------
    % Non-negative, normalized (optional but usually good):
    weight = max(weight(:), 0);               % enforce >= 0
    %if sum(weight) == 0
    %    weight = ones(size(weight));
    %end
    weight(end-2:end)=0; %%remove the t-2 to t, due to high correlation
    weight = weight ./ sum(weight);
    

    % ---------- 2. Build x from lag window and weight ----------
    N = b_dust - a_dust + 1;
    x0= NaN(N,1);
    x = NaN(N,1);

    %%%%x0 is just a movmean 3
    x0=movmean(obs_webdig,3);

    %for k = 1:(length(x0)/3)-1
    %observation_dates(3*k:(3*k)+2);
    %x0(3*k:(3*k)+2) = mean(obs_webdig(3*k:(3*k)+2));
    %end
    %x0(1:2) = mean(obs_webdig(1:2));
    %x0(end) = mean(obs_webdig(end));
 

   

    % Use your formula:
    for k = lagmonth:N
        idx = (k-lagmonth+1):k;
        x(k) = sum( (1 - MLT_avgnorm2(idx).^3) .* obs_webdig(idx)' .* weight' );

        %x(k) = sum( (1 - MLT_avgnorm2(idx).^3) .* x0(idx)' .* weight');

        
  
    end

    % Handle first lagmonth points
    idx0 = 1:lagmonth;
    x(idx0) = sum( (1 - MLT_avgnorm2(idx0).^3) .* obs_webdig(idx0) .* weight(1:numel(idx0)) );
    %x(idx0) = sum( (1 - MLT_avgnorm2(idx0).^3) .* x0(idx0) .* weight(1:numel(idx0)) );

    % ---------- 3. Only toty = 0 is used ----------
    toty = 0;
    CONC_LAG = toty;

    conc_previous = circshift(x, CONC_LAG);
    if CONC_LAG > 0
        conc_previous(1:CONC_LAG) = 0;
    else
        conc_previous(end+CONC_LAG+1:end) = 0;
    end

    % Extend conc_previous to match MLT_shifted length if needed:
    if numel(conc_previous) < numel(MLT_shifted)
        conc_previous = [conc_previous; ...
            zeros(numel(MLT_shifted) - numel(conc_previous), 1)];
    elseif numel(conc_previous) > numel(MLT_shifted)
        conc_previous = conc_previous(1:numel(MLT_shifted));
    end

    % ---------- 4. Compute seasonal MLT and normalized variables ----------
   

    % remaining_conc (your original code used a weird loop; here simplified)
    remaining_conc = conc_current_year;   % if this is already yearly conc
    if size(remaining_conc,2) == 144
        remaining_conc = remaining_conc';
    end

    % Normalizations (kept consistent with your code)
    %MLT_avgnorm1 = (MLT_avg - min(MLT_avg)) ./ (max(MLT_avg) - min(MLT_avg));
    % You had this line overwriting MLT_avgnorm1 with zeros; I assume that
    % was a test and not what you actually want, so I REMOVE that:
    MLT_avgnorm1 = zeros(144,1);  % <-- always 0

    MLT_avgnorm0 = (MLT_shifted - min(MLT_shifted)) ./ ...
                   (max(MLT_shifted) - min(MLT_shifted));

    DUST_avgnorm = (DUST - min(DUST)) ./ (max(DUST) - min(DUST));
    SST_avgnorm  = (SST_nonlinear - min(SST_nonlinear)) ./ ...
                   (max(SST_nonlinear) - min(SST_nonlinear));

    % ---------- 5. Build table for fitnlm ----------
    Training_ini = 1;
    Training_end = length(observation_dates);

    tbl = table(...
        MLT_shifted(Training_ini:Training_end), ...
        conc_previous(Training_ini:Training_end), ...
        DUST(Training_ini:Training_end), ...
        SST_nonlinear(Training_ini:Training_end), ...
        MLT_avgnorm1(Training_ini:Training_end), ...
        remaining_conc(Training_ini:Training_end), ...
        conc_current_year(Training_ini:Training_end), ...
        'VariableNames', ...
        {'MLT_shifted', 'conc_prev', 'DUST', 'SST', 'MLT_norm', 'rem_conc', 'conc_curr'} ...
        );

    % ---------- 6. Nonlinear model ----------
    modelfun = @(b,x) ( ...
        (x(:,1) >= 0) .* ( ...
            (b(1) .* (x(:,1)).^abs(b(5)) + ...
             b(2) .* (1 - x(:,5).^3) .* x(:,2) + ...
             b(3) .* x(:,3) + ...
             b(4) .* x(:,4)) .* ...
             sign(heaviside(b(1) .* (x(:,1)).^abs(b(5)) + ...
                            b(2) .* (1 - x(:,5).^3) .* x(:,2) + ...
                            b(3) .* x(:,3) + ...
                            b(4) .* x(:,4))) ...
        ) + ...
        (x(:,1) < 0) .* ( ...
            (-b(1) .* abs(x(:,1)).^abs(b(5)) + ...
              b(2) .* (1 - x(:,5).^3) .* x(:,2) + ...
              b(3) .* x(:,3) + ...
              b(4) .* x(:,4)) .* ...
              sign(heaviside(-b(1) .* abs(x(:,1)).^abs(b(5)) + ...
                             b(2) .* (1 - x(:,5).^3) .* x(:,2) + ...
                             b(3) .* x(:,3) + ...
                             b(4) .* x(:,4))) ...
        ) ...
    );

    beta0 = 5.0 .* [0.1 0.1 0.1 0.1 0.1];

    mdl_nl = fitnlm(tbl, modelfun, beta0);

    Bi = table2array(mdl_nl.Coefficients);
    b1 = Bi(1,1);
    b2 = Bi(2,1);
    b3 = Bi(3,1);
    b4 = Bi(4,1);  % you were missing this in your code
    b5 = Bi(5,1);

    % ---------- 7. Compute predictions yp1 ----------
    yp1 = zeros(length(MLT_shifted),1);
    for t = 1:length(MLT_shifted)
        if MLT_shifted(t) >= 0
            val =  b1.*(MLT_shifted(t)).^abs(b5) + ...
                   b2.*(1 - MLT_avgnorm1(t).^3).*conc_previous(t) + ...
                   b3.*DUST(t) + ...
                   b4.*SST_nonlinear(t);
            yp1(t) = val .* sign(heaviside(val));
        else
            val = -b1.*abs(MLT_shifted(t)).^abs(b5) + ...
                   b2.*(1 - MLT_avgnorm1(t).^3).*conc_previous(t) + ...
                   b3.*DUST(t) + ...
                   b4.*SST_nonlinear(t);
            yp1(t) = val .* sign(heaviside(val));
        end
    end

    % ---------- 8. Correlation with obs_webdig ----------
    % Make sure sizes align
    L = min(length(obs_webdig), length(yp1));
    R = corr(obs_webdig(1:L), yp1(1:L));

end

%%
function [R, yp1] = compute_R_for_weight2(weight, lagmonth, ...
                                  obs_webdig, ...
                                  MLT, MLT_shifted, ...
                                  DUST, SST_nonlinear, ...
                                  MLT_avgnorm2, ...
                                  conc_current_year, ...
                                  observation_dates, ...
                                  a_dust, b_dust)

    %%%%same function, but output a extra yp1
    % ---------- 1. Ensure weight is "nice" ----------
    % Non-negative, normalized (optional but usually good):
    weight = max(weight(:), 0);               % enforce >= 0
    %if sum(weight) == 0
    %    weight = ones(size(weight));
    %end
    weight(end-2:end)=0; %%remove the t-2 to t, due to high correlation
    weight = weight ./ sum(weight);
    

    % ---------- 2. Build x from lag window and weight ----------
    N = b_dust - a_dust + 1;
    x0= NaN(N,1);
    x = NaN(N,1);

    %%%%x0 is just a movmean 3
    x0=movmean(obs_webdig,3);

    %for k = 1:(length(x0)/3)-1
    %observation_dates(3*k:(3*k)+2);
    %x0(3*k:(3*k)+2) = mean(obs_webdig(3*k:(3*k)+2));
    %end
    %x0(1:2) = mean(obs_webdig(1:2));
    %x0(end) = mean(obs_webdig(end));
 

   

    % Use your formula:
    for k = lagmonth:N
        idx = (k-lagmonth+1):k;
        x(k) = sum( (1 - MLT_avgnorm2(idx).^3) .* obs_webdig(idx)' .* weight' );

        %x(k) = sum( (1 - MLT_avgnorm2(idx).^3) .* x0(idx)' .* weight');

        
  
    end

    % Handle first lagmonth points
    idx0 = 1:lagmonth;
    x(idx0) = sum( (1 - MLT_avgnorm2(idx0).^3) .* obs_webdig(idx0) .* weight(1:numel(idx0)) );
    %x(idx0) = sum( (1 - MLT_avgnorm2(idx0).^3) .* x0(idx0) .* weight(1:numel(idx0)) );

    % ---------- 3. Only toty = 0 is used ----------
    toty = 0;
    CONC_LAG = toty;

    conc_previous = circshift(x, CONC_LAG);
    if CONC_LAG > 0
        conc_previous(1:CONC_LAG) = 0;
    else
        conc_previous(end+CONC_LAG+1:end) = 0;
    end

    % Extend conc_previous to match MLT_shifted length if needed:
    if numel(conc_previous) < numel(MLT_shifted)
        conc_previous = [conc_previous; ...
            zeros(numel(MLT_shifted) - numel(conc_previous), 1)];
    elseif numel(conc_previous) > numel(MLT_shifted)
        conc_previous = conc_previous(1:numel(MLT_shifted));
    end

    % ---------- 4. Compute seasonal MLT and normalized variables ----------
   

    % remaining_conc (your original code used a weird loop; here simplified)
    remaining_conc = conc_current_year;   % if this is already yearly conc
    if size(remaining_conc,2) == 144
        remaining_conc = remaining_conc';
    end

    % Normalizations (kept consistent with your code)
    %MLT_avgnorm1 = (MLT_avg - min(MLT_avg)) ./ (max(MLT_avg) - min(MLT_avg));
    % You had this line overwriting MLT_avgnorm1 with zeros; I assume that
    % was a test and not what you actually want, so I REMOVE that:
    MLT_avgnorm1 = zeros(144,1);  % <-- always 0

    MLT_avgnorm0 = (MLT_shifted - min(MLT_shifted)) ./ ...
                   (max(MLT_shifted) - min(MLT_shifted));

    DUST_avgnorm = (DUST - min(DUST)) ./ (max(DUST) - min(DUST));
    SST_avgnorm  = (SST_nonlinear - min(SST_nonlinear)) ./ ...
                   (max(SST_nonlinear) - min(SST_nonlinear));

    % ---------- 5. Build table for fitnlm ----------
    Training_ini = 1;
    Training_end = length(observation_dates);

    tbl = table(...
        MLT_shifted(Training_ini:Training_end), ...
        conc_previous(Training_ini:Training_end), ...
        DUST(Training_ini:Training_end), ...
        SST_nonlinear(Training_ini:Training_end), ...
        MLT_avgnorm1(Training_ini:Training_end), ...
        remaining_conc(Training_ini:Training_end), ...
        conc_current_year(Training_ini:Training_end), ...
        'VariableNames', ...
        {'MLT_shifted', 'conc_prev', 'DUST', 'SST', 'MLT_norm', 'rem_conc', 'conc_curr'} ...
        );

    % ---------- 6. Nonlinear model ----------
    modelfun = @(b,x) ( ...
        (x(:,1) >= 0) .* ( ...
            (b(1) .* (x(:,1)).^abs(b(5)) + ...
             b(2) .* (1 - x(:,5).^3) .* x(:,2) + ...
             b(3) .* x(:,3) + ...
             b(4) .* x(:,4)) .* ...
             sign(heaviside(b(1) .* (x(:,1)).^abs(b(5)) + ...
                            b(2) .* (1 - x(:,5).^3) .* x(:,2) + ...
                            b(3) .* x(:,3) + ...
                            b(4) .* x(:,4))) ...
        ) + ...
        (x(:,1) < 0) .* ( ...
            (-b(1) .* abs(x(:,1)).^abs(b(5)) + ...
              b(2) .* (1 - x(:,5).^3) .* x(:,2) + ...
              b(3) .* x(:,3) + ...
              b(4) .* x(:,4)) .* ...
              sign(heaviside(-b(1) .* abs(x(:,1)).^abs(b(5)) + ...
                             b(2) .* (1 - x(:,5).^3) .* x(:,2) + ...
                             b(3) .* x(:,3) + ...
                             b(4) .* x(:,4))) ...
        ) ...
    );

    beta0 = 5.0 .* [0.1 0.1 0.1 0.1 0.1];

    mdl_nl = fitnlm(tbl, modelfun, beta0);

    Bi = table2array(mdl_nl.Coefficients);
    b1 = Bi(1,1);
    b2 = Bi(2,1);
    b3 = Bi(3,1);
    b4 = Bi(4,1);  % you were missing this in your code
    b5 = Bi(5,1);

    % ---------- 7. Compute predictions yp1 ----------
    yp1 = zeros(length(MLT_shifted),1);
    for t = 1:length(MLT_shifted)
        if MLT_shifted(t) >= 0
            val =  b1.*(MLT_shifted(t)).^abs(b5) + ...
                   b2.*(1 - MLT_avgnorm1(t).^3).*conc_previous(t) + ...
                   b3.*DUST(t) + ...
                   b4.*SST_nonlinear(t);
            yp1(t) = val .* sign(heaviside(val));
        else
            val = -b1.*abs(MLT_shifted(t)).^abs(b5) + ...
                   b2.*(1 - MLT_avgnorm1(t).^3).*conc_previous(t) + ...
                   b3.*DUST(t) + ...
                   b4.*SST_nonlinear(t);
            yp1(t) = val .* sign(heaviside(val));
        end
    end

    % ---------- 8. Correlation with obs_webdig ----------
    % Make sure sizes align
    L = min(length(obs_webdig), length(yp1));
    R = corr(obs_webdig(1:L), yp1(1:L));

end