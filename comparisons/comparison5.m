%% Comparison of performance MDICC using DS clustering and MDICC using DS clustering and feature selection
% influence of feature selection on MDICC using DS clustering

close all
clear all

clc

%% Load ARI and NMI scores per dataset

% LIHC
ari_lihc_fsvt = readmatrix('C:\Users\floor\OneDrive\Documents\MDICC\MDICC-main\MDICC-main\lihc\ari_scores_lihc_fsvt.csv'); 
nmi_lihc_fsvt = readmatrix('C:\Users\floor\OneDrive\Documents\MDICC\MDICC-main\MDICC-main\lihc\nmi_scores_lihc_fsvt.csv'); 

ari_lihc_fsk = readmatrix('C:\Users\floor\OneDrive\Documents\MDICC\MDICC-main\MDICC-main\lihc\ari_scores_lihc_fsk.csv'); 
nmi_lihc_fsk = readmatrix('C:\Users\floor\OneDrive\Documents\MDICC\MDICC-main\MDICC-main\lihc\nmi_scores_lihc_fsk.csv'); 

max_scores_lihc = readmatrix('C:\Users\floor\OneDrive\Documents\MDICC\MDICC-main\MDICC-main\lihc\opt_prec_scores_lihc_ds.txt');
max_scores_lihc = max_scores_lihc(:, 1:2);

% KIRC
ari_kirc_fsvt = readmatrix('C:\Users\floor\OneDrive\Documents\MDICC\MDICC-main\MDICC-main\kirc\ari_scores_kirc_fsvt.csv'); 
nmi_kirc_fsvt = readmatrix('C:\Users\floor\OneDrive\Documents\MDICC\MDICC-main\MDICC-main\kirc\nmi_scores_kirc_fsvt.csv'); 

ari_kirc_fsk = readmatrix('C:\Users\floor\OneDrive\Documents\MDICC\MDICC-main\MDICC-main\kirc\ari_scores_kirc_fsk.csv'); 
nmi_kirc_fsk = readmatrix('C:\Users\floor\OneDrive\Documents\MDICC\MDICC-main\MDICC-main\kirc\nmi_scores_kirc_fsk.csv');

max_scores_kirc = readmatrix('C:\Users\floor\OneDrive\Documents\MDICC\MDICC-main\MDICC-main\kirc\opt_prec_scores_kirc_ds.txt');
max_scores_kirc = max_scores_kirc(:, 1:2);

%% Set precision

precision = 0.0001:0.0001:0.04;

%% Feature selection SelectKBest
%% Calculate the maximum ARI and NMI and determine the corresponding precision value

% LIHC
% find maximum of ARI and NMI and their indices
[max_ari_lihc_fsk, idx1_lihc_fsk] = max(ari_lihc_fsk);
[max_nmi_lihc_fsk, idx2_lihc_fsk] = max(nmi_lihc_fsk);

if idx1_lihc_fsk == idx2_lihc_fsk % check if the indices of max ari and nmi are the same
    optimal_precision_lihc_fsk = precision(idx1_lihc_fsk); % determine the optimal precision value
    disp("LIHC - Precision of max ARI and NMI score is " + optimal_precision_lihc_fsk)
else
    disp("Maximum values do not have the same index.")
    optimal_precision_lihc_fsk = precision(idx1_lihc_fsk);
end

% write optimal scores and corresponding precision to file
writematrix([max_ari_lihc_fsk, max_nmi_lihc_fsk, optimal_precision_lihc_fsk],'C:\Users\floor\OneDrive\Documents\MDICC\MDICC-main\MDICC-main\lihc\opt_prec_scores_lihc_ds_fsk.txt','Delimiter','space')

% KIRC
% find maximum of ARI and NMI and their indices
[max_ari_kirc_fsk, idx1_kirc_fsk] = max(ari_kirc_fsk);
[max_nmi_kirc_fsk, idx2_kirc_fsk] = max(nmi_kirc_fsk);

if idx1_kirc_fsk == idx2_kirc_fsk % check if the indices of max ari and nmi are the same
    optimal_precision_kirc_fsk = precision(idx1_kirc_fsk); % determine the optimal precision value
    disp("KIRC - Precision of max ARI and NMI score is " + optimal_precision_kirc_fsk)
else
    disp("Maximum values do not have the same index.")
    optimal_precision_kirc_fsk = precision(idx1_kirc_fsk);
end

% write optimal scores and corresponding precision to file
writematrix([max_ari_kirc_fsk, max_nmi_kirc_fsk, optimal_precision_kirc_fsk],'C:\Users\floor\OneDrive\Documents\MDICC\MDICC-main\MDICC-main\kirc\opt_prec_scores_kirc_ds_fsk.txt','Delimiter','space')

%% calculate average of optimal precision values

opt_precs = [optimal_precision_lihc_fsk, optimal_precision_kirc_fsk];

% compute mean and round result to 4 decimals
avg_precision_fsk = round(mean(opt_precs), 4);

disp("The average precision is: " + avg_precision_fsk)

% find index average precision, account for floating-point numbers
if isfloat(avg_precision_fsk)
    tol = eps * 2; % floating point number tolerance
    idx_avg_prec_fsk = find(abs(precision - avg_precision_fsk) <= tol);
else
    idx_avg_prec_fsk = find(precision == avg_precision_fsk);
end

%% Compute ARI and NMI score per data set for the average precision

% LIHC
ap_ari_lihc_fsk = ari_lihc_fsk(idx_avg_prec_fsk);
ap_nmi_lihc_fsk = nmi_lihc_fsk(idx_avg_prec_fsk);

disp("lihc: ")
disp([ap_ari_lihc_fsk ap_nmi_lihc_fsk])

% write ARI and NMI for average precision to file
fileID = fopen('C:\Users\floor\OneDrive\Documents\MDICC\MDICC-main\MDICC-main\lihc\avg_prec_lihc_fsk.txt', 'w');
formatSpec = 'average precision is %6.4f \nari score is %6.8f \nnmi score is %6.8f';
fprintf(fileID, formatSpec, avg_precision_fsk, ap_ari_lihc_fsk, ap_nmi_lihc_fsk);

% KIRC
ap_ari_kirc_fsk = ari_kirc_fsvt(idx_avg_prec_fsk);
ap_nmi_kirc_fsk = nmi_kirc_fsvt(idx_avg_prec_fsk);

disp("kirc: ")
disp([ap_ari_kirc_fsk, ap_nmi_kirc_fsk])

% write ARI and NMI for average precision to file
fileID = fopen('C:\Users\floor\OneDrive\Documents\MDICC\MDICC-main\MDICC-main\kirc\avg_prec_kirc_fsk.txt', 'w');
formatSpec = 'average precision is %6.4f \nari score is %6.8f \nnmi score is %6.8f';
fprintf(fileID, formatSpec, avg_precision_fsk, ap_ari_kirc_fsk, ap_nmi_kirc_fsk);

%% Feature selection VarianceThreshold
%% Calculate the maximum ARI and NMI and determine the corresponding precision value

% LIHC
% find maximum of ARI and NMI and their indices
[max_ari_lihc_fsvt, idx1_lihc_fsvt] = max(ari_lihc_fsvt);
[max_nmi_lihc_fsvt, idx2_lihc_fsvt] = max(nmi_lihc_fsvt);

if idx1_lihc_fsvt == idx2_lihc_fsvt % check if the indices of max ari and nmi are the same
    optimal_precision_lihc_fsvt = precision(idx1_lihc_fsvt); % determine the optimal precision value
    disp("LIHC - Precision of max ARI and NMI score is " + optimal_precision_lihc_fsvt)
else
    disp("Maximum values do not have the same index.")
    optimal_precision_lihc_fsvt = precision(idx1_lihc_fsvt);
end

% write optimal scores and corresponding precision to file
writematrix([max_ari_lihc_fsvt, max_nmi_lihc_fsvt, optimal_precision_lihc_fsvt],'C:\Users\floor\OneDrive\Documents\MDICC\MDICC-main\MDICC-main\lihc\opt_prec_scores_lihc_ds_fsvt.txt','Delimiter','space')

% KIRC
% find maximum of ARI and NMI and their indices
[max_ari_kirc_fsvt, idx1_kirc_fsvt] = max(ari_kirc_fsvt);
[max_nmi_kirc_fsvt, idx2_kirc_fsvt] = max(nmi_kirc_fsvt);

if idx1_kirc_fsvt == idx2_kirc_fsvt % check if the indices of max ari and nmi are the same
    optimal_precision_kirc_fsvt = precision(idx1_kirc_fsvt); % determine the optimal precision value
    disp("KIRC - Precision of max ARI and NMI score is " + optimal_precision_kirc_fsvt)
else
    disp("Maximum values do not have the same index.")
    optimal_precision_kirc_fsvt = precision(idx1_kirc_fsvt);
end

% write optimal scores and corresponding precision to file
writematrix([max_ari_kirc_fsvt, max_nmi_kirc_fsvt, optimal_precision_kirc_fsvt],'C:\Users\floor\OneDrive\Documents\MDICC\MDICC-main\MDICC-main\kirc\opt_prec_scores_kirc_ds_fsvt.txt','Delimiter','space')

%% Calculate average precision of the optimal precision value of the datasets

opt_precs = [optimal_precision_lihc_fsvt, optimal_precision_kirc_fsvt];

% compute mean and round average to 4 decimals
avg_precision_fsvt = round(mean(opt_precs), 4);

disp("The average precision is: " + avg_precision_fsvt)

% find index average precision, account for floating-point numbers
if isfloat(avg_precision_fsvt)
    tol = eps * 2; % floating point number tolerance
    idx_avg_prec_fsvt = find(abs(precision - avg_precision_fsvt) <= tol);
else
    idx_avg_prec_fsvt = find(precision == avg_precision_fsvt);
end

%% Compute ARI and NMI score per data set for the average precision

% LIHC
ap_ari_lihc_fsvt = ari_lihc_fsvt(idx_avg_prec_fsvt);
ap_nmi_lihc_fsvt = nmi_lihc_fsvt(idx_avg_prec_fsvt);

disp("lihc: ")
disp([ap_ari_lihc_fsvt ap_nmi_lihc_fsvt])

% write ARI and NMI for average precision to file
fileID = fopen('C:\Users\floor\OneDrive\Documents\MDICC\MDICC-main\MDICC-main\lihc\avg_prec_lihc_fsvt.txt', 'w');
formatSpec = 'average precision is %6.4f \nari score is %6.8f \nnmi score is %6.8f';
fprintf(fileID, formatSpec, avg_precision_fsvt, ap_ari_lihc_fsvt, ap_nmi_lihc_fsvt);

% KIRC
ap_ari_kirc_fsvt = ari_kirc_fsvt(idx_avg_prec_fsvt);
ap_nmi_kirc_fsvt = nmi_kirc_fsvt(idx_avg_prec_fsvt);

disp("kirc: ")
disp([ap_ari_kirc_fsvt, ap_nmi_kirc_fsvt])

% write ARI and NMI for average precision to file
fileID = fopen('C:\Users\floor\OneDrive\Documents\MDICC\MDICC-main\MDICC-main\kirc\avg_prec_kirc_fsvt.txt', 'w');
formatSpec = 'average precision is %6.4f \nari score is %6.8f \nnmi score is %6.8f';
fprintf(fileID, formatSpec, avg_precision_fsvt, ap_ari_kirc_fsvt, ap_nmi_kirc_fsvt);

%% Plot ARI and NMI scores against precision 
% plot includes average precision and optimal precision and 
% plot includes scores of MDICC with DS clustering as horizontal line

fig = figure;

% LIHC - VarianceThreshold
subplot(2,2,1);
hold on
plot(precision, ari_lihc_fsvt, 'b')
plot(precision, nmi_lihc_fsvt, 'r')

% add vertical lines (average and optimal precision) to plot
xl1 = xline(avg_precision_fsvt,'-','Average');
xl2 = xline(optimal_precision_lihc_fsvt,'-','Optimal');
xl1.LabelVerticalAlignment = 'bottom';
xl1.LabelHorizontalAlignment = 'left';
xl2.LabelVerticalAlignment = 'bottom';
xl2.LabelHorizontalAlignment = 'left';

% add horizontal lines (scores using DS) to plot
yl1 = yline(max_scores_lihc(1,1), 'b', 'ARI DS');
yl2 = yline(max_scores_lihc(1,2), 'r', 'NMI DS');
yl1.LabelVerticalAlignment = 'top';
yl1.LabelHorizontalAlignment = 'left';
yl2.LabelVerticalAlignment = 'top';
yl2.LabelHorizontalAlignment = 'left';

legend('ARI DS with FS', 'NMI DS with FS', 'Location', 'southwest')
title('LIHC')
ylabel({'VarianceThreshold   '}, 'Rotation' , 0, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right');
hold off

% set size of plot
x=10;
y=10;
width=1000;
height=1000;

set(gcf, 'position', [x,y,width,height])
ylim([0, 1])

% KIRC - VarianceThreshold
subplot(2,2,2)
hold on
plot(precision, ari_kirc_fsvt, 'b')
plot(precision, nmi_kirc_fsvt, 'r')

% add vertical lines (average and optimal precision) to plot
xl1 = xline(avg_precision_fsvt,'-','Average');
xl2 = xline(optimal_precision_kirc_fsvt,'-','Optimal');
xl1.LabelVerticalAlignment = 'bottom';
xl1.LabelHorizontalAlignment = 'left';
xl2.LabelVerticalAlignment = 'bottom';
xl2.LabelHorizontalAlignment = 'right';

% add horizontal lines (scores using DS) to plot
yl1 = yline(max_scores_kirc(1,1), 'b', 'ARI DS');
yl2 = yline(max_scores_kirc(1,2), 'r', 'NMI DS');
yl1.LabelVerticalAlignment = 'top';
yl1.LabelHorizontalAlignment = 'left';
yl2.LabelVerticalAlignment = 'bottom';
yl2.LabelHorizontalAlignment = 'left';

legend('ARI DS with FS', 'NMI DS with FS', 'Location', 'southwest')
title('KIRC')
hold off

set(gcf, 'position', [x,y,width,height])
ylim([0, 1])

% LIHC - SelectKBest
subplot(2,2,3);
hold on
plot(precision, ari_lihc_fsk, 'b')
plot(precision, nmi_lihc_fsk, 'r')

% add vertical lines (average and optimal precision) to plot
xl1 = xline(avg_precision_fsk,'-','Average');
xl2 = xline(optimal_precision_lihc_fsk,'-','Optimal');
xl1.LabelVerticalAlignment = 'bottom';
xl1.LabelHorizontalAlignment = 'left';
xl2.LabelVerticalAlignment = 'bottom';
xl2.LabelHorizontalAlignment = 'right';

% add horizontal lines (scores using DS) to plot
yl1 = yline(max_scores_lihc(1,1), 'b', 'ARI DS');
yl2 = yline(max_scores_lihc(1,2), 'r', 'NMI DS');
yl1.LabelVerticalAlignment = 'top';
yl1.LabelHorizontalAlignment = 'left';
yl2.LabelVerticalAlignment = 'bottom';
yl2.LabelHorizontalAlignment = 'left';

legend('ARI DS with FS', 'NMI DS with FS', 'Location', 'southwest')
ylabel({'SelectKBest   '}, 'Rotation' , 0, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right');
hold off

set(gcf, 'position', [x,y,width,height])
ylim([0, 1])

% KIRC - SelectKBest
subplot(2,2,4)
hold on
plot(precision, ari_kirc_fsk, 'b')
plot(precision, nmi_kirc_fsk, 'r')

% add vertical lines (average and optimal precision) to plot
xl1 = xline(avg_precision_fsk,'-','Average');
xl2 = xline(optimal_precision_kirc_fsk,'-','Optimal');
xl1.LabelVerticalAlignment = 'bottom';
xl1.LabelHorizontalAlignment = 'right';
xl2.LabelVerticalAlignment = 'bottom';
xl2.LabelHorizontalAlignment = 'left';

% add horizontal lines (scores using DS) to plot
yl1 = yline(max_scores_kirc(1,1), 'b', 'ARI DS');
yl2 = yline(max_scores_kirc(1,2), 'r', 'NMI DS');
yl1.LabelVerticalAlignment = 'top';
yl1.LabelHorizontalAlignment = 'left';
yl2.LabelVerticalAlignment = 'bottom';
yl2.LabelHorizontalAlignment = 'left';

legend('ARI DS with FS', 'NMI  DS with FS', 'Location', 'southwest')
hold off

set(gcf, 'position', [x,y,width,height])
ylim([0, 1])

% set x label, y label and title for plot
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'evaluation score');
xlabel(han,'precision');
sgtitle('Comparison of MDICC with DS to MDICC with DS using feature selection')

% save figure as pdf file
exportgraphics(fig,'C:\Users\floor\OneDrive\Documents\MDICC\MDICC-main\MDICC-main\comparisons\comparison5.pdf','ContentType','vector')