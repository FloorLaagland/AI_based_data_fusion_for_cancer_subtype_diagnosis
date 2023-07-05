%% Comparison of performance original MDICC to MDICC using DS clustering method
% clustering methods that are used are Kmeans++ and DS clustering

close all
clear all

clc

%% Load ARI and NMI scores per dataset

% LIHC
ari_lihc = readmatrix('C:\Users\floor\OneDrive\Documents\MDICC\MDICC-main\MDICC-main\lihc\ari_scores_lihc.csv'); 
nmi_lihc = readmatrix('C:\Users\floor\OneDrive\Documents\MDICC\MDICC-main\MDICC-main\lihc\nmi_scores_lihc.csv'); 

% KIRC
ari_kirc = readmatrix('C:\Users\floor\OneDrive\Documents\MDICC\MDICC-main\MDICC-main\kirc\ari_scores_kirc.csv'); 
nmi_kirc = readmatrix('C:\Users\floor\OneDrive\Documents\MDICC\MDICC-main\MDICC-main\kirc\nmi_scores_kirc.csv'); 

%% Set precision

precision = 0.0:0.0001:0.04;

%% Calculate the maximum of ARI and NMI and determine the corresponding precision value

% LIHC
% find max values of ARI and NMI and their indices
[max_ari_lihc, idx1_lihc] = max(ari_lihc);
[max_nmi_lihc, idx2_lihc] = max(nmi_lihc);

if idx1_lihc == idx2_lihc % check if the indices of max ARI and NMI are equal
    optimal_precision_lihc = precision(idx1_lihc); % determine the optimal precision value
    disp("LIHC - Precision of max ARI and NMI score is " + optimal_precision_lihc)
else
    disp("Maximum values do not have the same index.")
    optimal_precision_lihc = precision(idx1_lihc);
end

% write optimal scores and corresponding precision to file
writematrix([max_ari_lihc, max_nmi_lihc, optimal_precision_lihc],'C:\Users\floor\OneDrive\Documents\MDICC\MDICC-main\MDICC-main\lihc\opt_prec_scores_lihc_ds.txt','Delimiter','space')

% KIRC
% find max values of ARI and NMI and their indices
[max_ari_kirc, idx1_kirc] = max(ari_kirc);
[max_nmi_kirc, idx2_kirc] = max(nmi_kirc);

if idx1_kirc == idx2_kirc % check if the indices of max ARI and NMI are equal
    optimal_precision_kirc = precision(idx1_kirc); % determine the optimal precision value
    disp("KIRC - Precision of max ARI and NMI score is " + optimal_precision_kirc)
else
    disp("Maximum values do not have the same index.")
    optimal_precision_kirc = precision(idx1_kirc);
end

% write optimal scores and corresponding precision to file
writematrix([max_ari_kirc, max_nmi_kirc, optimal_precision_kirc],'C:\Users\floor\OneDrive\Documents\MDICC\MDICC-main\MDICC-main\kirc\opt_prec_scores_kirc_ds.txt','Delimiter','space')

%% Calculate average of optimal precision values

opt_precs = [optimal_precision_lihc, optimal_precision_kirc];

% compute mean and round average to 4 decimals
avg_precision = round(mean(opt_precs), 4);

disp("The average precision is: " + avg_precision)

% find index average precision, account for floating-point numbers
if isfloat(avg_precision)
    tol = eps * 2; % floating point number tolerance
    idx_avg_prec = find(abs(precision - avg_precision) <= tol);
else
    idx_avg_prec = find(precision == avg_precision);
end

%% Compute ARI and NMI score per data set for the average precision

% LIHC
% find precision value of average precision index
ap_ari_lihc = ari_lihc(idx_avg_prec);
ap_nmi_lihc = nmi_lihc(idx_avg_prec);

disp("lihc: ")
disp([ap_ari_lihc ap_nmi_lihc])

% write ARI and NMI for average precision to file
fileID = fopen('C:\Users\floor\OneDrive\Documents\MDICC\MDICC-main\MDICC-main\lihc\avg_prec_lihc_ds.txt', 'w');
formatSpec = 'average precision is %6.4f \nari score is %6.8f \nnmi score is %6.8f';
fprintf(fileID, formatSpec, avg_precision, ap_ari_lihc, ap_nmi_lihc);

% KIRC
% find precision value of average precision index
ap_ari_kirc = ari_kirc(idx_avg_prec);
ap_nmi_kirc = nmi_kirc(idx_avg_prec);

disp("kirc: ")
disp([ap_ari_kirc, ap_nmi_kirc])

% write ARI and NMI for average precision to file
fileID = fopen('C:\Users\floor\OneDrive\Documents\MDICC\MDICC-main\MDICC-main\kirc\avg_prec_kirc_ds.txt', 'w');
formatSpec = 'average precision is %6.4f \nari score is %6.8f \nnmi score is %6.8f';
fprintf(fileID, formatSpec, avg_precision, ap_ari_kirc, ap_nmi_kirc);

%% Plot ARI and NMI scores against precision
% plot includes average precision and optimal precision
% plot includes scores of MDICC with kmeans

fig = figure;

% LIHC
kmeans_lihc = readmatrix('C:\Users\floor\OneDrive\Documents\MDICC\MDICC-main\MDICC-main\lihc\scores_km_lihc.csv');
kmeans_lihc = kmeans_lihc(2:3,2); % remove other scores so that only ARI and NMI scores remain

subplot(2,1,1);
hold on
plot(precision, ari_lihc, 'b')
plot(precision, nmi_lihc, 'r')

% add vertical lines (average and optimal precision) to plot
xl1 = xline(avg_precision,'-','Average');
xl2 = xline(optimal_precision_lihc,'-','Optimal');
xl1.LabelVerticalAlignment = 'bottom';
xl1.LabelHorizontalAlignment = 'left';
xl2.LabelVerticalAlignment = 'bottom';
xl2.LabelHorizontalAlignment = 'left';

% add horizontal lines (scores using kmeans) to plot
yl1 = yline(kmeans_lihc(1,1), 'b', 'ARI Kmeans');
yl2 = yline(kmeans_lihc(2,1), 'r', 'NMI Kmeans');
yl1.LabelVerticalAlignment = 'top';
yl1.LabelHorizontalAlignment = 'left';
yl2.LabelVerticalAlignment = 'bottom';
yl2.LabelHorizontalAlignment = 'left';

legend('ARI DS', 'NMI DS', 'Location', 'southwest')
title('LIHC')
hold off

% set size of plot
x=10;
y=10;
width=550;
height=1000;

set(gcf, 'position', [x,y,width,height])
ylim([0, 1])

% KIRC
kmeans_kirc = readmatrix('C:\Users\floor\OneDrive\Documents\MDICC\MDICC-main\MDICC-main\kirc\scores_km_kirc.csv');
kmeans_kirc = kmeans_kirc(2:3,2); % remove data so that only ARI and NMI score remain

subplot(2,1,2)
hold on
plot(precision, ari_kirc, 'b')
plot(precision, nmi_kirc, 'r')

% add vertical lines (average and optimal precision) to plot
xl1 = xline(avg_precision,'-','Average');
xl2 = xline(optimal_precision_kirc,'-','Optimal');
xl1.LabelVerticalAlignment = 'bottom';
xl1.LabelHorizontalAlignment = 'left';
xl2.LabelVerticalAlignment = 'bottom';
xl2.LabelHorizontalAlignment = 'left';

% add horizontal lines (scores using kmeans) to plot
yl1 = yline(kmeans_kirc(1,1), 'b', 'ARI Kmeans');
yl2 = yline(kmeans_kirc(2,1), 'r', 'NMI Kmeans');
yl1.LabelVerticalAlignment = 'top';
yl1.LabelHorizontalAlignment = 'left';
yl2.LabelVerticalAlignment = 'bottom';
yl2.LabelHorizontalAlignment = 'left';

legend('ARI DS', 'NMI DS', 'Location', 'southwest')
title('KIRC')
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

sgtitle('Comparison of original MDICC to MDICC using DS clustering')

% save figure as pdf file
exportgraphics(fig,'C:\Users\floor\OneDrive\Documents\MDICC\MDICC-main\MDICC-main\comparisons\fs_to_fs_ds.pdf','ContentType','vector')
