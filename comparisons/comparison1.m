%% Comparison of MDICC with and without feature selection
% plot ARI and NMI scores of MDICC with and without feature selection method
% the used feature selection methods are VarianceThreshold and SelectKBest

close all
clear all

clc

%%
fig = figure;

%% LIHC
% load data
mdicc_lihc = readmatrix('C:\Users\floor\OneDrive\Documents\MDICC\MDICC-main\MDICC-main\lihc\scores_km_lihc.csv');
mdicc_fsvt_lihc = readmatrix('C:\Users\floor\OneDrive\Documents\MDICC\MDICC-main\MDICC-main\lihc\featureselection\scores_km_lihc_fsvt.csv'); 
mdicc_fsk_lihc = readmatrix('C:\Users\floor\OneDrive\Documents\MDICC\MDICC-main\MDICC-main\lihc\featureselection\scores_km_lihc_fsk.csv'); 

% select ARI and NMI score
mdicc_lihc = mdicc_lihc(2:3,2);
mdicc_fsvt_lihc = mdicc_fsvt_lihc(2:3,2);
mdicc_fsk_lihc = mdicc_fsk_lihc(2:3,2);

X = categorical({'ARI', 'NMI'});
y_lihc = [mdicc_lihc mdicc_fsvt_lihc mdicc_fsk_lihc];

% create barplot of the scores for LIHC
subplot(2,1,1);
b = bar(X, y_lihc);
legend('Original', 'VarianceThreshold', 'SelectKBest', 'Location', 'northeastoutside')
title('LIHC')

%% KIRC
% load data
mdicc_kirc = readmatrix('C:\Users\floor\OneDrive\Documents\MDICC\MDICC-main\MDICC-main\kirc\scores_km_kirc.csv');
mdicc_fsvt_kirc = readmatrix('C:\Users\floor\OneDrive\Documents\MDICC\MDICC-main\MDICC-main\kirc\featureselection\scores_km_kirc_fsvt.csv'); 
mdicc_fsk_kirc = readmatrix('C:\Users\floor\OneDrive\Documents\MDICC\MDICC-main\MDICC-main\kirc\featureselection\scores_km_kirc_fsk.csv'); 

% select ARI and NMI score
mdicc_kirc = mdicc_kirc(2:3,2);
mdicc_fsvt_kirc = mdicc_fsvt_kirc(2:3,2);
mdicc_fsk_kirc = mdicc_fsk_kirc(2:3,2);

y_kirc = [mdicc_kirc mdicc_fsvt_kirc mdicc_fsk_kirc];

% create barplot of the scores for KIRC
subplot(2,1,2);
b = bar(X, y_kirc);
legend('Original', 'VarianceThreshold', 'SelectKBest', 'Location', 'northeastoutside')
title('KIRC')

% settings for plots
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'evaluation score');
sgtitle('Comparison of original MDICC to MDICC with feature selection')

% save figure as 
exportgraphics(fig,'C:\Users\floor\OneDrive\Documents\MDICC\MDICC-main\MDICC-main\comparisons\influence_fs_omdicc.pdf','ContentType','vector')