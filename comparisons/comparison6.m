%% Comparisons of performance original MDICC and our approach
% our approach includes adding feature selection and DS clustering to
% original MDICC

close all
clear all

clc

%% Create a new figure
fig = figure;

%% MDICC with and without feature selection

% LIHC
mdicc_lihc = readmatrix('C:\Users\floor\OneDrive\Documents\MDICC\MDICC-main\MDICC-main\lihc\scores_km_lihc.csv');
mdicc_fsvt_lihc = readmatrix('C:\Users\floor\OneDrive\Documents\MDICC\MDICC-main\MDICC-main\lihc\featureselection\scores_km_lihc_fsvt.csv'); 
mdicc_fsk_lihc = readmatrix('C:\Users\floor\OneDrive\Documents\MDICC\MDICC-main\MDICC-main\lihc\featureselection\scores_km_lihc_fsk.csv'); 

mdicc_lihc = mdicc_lihc(2:3,2);
mdicc_fsvt_lihc = mdicc_fsvt_lihc(2:3,2);
mdicc_fsk_lihc = mdicc_fsk_lihc(2:3,2);

X = categorical({'ARI', 'NMI'});
y_lihc = [mdicc_lihc mdicc_fsvt_lihc mdicc_fsk_lihc];

% set size of plott
x=10;
y=10;
width=1000;
height=800;

subplot(3,2,1);
b = bar(X, y_lihc);
title('LIHC')
set(gcf, 'position', [x,y,width,height])

ylim([0, 1])
ylabel({'Feature Selection   '}, 'Rotation' , 0, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right');

% KIRC
mdicc_kirc = readmatrix('C:\Users\floor\OneDrive\Documents\MDICC\MDICC-main\MDICC-main\kirc\scores_km_kirc.csv');
mdicc_fsvt_kirc = readmatrix('C:\Users\floor\OneDrive\Documents\MDICC\MDICC-main\MDICC-main\kirc\featureselection\scores_km_kirc_fsvt.csv'); 
mdicc_fsk_kirc = readmatrix('C:\Users\floor\OneDrive\Documents\MDICC\MDICC-main\MDICC-main\kirc\featureselection\scores_km_kirc_fsk.csv'); 

mdicc_kirc = mdicc_kirc(2:3,2);
mdicc_fsvt_kirc = mdicc_fsvt_kirc(2:3,2);
mdicc_fsk_kirc = mdicc_fsk_kirc(2:3,2);

y_kirc = [mdicc_kirc mdicc_fsvt_kirc mdicc_fsk_kirc];

subplot(3,2,2);
b = bar(X, y_kirc);
legend('Original', 'VarianceThreshold', 'SelectKBest', 'Location', 'southeast')
title('KIRC')

set(gcf, 'position', [x,y,width,height])
ylim([0, 1])

%% MDICC using Kmeans and using DS clustering

% LIHC
mdicc_lihc = readmatrix('C:\Users\floor\OneDrive\Documents\MDICC\MDICC-main\MDICC-main\lihc\scores_km_lihc.csv');
ds_lihc = readmatrix('C:\Users\floor\OneDrive\Documents\MDICC\MDICC-main\MDICC-main\lihc\opt_prec_scores_lihc_ds.txt');

mdicc_lihc = mdicc_lihc(2:3,2);
ds_lihc = ds_lihc(:,1:2)';

y_lihc = [mdicc_lihc ds_lihc];

subplot(3,2,3);
b = bar(X, y_lihc);

set(gcf, 'position', [x,y,width,height])
ylabel({'Clustering           '}, 'Rotation' , 0, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right');
ylim([0, 1])

% KIRC
mdicc_kirc = readmatrix('C:\Users\floor\OneDrive\Documents\MDICC\MDICC-main\MDICC-main\kirc\scores_km_kirc.csv');
ds_kirc = readmatrix('C:\Users\floor\OneDrive\Documents\MDICC\MDICC-main\MDICC-main\kirc\opt_prec_scores_kirc_ds.txt');

mdicc_kirc = mdicc_kirc(2:3,2);
ds_kirc = ds_kirc(:,1:2)';

y_kirc = [mdicc_kirc ds_kirc];

subplot(3,2,4);
b = bar(X, y_kirc);
legend('Original', 'DS clustering', 'Location', 'southeast')

set(gcf, 'position', [x,y,width,height])
ylim([0, 1])

%% Original MDICC vs our approach

% LIHC - original MDICC vs our approach
mdicc_lihc = readmatrix('C:\Users\floor\OneDrive\Documents\MDICC\MDICC-main\MDICC-main\lihc\scores_km_lihc.csv');
ds_fsvt_lihc = readmatrix('C:\Users\floor\OneDrive\Documents\MDICC\MDICC-main\MDICC-main\lihc\opt_prec_scores_lihc_ds_fsvt.txt'); 
ds_fsk_lihc = readmatrix('C:\Users\floor\OneDrive\Documents\MDICC\MDICC-main\MDICC-main\lihc\opt_prec_scores_lihc_ds_fsk.txt'); 

mdicc_lihc = mdicc_lihc(2:3,2);
ds_fsvt_lihc = ds_fsvt_lihc(:,1:2)';
ds_fsk_lihc = ds_fsk_lihc(:,1:2)';

X = categorical({'ARI', 'NMI'});
y_lihc = [mdicc_lihc ds_fsvt_lihc ds_fsk_lihc];

subplot(3,2,5);
b = bar(X, y_lihc);

set(gcf, 'position', [x,y,width,height])
ylabel({'Final Model   '}, 'Rotation' , 0, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right');
ylim([0, 1])

% KIRC - original MDICC vs our approach
mdicc_kirc = readmatrix('C:\Users\floor\OneDrive\Documents\MDICC\MDICC-main\MDICC-main\kirc\scores_km_kirc.csv');
ds_fsvt_kirc = readmatrix('C:\Users\floor\OneDrive\Documents\MDICC\MDICC-main\MDICC-main\kirc\opt_prec_scores_kirc_ds_fsvt.txt'); 
ds_fsk_kirc = readmatrix('C:\Users\floor\OneDrive\Documents\MDICC\MDICC-main\MDICC-main\kirc\opt_prec_scores_kirc_ds_fsk.txt');  

mdicc_kirc = mdicc_kirc(2:3,2);
ds_fsvt_kirc = ds_fsvt_kirc(:,1:2)';
ds_fsk_kirc = ds_fsk_kirc(:,1:2)';

y_kirc = [mdicc_kirc ds_fsvt_kirc ds_fsk_kirc];

subplot(3,2,6);
b = bar(X, y_kirc);

legend('Original', 'VarianceThreshold DS', 'SelectKBest DS', 'Location', 'southeast')
set(gcf, 'position', [x,y,width,height])
ylim([0, 1])

% set x label, y label and title for plot
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.YLabel.Visible='on';
han.XLabel.Visible='on';
ylabel(han,'evaluation score');
xlabel(han,'performance metric');
sgtitle('Comparison of original MDICC to adapted MDICC')

% save figure as pdf file
exportgraphics(fig,'C:\Users\floor\OneDrive\Documents\MDICC\MDICC-main\MDICC-main\comparisons\comparison6.pdf','ContentType','vector')
