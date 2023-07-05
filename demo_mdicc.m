close all
clear all

clc

%% Create the edge-weighted graph (affinity matrix)

A = readmatrix('your_path\afffus_lihc_fsvt.csv'); 
A=A.*not(eye(size(A))); % the graph should not have self−loops
A = A(2:size(A,1), 2:size(A,2));

dynType=1; % 0=Replicator Dynamics, 1=InfectionImmunization, 2=Exponential replicator dynamics
precision=0.0211; % the precision required from the dynamical system
maxIters=1000; % number of maximum iteration of the dynamical system
x=ones(size(A,1),1)./size(A,1); % starting point of the dynamical system
supportThreshold=1e-6; % threshold used to extract the support from x.

%% Clustering
% cluster for one precision value
% [C] = dominantset(A,x,supportThreshold,precision,maxIters,dynType);

% Cluster for wide range of precision values
C = zeros(size(A,1), length(0.0001:0.0001:0.04));
i = 1;

for p = 0.0001:0.0001:0.04
    precision = p; 
    C(:,i) = dominantset(A,[],[],precision,[],[]);
    i = i + 1;
end

% save all clusterings in a .txt file
writematrix(C,'your_path\full_clust_lihc_fsvt.txt','Delimiter','space')