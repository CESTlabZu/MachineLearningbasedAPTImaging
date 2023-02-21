%% Training Dataset
clear;
clc;
tic
load("train_data.mat");
k =1:39;
% Combining clean and noisy dataset for robustness
matrix_inp1 = matrix_input(k,:);
matrix_inp2 = matrix_input_noise(k,:);
matrix_inp =[matrix_inp1, matrix_inp2];

%Target matrix
tar3 = matrix_MTR_output.';
rawDataTarget = [tar3,tar3];
rawDataTrain = matrix_inp;

n_timesteps = size(matrix_inp,2);
n_features = size(matrix_inp,1);
n_outputs = size(tar3,1);
idx = 1:n_timesteps;    

%Training and Validation Data
DataTrain = rawDataTrain(:,idx(1:round(n_timesteps*0.75)));  
ValidDataTrain = rawDataTrain(:,idx(round(n_timesteps*0.75)+1:end));
DataTarget = [rawDataTarget(1, idx(1:round(n_timesteps*0.75))); ...
    log2(rawDataTarget(2,  idx(1:round(n_timesteps*0.75))))/150]; % scale target because of high SD
ValidDataTarget = [rawDataTarget(1, idx(round(n_timesteps*0.75)+1:end));...
    log2(rawDataTarget(2,  idx(round(n_timesteps*0.75)+1:end)))/150]; % scale target because of high SD

layers = [...
    sequenceInputLayer(n_features, Normalization="zerocenter")

    fullyConnectedLayer(300)
    reluLayer

    fullyConnectedLayer(300)
    reluLayer 


    fullyConnectedLayer(n_outputs)
    regressionLayer
    
    ];


options = trainingOptions('sgdm', ...
    MaxEpochs = 3000, ...
    miniBatchSize = 512,...
    InitialLearnRate = 0.000001, ...
    ValidationData = {ValidDataTrain, ValidDataTarget}, ...
    ValidationFrequency = 50, ...
    Shuffle='every-epoch',...
    Plots = 'training-progress');

net = trainNetwork(DataTrain, DataTarget, layers, options);	
toc

