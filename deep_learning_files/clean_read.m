%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Training for Deep Learning based APT imaging using partially
% synthetic data
%
% Authors: Malvika Viswanathan, Zhongliang Zu
%
%
%  Objective: DL model to be trained for predicting amplitude and width
%  from input Z spectrum. 
%
% Please contact zhongliang.zu@vumc.org incase you have any doubts with the
% code.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;
tic

% Upload noisy dataset for robustness // requires input and target files
% Requires Deep Learning Toolbox 
load(' ');
% Update: use noisy data with different levels of SNR for testing data
matrix_inp = matrix_input;

%Target matrix
tar3 = matrix_MTR_output1.';
rawDataTarget = tar3;
rawDataTrain = matrix_inp;
n_timesteps = size(matrix_inp,2);
n_features = size(matrix_inp,1);
n_outputs = size(tar3,1);
idx = 1:n_timesteps;    

%Training and Validation Data
DataTrain = rawDataTrain(:,idx(1:round(n_timesteps*0.75)));  
ValidDataTrain = rawDataTrain(:,idx(round(n_timesteps*0.75)+1:end));
DataTarget = [rawDataTarget(1, idx(1:round(n_timesteps*0.75))); ...
    (rawDataTarget(2,  idx(1:round(n_timesteps*0.75))))/150]; % scale target because of high SD
ValidDataTarget = [rawDataTarget(1, idx(round(n_timesteps*0.75)+1:end));...
    (rawDataTarget(2,  idx(round(n_timesteps*0.75)+1:end)))/150]; % scale target because of high SD
% rescale width after prediction --> Predicted DL width*150

% Defining the model.
layers = [...
    sequenceInputLayer(n_features, Normalization="zerocenter")

    fullyConnectedLayer(100) 
    reluLayer

    fullyConnectedLayer(100)
    reluLayer 


    fullyConnectedLayer(n_outputs)
    regressionLayer %MSE Loss
    
    ];

% most of these hyper parameters you have to try and choose the best one
% according to your dataset. 
% if you get spike or Nan training, lower the Initial Learn Rate.
% Hyper-parameter selection showed in Supporting information
options = trainingOptions('adam', ...
    MaxEpochs = 4000, ...
    miniBatchSize = 64,...
    InitialLearnRate = 0.001, ...
    ValidationData = {ValidDataTrain, ValidDataTarget}, ...
    ValidationFrequency = 100, ...
    Shuffle='every-epoch',...
    Plots = 'training-progress');

net = trainNetwork(DataTrain, DataTarget, layers, options);	
toc

% Use predict() to perform prediction on testing data(from tissue mimicking
% folder or in vivo Z spectrum
%
% input = testing_input;
% predicted = predict(net, input);
