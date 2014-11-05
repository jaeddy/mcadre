
clear
clc

%% Load input data
load('testInputs');
load('humanModel');
load('precursorMets');
changeCobraSolver('glpk');

%% Set parameters for mcadre

% method = 1; % fastFVA
method = 2; % fastcc
salvageCheck = 1;

%% Run mcadre

[PM, GM, C, NC, Z, model_C, pruneTime, cRes] = ...
    mcadre(model, G, U, confidenceScores, salvageCheck, C_H_genes, method);
