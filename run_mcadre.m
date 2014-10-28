
clear
clc

%% Load input data
load('U87mCADREInputs');
load('HR1_CbModel');
load('confidenceScores');
load('precursorMets');
changeCobraSolver('glpk');

%% Set parameters for mcadre

method = 1; % fastFVA
% method = 2; % fastcc

%% Run mcadre

[GM, C, NC, PM, Z, model_C, pruneTime, cRes] = ...
    mcadre(model, G, U, confidenceScores, C_H_genes, method)
