
load precursorMets
changeCobraSolver('glpk');

%% Test set_media_ex_bounds

% mediaDef option not implemented in this version


%% Test set_organic_met_bounds


%% Test

%% Test check_model_function

% Test with default inputs
genericStatus = check_model_function(model, 'requiredMets', precursorMets);

% Test with 'biomass' flag
genericStatus = check_model_function(model, 'requiredMets', precursorMets, ...
    'biomass', 'dummy_rxn');

% Test with 'biomass' flag
genericStatus = check_model_function(model, 'requiredMets', precursorMets, ...
    'media', 'dummy_media');
