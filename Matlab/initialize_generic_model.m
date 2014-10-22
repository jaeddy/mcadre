function [GM, E_X, E_L] = initialize_generic_model(model, inactiveRxns, confidenceScores)

% This function creates a consistent generic model by removing all inactive
% reactions. It will also return adjusted vectors for expression-based and
% literature-based evidence, corresponding to the subset of reactions in GM.


% Define literature-based evidence from confidenceScores
if numel(confidenceScores)
    E_L = confidenceScores;
    % Reactions with no confidence information should not be ranked higher than
    % those with non-zero confidence
    E_L(isnan(E_L)) = 0;
else
    E_L = zeros(size(model.rxns));
end

GM = removeRxns(model, inactiveRxns);
R_G = GM.rxns; 
[~, NC_idx] = setdiff(R_G, C);

% Update list of expression-based evidence
[~, GM_idx, model_idx] = intersect(GM.rxns, model.rxns);
E_X_GM = zeros(size(GM.rxns)); E_X_GM(GM_idx) = E_X(model_idx);
E_X = E_X_GM;

% Update list of confidence level-based evidence
E_L_GM = zeros(size(GM.rxns));
E_L_GM(GM_idx) = E_L(model_idx);
E_L = E_L_GM;