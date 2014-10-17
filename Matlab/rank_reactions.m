function [GM,C,NC,P,Z,inactiveRxns,model_C] = rank_reactions(model, G, U, confidenceScores, C_H_genes)

% Inputs:
% - model
% - gene IDs from expression data
% - gene ubiquity scores (i.e., from mas5callToExpression)
% - C_H_genes: high confidence genes (optional)

% Outputs:
% - GM: generic model with inactive reactions removed
% - C: core reactions
% - NC: non-core reactions
% - P: removal order of non-core reactions


%% Parse GPRs
[GPRrxns, GPRmat] = parse_gprs(model);


%% Map high confidence genes to reactions
if nargin > 4
    is_C_H = map_high_conf_to_rxns(model, GPRmat, GPRrxns, C_H_genes);
else is_C_H = false(size(model.rxns));
end


%% Map gene ubiqiuty scores to reactions
U_GPR = map_gene_scores_to_rxns(model, G, U, GPRrxns, GPRmat);


%% Determine confidence level-based evidence
E_L = confidenceScores;
% Reactions with no confidence information should not be ranked higher than
% those with non-zero confidence
E_L(isnan(E_L)) = 0;


%% Calculate expression-based evidence
E_X = calc_expr_evidence(model, GPRrxns, U_GPR, is_C_H);

C = model.rxns(E_X > 0.9);
model_C = C;


%% Find inactive/blocked reactions
% I'm not 100% convinced that we should do this before calculating the
% connectivity-based evidence, but I'll look into this more later.
method = 1; % fastFVA
inactiveRxns = check_model_consistency(model, method);

inactive_C = intersect(inactiveRxns, C);
C = setdiff(C, inactiveRxns);

GM = removeRxns(model,inactiveRxns);
R_G = GM.rxns; [NC, NC_idx] = setdiff(R_G, C);

% Update list of expression-based evidence
[rxnsInt, GM_idx, model_idx] = intersect(GM.rxns, model.rxns);
E_X_GM = zeros(size(GM.rxns)); E_X_GM(GM_idx) = E_X(model_idx);
E_X = E_X_GM;

% Update list of confidence level-based evidence
E_L_GM = zeros(size(GM.rxns));
E_L_GM(GM_idx) = E_L(model_idx);
E_L = E_L_GM;


%% Calculate connectivity-based evidence
E_C = calc_conn_evidence(GM);


%% Rank non-core reactions
E_X_NC = E_X(NC_idx); E_C_NC = E_C(NC_idx); E_L_NC = E_L(NC_idx);
[E_NC,NC_order] = sortrows([E_X_NC, E_C_NC, E_L_NC], [1, 2, 3]);
P = NC(NC_order);


%% Identify zero-expression reactions
Z = P(E_NC(:,1) == -1e-6);
