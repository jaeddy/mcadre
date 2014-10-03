function [R_P, PM, P, NC_removed, num_removed] = update_pruned_model(R_P, PM, P, NC_removed, )

R_P = setdiff(R_P, inactive_NC);
PM = removeRxns(PM, inactive_NC);
P(ismember(P, inactive_NC)) = [];
NC_removed = NC_removed + numel(inactive_NC);
num_removed = NC_removed + C_removed;
