function is_C_H = map_high_conf_to_rxns(model, GPRmat, GPRrxns, C_H_genes)

C_H_GPR = double(ismember(GPRmat, C_H_genes));
C_H_GPR(cellfun('isempty', GPRmat)) = nan;
C_H_GPR_min = min(C_H_GPR, [], 2);
is_C_H = zeros(size(model.rxns));
for i = 1:numel(model.rxns)
    rxn_GPRs = strmatch(model.rxns{i}, GPRrxns, 'exact');
    if numel(rxn_GPRs)
        is_C_H(i) = max(C_H_GPR_min(rxn_GPRs));
    end
end
is_C_H = logical(is_C_H);
