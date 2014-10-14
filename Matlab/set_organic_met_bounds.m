function model = set_organic_met_bounds(model, exRxns)

% Identify all exchange reactions that include organic metabolites
organicExRxns = find_organic_ex_rxns(model, exRxns);

% Reset all lower bounds for organic reactions to 0 to turn off uptake
model = changeRxnBounds(model, organicExRxns, 0, 'l');

% Allow uptake of glucose and CO2
model = changeRxnBounds(model, 'EX_glc(e)', -5, 'l');
model = changeRxnBounds(model, 'EX_co2(e)', -1000, 'l');
