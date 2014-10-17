
load('U87mCADREInputs');
load('HR1_CbModel');
changeCobraSolver('glpk');

%% Test parse_gprs

% Test with different GPR types
display('## Testing parse_gprs with test reactions...');

none = 1;
iso = 25;
or = 54;
and = 805;
both = 1139;
complex = 1140;

gprTest = [none, iso, or, and, both, complex];
gprTestRxns = model.rxns(gprTest);
testModel = extractSubNetwork(model, gprTestRxns);

try
    [GPRrxns, GPRmat] = parse_gprs(testModel);
    display(['PASS...', ...
        'Function parse_gprs ran without error']);
catch err
    display(['FAIL...', ...
        'Function parse_gprs was terminated with the error:']);
    display(['> ', err.message]);
end
display('---');

% Check outputs with pre-defined expected result
display('## Checking output of parse_gprs for test reactions...');

testGPRmat = {'1591', ''; ... % iso
    '221', ''; '222', ''; ... % or
    '23657', '6520'; ... % and
    '1892', ''; '3030', '3032'; ... % both
    '549', ''; '3030', '3032'; '1892', ''}; % complex
    
if exist('GPRmat', 'var')
    if sum(any(~strcmp(GPRmat, testGPRmat))) == 0
        display(['PASS...', ...
            'Function parse_gprs returns the expected result']);
    else
        display(['FAIL...', ...
            'Function parse_gprs returns unexpected result']);
    end
else
    display(['FAIL...', ...
        'Function parse_gprs encountered error, cannot check output']);
end

%% Test map_high_conf_to_rxns

% Test with results from parse_gprs and high-confidence core from U87
display('## Testing map_high_conf_to_rxns...');

if exist('GPRmat', 'var')
    try
        is_C_H = map_high_conf_to_rxns(model, GPRmat, GPRrxns, C_H_genes);
        display(['PASS...', ...
            'Function map_high_conf_to_rxns ran without error']);
    catch err
        display(['FAIL...', ...
            'Function map_high_conf_to_rxns was terminated with the error:']);
        display(['> ', err.message]);
    end
else
    display(['FAIL...', ...
        'GPRmat missing, cannot check map_high_conf_to_rxns']);
end
display('---');


% Check outputs for small GPRmat
display('## Checking output of map_high_conf_to_rxns...');

if exist('is_C_H', 'var')
    if sum(is_C_H) == 2
        display(['PASS...', ...
            'Function map_high_conf_to_rxns returns the expected result']);
    else
        display(['FAIL...', ...
            'Function map_high_conf_to_rxns returns unexpected result']);
    end
else
    display(['FAIL...', ...
        'No output to check from function map_high_conf_to_rxns']);
end


%% Test map_gene_scores_to_rxns