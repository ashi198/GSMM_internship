%% this code tells the similarity of the reactions and metabolites between the models processed by MetaNetX. The original curated model,
% and models from Kbase, ModelSEED, CarveMe, and Raven are analysed here

% initialize cobra toolbox 
initCobraToolbox('false') %false to avoid updating the toolbox 
solverOK=changeCobraSolver('gurobi','LP'); % set the solver here 

%% upload models 
% keep the original (true) model in the first position
model_names= {'original_model.xml'; 'carve_A1.xml'; 'ModelSeed_A1.xml';....
    'Kbase_A1.xml'; 'raven_refined.xml'}; % change this with the names of your models

% stores all the models in the 'models' cell array 
for i=1:length(model_names)
    models{i, 1}= readCbModel(model_names{i});
end 


%% extract IDs, reaction equations, and metabolite formulas   

for i=1:length(models)
    reactions_ids{i,1}=models{i,1}.rxns;
    rxn_formulas{i,1}=printRxnFormulaOri(models{i});
    metabolite_ids{i,1}=models{i,1}.mets;
    metabolite_formulas{i,1}=models{i,1}.metFormulas;
end 

%% Find the common reactions/metabolites between the original and draft models...
% This will result in a array called "common_rxns_metabolites. This has 4
% rows and 5 columns. The rows denote 'common rxn ids';'common rxn
% formulas'; 'common metabolite ids'; 'common metabolite formulas' between
% the original model and the draft model. The column denotes draft models
% made from carveMe, Kbase, ModelSEED, RAVEN. 

for i=1:length(reactions_ids)
    for j=1:5
       if j~=1 %exclude comparing original model with itself
    common_rxns_metabolites{1, j} = reactions_ids{1}(ismember(reactions_ids{1}, reactions_ids{j}));
    common_rxns_metabolites{2, j} = rxn_formulas{1}(ismember(rxn_formulas{1}, rxn_formulas{j}));
    common_rxns_metabolites{3, j} = metabolite_ids{1}(ismember(metabolite_ids{1}, metabolite_ids{j}));
    common_rxns_metabolites{4, j} = metabolite_formulas{1}(ismember(metabolite_formulas{1}, metabolite_formulas{j}));
        else 
    common_rxns_metabolites{1, j} = [];
    common_rxns_metabolites{2, j} = [];
    common_rxns_metabolites{3, j} = [];
    common_rxns_metabolites{4, j} = [];
        end 
    end 
end 

%% calculate true positives, false positives and negatives. 
%The rows denote 'common rxn ids';'common rxn
% formulas'; 'common metabolite ids'; 'common metabolite formulas' between
% the original model and the draft model. The column (starting from 2) denotes draft models
% made from carveMe, Kbase, ModelSEED, RAVEN.
% Initialize counters for TP, FP, and FN

TP = 0;
FP = 0;
FN = 0;

% Iterate over the columns of "common_rxns_metabolites" excluding the first column 
for j = 2:5 %change this with the number of columns in common_rxns_metabolites
    for i= 1:4 % change this with appropriate number of models
    % Calculate true positives
    TP(i, j) = numel(common_rxns_metabolites{i, j});
    TP(i, j) = numel(common_rxns_metabolites{i, j});
    TP(i, j) = numel(common_rxns_metabolites{i, j});
    TP(i, j) = numel(common_rxns_metabolites{i, j});
    end 
end

% for FPs 
for j = 2:5
    FP (1, j)= numel(setdiff(reactions_ids{j}, reactions_ids{1}));
    FP (2, j)=numel(setdiff(rxn_formulas{j}, rxn_formulas{1}));
    FP (3, j)=numel(setdiff(metabolite_ids{j}, metabolite_ids{1}));
    FP (4, j)=numel(setdiff(metabolite_formulas{j}, metabolite_formulas{1}));
end 
    
% for FNs
for j = 2:5
    FN (1, j)= numel(setdiff(reactions_ids{1}, reactions_ids{j}));
    FN (2, j)=numel(setdiff(rxn_formulas{1}, rxn_formulas{j}));
    FN (3, j)=numel(setdiff(metabolite_ids{1}, metabolite_ids{j}));
    FN (4, j)=numel(setdiff(metabolite_formulas{1}, metabolite_formulas{j}));
end 
    
% create confusion matrix 
% Calculate precision
precision= 0;
recall= 0;
f1_score= 0;
for i = 2:5
 for j= 1: 4
precision (j, i) = TP(j, i) / (TP(j, i) + FP(j, i));
% Calculate recall
recall(j, i) = TP(j, i) / (TP(j, i) + FN(j, i));

% Calculate F1 score
 f1_score (j, i)= 2 * (precision(j, i) * recall(j, i)) / (precision(j, i) + recall(j, i));
 end 
end 

%% Plot Venn diagrams for analysis
% download the add on Venn Euler diagram before executing this code 
%https://www.mathworks.com/matlabcentral/fileexchange/98974-venn-euler-diagram

Data {1, 1}= {reactions_ids{2}, reactions_ids{3}, reactions_ids{4}, reactions_ids{5}};
Data {2, 1}= {rxn_formulas{2}, rxn_formulas{3}, rxn_formulas{4}, rxn_formulas{5}};
Data {3, 1}= {metabolite_ids{2}, metabolite_ids{3}, metabolite_ids{4}, metabolite_ids{5}};
Data {4, 1}= {metabolite_formulas{2}, metabolite_formulas{3}, metabolite_formulas{4}, metabolite_formulas{5}};
setName={'CarveMe', 'ModelSeed', 'Kbase', 'Raven'};
h = vennEulerDiagram(Data{4, 1}, setName, 'drawProportional', true, 'PositionConstraint', 'outerposition');
h.ShowIntersectionCounts = true;


