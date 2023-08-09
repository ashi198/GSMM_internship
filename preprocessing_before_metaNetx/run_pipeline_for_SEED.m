%% Run pipeline for preprocessing MetaNetX
% this pipeline is specifically for models derived from SEED database (Kbase, modelSEED)
function [model_after_rules, deleted_reactions] = run_pipeline_for_SEED(old_model)

%% Add reversible/irreversible reaction status for model "rev" field 
% add rev field in the main model 
rev=createRev(old_model);
old_model.rev=rev;

%% Add metabolite formulas 
% Models created using SEED database (ModelSeed and Kbase) files do not contain the following fields:
% MetFormulas. For this purpose, the compound .tsv file was downloaded from (https://github.com/ModelSEED/ModelSEEDDatabase/tree/master/Biochemistry)
% and use for extracting metabolite formulas

compounds = readtable('output.csv');

% use createMetaboliteFormulas to add metabolite formulas 
[model_with_formulas, ~] = createMetaboliteFormulas(old_model, compounds);

%% remove bad reactions
% This function removes reactions enable production/consumption of metabolites even when exchange reactions aren't used.
% this is a particular demand of MetaNetX for these reactions to be
% removed. They may be important for other analysis 

[drxns, deleted_reactions]= removeBadRxns(model_with_formulas, 3);

%% Solve imbalance btw dimensions of drxns and rules field 
% remove certain reactions from the model doesn't automatically remove
% rules regarding those reactions. Run this function to solve this problem 
[model_after_rules]= create_rules(model_with_formulas, drxns);

