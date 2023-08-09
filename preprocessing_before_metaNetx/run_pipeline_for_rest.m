%% Run pipeline for preprocessing MetaNetX 
function [model_after_rules, deleted_reactions] = run_pipeline_for_rest(old_model)
%% Add reversible/irreversible reaction status for model "rev" field 
% add rev field in the main model 
rev=createRev(old_model);
old_model.rev=rev;

%% remove bad reactions
% This function removes reactions enable production/consumption of metabolites even when exchange reactions aren't used.
% this is a particular demand of MetaNetX for these reactions to be
% removed. They may be important for other analysis 

[drxns, deleted_reactions]= removeBadRxns(old_model, 3);

%% Solve imbalance btw dimensions of drxns and rules field 
% remove certain reactions from the model doesn't automatically remove
% rules regarding those reactions. Run this function to solve this problem 
[model_after_rules]= create_rules(old_model, drxns);

