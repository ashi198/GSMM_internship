%% initialize the toolbox  
initCobraToolbox('false');
solverOK=changeCobraSolver('gurobi','LP');

% RAVEN toolbox's fillGap function was used to find all blocked reactions and dead-end metabolites. 
% KEGG, a reference database, helps to link these blocked reactions and fill in the gaps. 
% The fillGap function also checks model connectivity with the gapReport function. 

% import models
model_names = {'raven.xml'; '1A01_carveMe.xml'; 'ModelSeed_1A01.sbml'; '1A01_Kbase.sbml'};
for i =1:length(model_names)
    models{i, 1}= importModel(model_names{i});
end 

%% find the number of blocked reactions and deadend metabolites 

for j= 1:length(models)
BlockedReaction{j} = findBlockedReaction(models{j});
blockedIrr{j}=findBlockedIrrRxns(models{j});
mets{j} = detectDeadEnds(models{j});
end  
% for raven = 1273; CarveMe = 86; ModelSEED: 653; Kbase: 573

%% only carveMe and RAVEN models are used for gap-filling
%import reference KEGG model 
keggModel=getModelFromKEGG([],false,false,false,false);

% The KEGG model is associated to more than 6,400,000 genes. They will not
% be used for the gapfilling, so they are removed to make this a little
% faster
keggModel=rmfield(keggModel,'genes');
keggModel=rmfield(keggModel,'rxnGeneMat');

%% GAP FILLING FOR RAVEN
% It is already known from the previous commands that there are some
% unbalanced reactions in KEGG. Only use the balanced ones for the gap
% filling
balanceStructure=getElementalBalance(keggModel);
keggModel=removeReactions(keggModel,balanceStructure.balanceStatus~=1,true,true);

% The function fillGaps with these settings will try to include reactions in
% order to have flux through all reactions in the model. There are other
% settings as well. The first flag says that production of all metabolites
% should be allowed.
params.relGap=0.6; %Lower number for a more exhaustive search
params.printReport=true;

[newConnected, cannotConnect, addedRxns, newModel, exitFlag]=fillGaps(models{1},keggModel,true,false,false,[],params);
writeCbModel(model, 'format', 'sbml')

%% PLEASE NOTE: THE OUTPUT OF THIS FILE WILL GIVE YOU A MODEL WHICH MIGHT NECESSARILY NOT ADHERE TO THE SBML STANDARDS. VALIDATE THE MODEL FIRST USING
% SBML VALIDATOR AND THEN DO FURTHER ANALYSIS
% (https://synonym.caltech.edu/validator_servlet/)

%% FURTHER MODEL IMPROVMENT (MIGHT TAKE SEVERAL HOURS TO RUN)
% Continue to improve the connectivity of the model by identifying
% metabolites that should be connected. A convenient way to get an overview
% of how connected the model is, and at the same time getting a lot of
% useful data, is to use gapReport. Note that it can take several to many
% hours to run, depending on the number of gaps in the model.
[noFluxRxns, noFluxRxnsRelaxed, subGraphs, notProducedMets, minToConnect,...
    neededForProductionMat]=gapReport(newModel);

% Rerun gapReport and use the output for targeting the gap-filling efforts.
% Note that only some info is printed; most of it is available in the output
% structures. Work like this in an iterative manner until the model is of
% sufficient quality.