%% RAVEN RECONSTRUCTION 

% This code creates a model for V.splendidus. The parameters are set to exclude
% general or unclear reactions and reactions with undefined stoichiometry.
% Type "help getKEGGModelForOrganism" to see what the different parameters
% are for. This process takes up to 20-35 minutes in macOS, Unix systems and
% 40-55 minutes in Windows, depending on your hardware and the size of
% target organism proteome
% the .fa was exported from https://www.ncbi.nlm.nih.gov/genome/933?genome_assembly_id=348374

%% initialize the toolbox  
initCobraToolbox('false');
solverOK=changeCobraSolver('gurobi','LP');

%%
model=getKEGGModelForOrganism('vsl','vsl.fna','prok90_kegg102','output',false,false,false,false,10^-30,0.8,0.3,-1);

% The resulting model should contain around 2156 reactions, 2276
% metabolites and 937 genes. Small variations are possible since it is an
% heuristic algorithm and different KEGG versions will give slightly
% different results.
disp(model);

% A first control is that the model should not be able to produce any
% metabolites without uptake of some metabolites. This commonly happens when
% metabolites have a different meaning in different reactions. The best way
% to find such reactions is to run makeSomething and analyze the resulting
% solution for bad reactions. An automated approach is to use removeBadRxns,
% which tries to do the same thing in an automated manner. Type
% "help removeBadRxns" for details.
[newModel, removedRxns]=removeBadRxns(model);

% One can see an error about that H+ can be made even if no reactions were
% unbalanced. Protons are particularly problematic since it is rather
% arbitary at which pH the formulas are written for. For the purpose of this
% analysis, the protons can be ignored and fixed later.
[newModel, removedRxns]=removeBadRxns(model,1,{'H+'},true);

% Only one reaction was removed because it enabled the model to produce
% something from nothing. Since it is only one reaction, it might be
% worthwhile to look into this in more detail.
disp(removedRxns);

% According to the information in KEGG about this reaction is
% for 1,4-alpha-D-Glucan maltohydrolase. One might want to look at the flux distributions
% in more detail to try to find out if there is any better alternative to
% delete. Use makeSomething to do this.
[fluxes, metabolite]=makeSomething(model,{'H+'},true);
model.metNames(metabolite)

% The model could produce H2O using the following reactions
printFluxes(model, fluxes, false, [], [],'%rxnID (%rxnName):\n\t%eqn: %flux\n')

% That resulted in quite a lot of fluxes to look through. Maybe it is easier
% if the elementally balanced reactions are excluded.
balanceStructure=getElementalBalance(model);

% The hydrogen balancing is a bit tricky, so this time it seems more
% relevant to only look at the ones unbalanced for oxygen (since water was
% produced)
goodOnes=balanceStructure.leftComp(:,6)==balanceStructure.rightComp(:,6);
printFluxes(removeReactions(model,goodOnes), fluxes(~goodOnes), false, [], [],'%rxnID (%rxnName):\n\t%eqn: %flux\n');

% There is still a good number of reactions. Leave only the reactions which
% involve amylose or starch, from one of the problematic reactions that was
% identified earlier.
printFluxes(model, fluxes, false, [], [],'%rxnID (%rxnName):\n\t%eqn: %flux\n',{'Amylose';'Starch'});

% There are two elementally unbalanced reactions, including, the reaction
% which was identified by removeBadRxns. When looking to these reactions
% closer, one can notice the contradiction between the reactions. The first
% one shows that amylose and starch are interconvertible, while the second
% reaction shows that amylose contains one less glucose unit than starch.
% This type of general reactions are problematic and should be fixed
% manually. It is hereby chosen to trust removeBadRxns and delete R02110.
model=removeReactions(model,'R02112');

% The model can no longer make something from nothing. Can it consume
% something without any output?
[solution, metabolite]=consumeSomething(model,{'H+'},true);
model.metNames(metabolite)

% Nope, so that was good. Add some uptakes and see what it can produce.
[~, J]=ismember({'D-Glucose';'H2O';'Orthophosphate';'Oxygen';'NH3';'Sulfate'},model.metNames);
[model, addedRxns]=addExchangeRxns(model,'in',J);

% Check which metabolites can be produced given these uptakes. The
% canProduce function allows for output of all metabolites. This will not
% happen in the real cell, but it is very useful for functionality testing
% of the model. Once it is functional, the excretion reactions based on
% evidence can be added as well.
I=canProduce(model);

fprintf('%d%%\n', round(sum(I)/numel(model.mets)*100));
% It seems that around 31% of the metabolites could be synthesized. It is
% not directly clear whether this is a high or low number, many metabolites
% should not be possible to synthesize from those simple precursors.

% Try to fill gaps using the full KEGG model to see if that gives a
% significantly higher number
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
[newConnected, cannotConnect, addedRxns, newModel, exitFlag]=fillGaps(model,keggModel,true,false,false,[],params);
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
