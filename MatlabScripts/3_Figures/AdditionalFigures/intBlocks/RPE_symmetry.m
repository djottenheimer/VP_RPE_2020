Colors = load_colors();
load('ModData_intBlocks.mat');
load ('R_intBlocks.mat');
load ('intBlocks_MLEfits.mat');
bm_RD=select_RPEmods(os, 'RD', 'particularModels', {'mean','curr','base'},'plotmodels_Flag',false);
    
masks=cat(2,bm_RD.mask_base',bm_RD.mask_curr',bm_RD.mask_mean');
region = strcmp(R_blocks.Region,'VP');
task = R_blocks.Blocks==0;

figure;

%% activity plots
ymax=4;
alph=0.3;
Sel=region&task&masks(:,1);

neuronIDs = find(Sel);
pos_error=[];
neg_error=[];
for neuron = 1:sum(Sel)
    neuronID = neuronIDs(neuron);
    
    SSHz = mean(CS.RDHz{neuronID}(CS.Predictors{neuronID}(:,1)==1 & CS.Predictors{neuronID}(:,2)==1));
    SMHz = mean(CS.RDHz{neuronID}(CS.Predictors{neuronID}(:,1)==1 & CS.Predictors{neuronID}(:,2)==0));
    pos_error(neuron,1) = (SMHz - SSHz) / std(CS.RDHz{neuronID});
    
    MSHz = mean(CS.RDHz{neuronID}(CS.Predictors{neuronID}(:,1)==0 & CS.Predictors{neuronID}(:,2)==1));
    MMHz = mean(CS.RDHz{neuronID}(CS.Predictors{neuronID}(:,1)==0 & CS.Predictors{neuronID}(:,2)==0));
    neg_error(neuron,1) = (MSHz - MMHz) / std(CS.RDHz{neuronID});
    
end
    
subplot(1,1,1);
hold on;
scatter(pos_error,neg_error);
plot([-1.5 1.5],[0 0],'color','k','linewidth',1);
plot([0 0],[-1.5 1.5],'color','k','linewidth',1);
xlabel('(suc after mal) - (suc after suc)');
ylabel('(mal after suc) - (mal after mal)');
title('Negative prediction error vs. positive prediction error');


