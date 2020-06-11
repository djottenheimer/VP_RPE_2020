function bestModStruct = select_RPEmods(os, timePeriod, varargin)

% timePeriod should be 'RD', 'cue', or 'PE'
% scoreToUse should be 'AIC', 'BIC', or 'LR'
% pValue is used for LR tests

p = inputParser();
p.addParameter('scoreToUse', 'AIC')
p.addParameter('plotModels_Flag', true)
p.addParameter('particularModels', '')
p.addParameter('normalizeHistogram_Flag', false)
p.addParameter('pValue',0.05)
p.parse(varargin{:});

bestModStruct = struct();
bestModStruct.all_scores = [];

modToUse = ['mod_' timePeriod];

allMods = fields(os(1).(modToUse))';

if ~isempty(p.Results.particularModels)
    allMods = intersect(allMods, p.Results.particularModels, 'stable');
end

switch p.Results.scoreToUse
    % for AIC or BIC, grab the pre-computed score
    case {'AIC','BIC'}
        for i = 1:length(os)
            tmp_scores = [];
            for curr_mod = allMods
                curr_mod = curr_mod{:};
                tmp_scores = [tmp_scores; os(i).(modToUse).(curr_mod).(p.Results.scoreToUse)];
            end
            bestModStruct.all_scores = [bestModStruct.all_scores tmp_scores];
        end

        [~, bestMod] = min(bestModStruct.all_scores);
    case 'LR'
        warning('off','econ:lratiotest:RLLExceedsULL') % if LH values are very close, you'll get a ton of warnings; turn off
        
        % for LR test, a touch more involved; start with unrestricted model and go down stepwise
        % first find the dof for each model and sort in descending order
        LR_dof = [];
        for curr_mod = allMods
            curr_mod = curr_mod{:};
            LR_dof = [LR_dof length(os(1).(modToUse).(curr_mod).bestParams)];
        end
        [LR_dof, LR_order] = sort(LR_dof,'descend');
        allMods = allMods(LR_order); % reorder the models if necessary

        bestMod = [];
        for i = 1:length(os)
            modFound_flag = false;
            modInd = 1;
            while modFound_flag == false
                h = lratiotest(os(i).(modToUse).(allMods{modInd}).LH, ...
                               os(i).(modToUse).(allMods{modInd + 1}).LH, ...
                               LR_dof(modInd) - LR_dof(modInd + 1), ...
                               p.Results.pValue);
                if h == 1 % if the improvement is significant, keep it and move on
                    modFound_flag = true;
                    bestMod = [bestMod modInd];
                else % if the improvement is insigificant, mov down
                    modInd = modInd + 1; % iterate to the next model
                    if modInd == length(allMods) % if this is the last model, then it's best
                        modFound_flag = true;
                        bestMod = [bestMod modInd];
                    end
                end
            end
        end
        warning('on','econ:lratiotest:RLLExceedsULL') % turn warnings back on
end

bestModStruct.bestMod = bestMod;
bestModStruct.bestMod_name = allMods(bestMod);

for curr_mod = allMods
    curr_mod = curr_mod{:};
    bestModStruct.(['mask_' curr_mod]) = strcmp(bestModStruct.bestMod_name, curr_mod);
end

if p.Results.plotModels_Flag == true
    figure
    bins = 0.5:length(allMods) + 0.5;
    if p.Results.normalizeHistogram_Flag == true
        histogram(bestMod, bins, 'Normalization', 'probability')
    else
        histogram(bestMod, bins)
    end

    xlim([min(bins)-0.5 max(bins)+0.5])
    set(gca,'tickdir','out', 'xtick', 1:length(allMods) + 0.5,...
        'xticklabel', strrep(allMods, '_', '-'), ...
        'xticklabelrotation',60)
    ylabel('Number of neurons')
    title(['Lowest (best) ' p.Results.scoreToUse])
end