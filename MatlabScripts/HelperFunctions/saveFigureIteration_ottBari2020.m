function saveFigureIteration_ottBari2020(fHandle, saveFigLoc, figureName, varargin)

answer = questdlg('Save figure?');
switch answer
    case 'Yes'
        p = inputParser;
        p.addParameter('figureSubName', '')
        p.addParameter('FigureSize','')
        p.parse(varargin{:});

        if strcmp(p.Results.FigureSize, 'max')
        %     set(fHandle, 'Position', get(0,'Screensize'))
            set(fHandle, 'units', 'normalized', 'outerposition', [0 0 1 1]);
        elseif strcmp(p.Results.FigureSize, 'half')
            set(fHandle, 'units','normalized','outerposition',[0 0 0.5 1]);
        end
        fHandle.Renderer = 'Painters';

        tmp = dir(saveFigLoc);
        nextVersion = sum(contains({tmp.name}, '.pdf') & contains({tmp.name}, [figureName '_' p.Results.figureSubName]));
        if nextVersion >= 10
            version = int2str(nextVersion);
        else
            version = ['0' int2str(nextVersion)];
        end
        if isempty(p.Results.figureSubName)
            saveName = fullfile(saveFigLoc, ['ottBari2019_' figureName '_v' version]);
        else
            saveName = fullfile(saveFigLoc, ['ottBari2019_' figureName '_' p.Results.figureSubName '_v' version]);
        end

        saveFigurePDF(fHandle, [saveName '.pdf']);
        savefig(fHandle, saveName);

        fprintf('\nSaved as %s\n', saveName);
    otherwise
        fprintf('Figure not saved\n')
end