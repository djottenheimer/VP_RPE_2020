function hOut = plotRegressionWithCI(model, coefInds, figHandle, varargin)

p = inputParser;
p.addParameter('LineWidth', 2)
p.addParameter('LineStyle', '-')
p.addParameter('XOffset', 0)
p.addParameter('Color', [0 0 0]);
p.parse(varargin{:});

if nargin < 3
    figure
else
    subplot(figHandle)
end

coef = model.Coefficients.Estimate(coefInds);
coefCI = model.coefCI;
coefCI = coefCI(coefInds, 2) - coef;
hOut = errorbar((1:length(coefInds)) + p.Results.XOffset, coef, coefCI, p.Results.LineStyle, 'Color', p.Results.Color, 'linewidth', p.Results.LineWidth);