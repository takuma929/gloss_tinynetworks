%--------------------------------------------------------------------------
% This script generates Figure 6 for the manuscript.
% It visualizes model predictions for highlight and surface manipulations by:
% - Loading human gloss labels and predicted responses from one-layer and 
%   three-layer neural network models under various manipulations,
% - Normalizing model predictions relative to the baseline condition,
% - Plotting the effect of highlight rotation, highlight translation, surface roughness,
%   and surface contrast manipulations,
% - Drawing error bars and overlays for both models, with consistent color and layout,
% - Saving each panel as a vector graphic in the 'figs' directory.
%
% Author: TM, 2025
%--------------------------------------------------------------------------
clearvars; close all;

disp('Generating figure 6...')

%% Draw a figure to show manipulated highlights
load(fullfile('data','fig_parameters'))

imgN = 3888;
label_human = table2array(readtable(fullfile('data','humanlabel.csv')));
label_gt    = table2array(readtable(fullfile('data','groundtruthlabel.csv')));

lw = 0.2;
symbolsize = 30;

% Load predictions
pred.oneLayer    = load(fullfile('data','networks','human_threelayer_kernelN64_trainedby3888imgs','responses_highlight_manipulation.mat'));
pred.threeLayers = load(fullfile('data','networks','human_onelayer_kernelN1_trainedby3888imgs','responses_highlight_manipulation.mat'));

% Normalize responses
norm_by_col = @(X, col) (X-X(:,col))./repmat(X(:,col),1,size(X,2))*100;

pred.oneLayer.pred_highlightRotated_normalized    = norm_by_col(pred.oneLayer.pred_highlightRotated,4);
pred.threeLayers.pred_highlightRotated_normalized = norm_by_col(pred.threeLayers.pred_highlightRotated,4);

pred.oneLayer.pred_highlightTranslated_normalized    = norm_by_col(pred.oneLayer.pred_highlightTranslated,4);
pred.threeLayers.pred_highlightTranslated_normalized = norm_by_col(pred.threeLayers.pred_highlightTranslated,4);

pred.oneLayer.pred_roughness_normalized    = (pred.oneLayer.pred_roughness - pred.oneLayer.pred_roughness(:,2)) ./ repmat(pred.oneLayer.pred_roughness(:,1),1,5)*100;
pred.threeLayers.pred_roughness_normalized = (pred.threeLayers.pred_roughness - pred.threeLayers.pred_roughness(:,2)) ./ repmat(pred.threeLayers.pred_roughness(:,1),1,5)*100;

pred.oneLayer.pred_contrast_normalized    = norm_by_col(pred.oneLayer.pred_contrast,5);
pred.threeLayers.pred_contrast_normalized = norm_by_col(pred.threeLayers.pred_contrast,5);

% Plot parameters for each manipulation
params = struct( ...
    'pred_highlightRotated_normalized',    {struct('xlim',[0.5,7.5],'xlabel','Rotation angle [deg]', ...
        'xticks',1:7, 'xticklabels',{{'-90','-60','-30','0','+30','+60','+90'}}, ...
        'yticks',-14:7:0, 'ylim',[-14 1], 'yticklabels',{{'-14','-7','0'}} )}, ...
    'pred_highlightTranslated_normalized', {struct('xlim',[0.5,7.5],'xlabel','Translation [pixel]', ...
        'xticks',1:1.5:7, 'xticklabels',{{'-30','-15','0','+15','+30'}}, ...
        'yticks',-8:4:0, 'ylim',[-8 1], 'yticklabels',{{'-8','-4','0'}} )}, ...
    'pred_roughness_normalized',           {struct('xlim',[0.7,5.3],'xlabel','Surface roughness [-]', ...
        'xticks',1:5, 'xticklabels',{{'ε','0.05','0.10','0.15','0.20'}}, ...
        'yticks',-12:6:6,'ylim',[-12 10],'yticklabels',{{'-12','-6','0','6'}} )}, ...
    'pred_contrast_normalized',            {struct('xlim',[0.7,5.3],'xlabel','Surface contrast [%]', ...
        'xticks',1:5, 'xticklabels',{{'0','25','50','75','100'}}, ...
        'yticks',-40:20:0,'ylim',[-40 4],'yticklabels',{{'-40','20','0'}} )} ...
);

colors = struct('threeLayers',[78 85 246]/255, 'oneLayer',[68 84 106]/255);

for manipulation = fieldnames(params)'
    man = manipulation{1};
    fig = figure; hold on

    for model = {'threeLayers','oneLayer'}
        m = model{1};
        data = pred.(m).(man);
        y = mean(data);
        y_sd = std(data)/sqrt(size(data,1));
        errorbar(1:length(y), y, y_sd, 'LineStyle','none','Color','k','capsize',0);
        marker = 'o';
        if strcmp(m,'threeLayers'), marker = 'd'; end
        plot(y, [marker,'-'], 'MarkerFaceColor','w', 'MarkerSize',7, ...
            'Color', colors.(m), 'LineWidth',0.6);
    end

    axis square;
    set(gca, ...
        'XLim',params.(man).xlim, 'YLim',params.(man).ylim, ...
        'XTick',params.(man).xticks, 'XTickLabel',params.(man).xticklabels, ...
        'YTick',params.(man).yticks, 'YTickLabel',params.(man).yticklabels, ...
        'FontName','Arial', 'FontSize',figp.fontsize, 'Color',[.97 .97 .97], ...
        'XColor','k', 'YColor','k', 'LineWidth',0.5, ...
        'Units','centimeters', 'Position',[0.9 0.7 3.5 3.5]);
    xlabel(params.(man).xlabel);
    ylabel('Change in gloss level [%]');
    grid on; box off;
    fig.Units = 'centimeters';
    fig.Position = [10,10,figp.twocolumn/4,figp.twocolumn/4];
    fig.Color = 'w'; fig.InvertHardcopy = 'off';

    ticklengthcm(gca,0.0)
    exportgraphics(fig, fullfile('figs', ['fig6(highlightManipulation_',man,').pdf']), 'ContentType','vector')
end

disp('Done.')
close all