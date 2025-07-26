%--------------------------------------------------------------------------
% This script generates Figure 3 for the manuscript.
% It compares computational model performance in gloss perception by:
% - Loading human and model response data,
% - Computing and comparing correlation coefficients for one-layer, three-layer,
%   and ResNet18-based networks (with/without additional training data),
% - Including image statistics and specular metrics models,
% - Visualizing the relationship between correlation-to-ground-truth and
%   correlation-to-human for all approaches in a single comprehensive plot.
% The script saves the resulting figure in the 'figs' directory.
%
% Author: TM, 2025
%--------------------------------------------------------------------------

% Source code to generate figures
clearvars; close all; % cleaning

disp('Generating figure 3...')

%% Load Data
load(fullfile('data','onlineData'))
load(fullfile('data','imageStats_corrCoeff'))
load(fullfile('data','imgStats_multiRegression_corrCoeff'))
load(fullfile('data','fig_parameters'))

%% One-layer and three-la models: load correlation coefficient
kernelN_list.onelayer = [1 2 4 9];
kernelN_list.threelayer = [9 16 32 64];

filedir.onelayer = fullfile('data','networks','onelayer_models');
filedir.threelayer = fullfile('data','networks','threelayer_models');

cvtypeList = {'shape','lighting'};
traininglabelList = {'human','groundtruth'};

for model = {'onelayer','threelayer'}
    for kernelN_idx = 1:length(kernelN_list.(model{1}))
        kernelN = kernelN_list.(model{1})(kernelN_idx);
        for traininglabel = traininglabelList
            for cvtype = cvtypeList % validation type (shape-based or lighting-based cross validation)
                for ii = 1:12

                    % Determine label ordering for correlation
                    if strcmp(traininglabel{1}, 'human')
                        % correlation to human response (human), or correlation to
                        % physical ground-truth (gt)
                        corrlabelList = {'human', 'gt'};
                    else
                        corrlabelList = {'gt', 'human'};
                    end

                    for corr_label = corrlabelList
                        fname = fullfile(filedir.(model{1}), [traininglabel{1},'_kernelN',num2str(kernelN), ...
                            '_',cvtype{1}, num2str(ii),'/corrs_',corr_label{1},'.csv']);

                        temp = readmatrix(fname);

                        % Pick the max correlation row according to label type
                        % and use the same id to get the correlation for the
                        % other objective
                        if (strcmp(traininglabel{1},'human') && strcmp(corr_label{1},'human')) || ...
                           (strcmp(traininglabel{1},'groundtruth') && strcmp(corr_label{1},'gt'))
                            [~,maxid] = max(temp(:,1));
                        end

                        corr_all.(model{1}).(traininglabel{1}).(corr_label{1})(ii,strcmp(cvtype{1},cvtypeList),kernelN_idx) = temp(maxid,1);
                    end
                end
            end
        end
    end
end

% Reshape and average correlation results for plotting
for traininglabel = {'human','groundtruth'}
    for corr_label = {'human','gt'}
        oneLayer_corr.(traininglabel{1}).(corr_label{1}) = mean(reshape(corr_all.onelayer.(traininglabel{1}).(corr_label{1}),24,length(kernelN_list.onelayer)));
        threeLayers_corr.(traininglabel{1}).(corr_label{1}) = mean(reshape(corr_all.threelayer.(traininglabel{1}).(corr_label{1}),24,length(kernelN_list.threelayer)));
    end
end

%% Visualization & Figure generation
cnt = 0;lw = 0.2;symbolsize = 50;

% set color codes
c_human = [90 152 152]/255;
c_human_mean = [160 212 212]/255;
c_oneLayer_human = [78 85 246]/255;
c_oneLayer_gt = [184 0 127]/255;
c_specularStr = [151 217 92]/255;
c_ResNet18 = c_oneLayer_gt;
c_additionalImage = [248 92 1]/255;

basedir_additional = fullfile('data','networks','additional_trainingimgs');

% load ResNet18 trained on ground-truth (additional training)

% additional 500,000 imgs
temp = readmatrix([basedir_additional,'/groundtruth_ResNet18_kernelN64_addtrainingN_500000/corrs_gt.csv']);
[ResNet18_additional500000.gt,maxid] = max(temp(:,1));
temp = readmatrix(fullfile(basedir_additional,'groundtruth_ResNet18_kernelN64_addtrainingN_500000/corrs_human.csv'));
ResNet18_additional500000.human = temp(maxid,1);

% additional 100,000 imgs
temp = readmatrix([basedir_additional,'/groundtruth_ResNet18_kernelN64_addtrainingN_100000/corrs_gt.csv']);
[ResNet18_additional100000.gt,maxid] = max(temp(:,1));
temp = readmatrix([basedir_additional,'/groundtruth_ResNet18_kernelN64_addtrainingN_100000/corrs_human.csv']);
ResNet18_additional100000.human = temp(maxid,1);

% additional 10,000 imgs
temp = readmatrix([basedir_additional,'/groundtruth_ResNet18_kernelN64_addtrainingN_10000/corrs_human.csv']);
[ResNet18_additional10000.human,maxid] = max(temp(:,1)); 
temp = readmatrix([basedir_additional,'/groundtruth_ResNet18_kernelN64_addtrainingN_10000/corrs_gt.csv']);
ResNet18_additional10000.gt = temp(maxid,1);

% load three-layer models trained on ground-truth (additional training)

% additional 500,000 imgs
temp = readmatrix([basedir_additional,'/groundtruth_threelayer_kernelN64_addtrainingN_500000/corrs_gt.csv']);
[twoArea_additional500000.gt,maxid] = max(temp(:,1)); 
temp = readmatrix([basedir_additional,'/groundtruth_threelayer_kernelN64_addtrainingN_500000/corrs_human.csv']);
twoArea_additional500000.human =temp(maxid,1); 

% additional 100,000 imgs
temp = readmatrix([basedir_additional,'/groundtruth_threelayer_kernelN64_addtrainingN_100000/corrs_gt.csv']);
[twoArea_additional100000.gt,maxid] = max(temp(:,1)); 
temp = readmatrix([basedir_additional,'/groundtruth_threelayer_kernelN64_addtrainingN_100000/corrs_human.csv']);
twoArea_additional100000.human = temp(maxid,1); 

% additional 10,000 imgs
temp = readmatrix([basedir_additional,'/groundtruth_threelayer_kernelN64_addtrainingN_10000/corrs_gt.csv']);
[twoArea_additional10000.gt,maxid] = max(temp(:,1)); 
temp = readmatrix([basedir_additional,'/groundtruth_threelayer_kernelN64_addtrainingN_10000/corrs_human.csv']);
twoArea_additional10000.human = temp(maxid,1); 

% load ResNet18 trained on ground-truth (no additional training)
basedir_ResNet18 = fullfile('data','networks');
temp = readmatrix([basedir_ResNet18,'/groundtruth_ResNet18_kernelN64/corrs_gt.csv']);
[ResNet18.gt,maxid] = max(temp(:,1)); 
temp = readmatrix([basedir_ResNet18,'/groundtruth_ResNet18_kernelN64/corrs_human.csv']);
ResNet18.human = temp(maxid,1); 

%% compute correlation to groundtruth and correlation to other participants for each participant
groupN_summary = zeros(54,1);
for N = 1:length(data)
    onlineData(:,:,N) = data(N).response_Pellacini_c;
    groupN_summary(N) = data(N).groupN;
end

for groupN = 1:54
    idx = find(groupN_summary == groupN);
    for N = 1:length(idx)
        response_group.(['group',num2str(groupN)])(:,N) = mean(data(idx(N)).response_Pellacini_c,2);
        if N == 1
            gt_group.(['group',num2str(groupN)])(:,1) = gt(idx(N)).Pellacini_c(1:84);
        end
    end
end

for groupN = 1:54
    for subjectN = 1:size(response_group.(['group',num2str(groupN)]),2)
        cnt = cnt + 1;
        obs1 = response_group.(['group',num2str(groupN)])(:,subjectN);
        obs_rest = mean(response_group.(['group',num2str(groupN)])(:,[1:subjectN-1,subjectN+1:end]),2);
        gt_temp = gt_group.(['group',num2str(groupN)]);

        human_human_corr(cnt) = corr(obs1,obs_rest);
        human_gt_corr(cnt) = corr(obs1,gt_temp);
    end
end

%% Generate figure
fig = figure;
ax = gca;

% fill the figure panel with different colors and draw a diagonal line
fill([0,0,1],[0,1,1],c_oneLayer_human,'FaceAlpha',0.05,'EdgeColor','none');hold on;
fill([0,1,1],[0,0,1],c_oneLayer_gt,'FaceAlpha',0.05,'EdgeColor','none');hold on;
line([0,100],[0,100],'Color','k','LineWidth',0.5)

%%%%%% plot human participant %%%%%%
scatter(human_gt_corr,human_human_corr,20,c_human,'o','filled','MarkerEdgeColor','none','LineWidth',lw,'MarkerFaceAlpha',0.4);hold on;
scatter(mean(human_gt_corr),mean(human_human_corr),50,c_human_mean,'x','filled','MarkerEdgeColor',[0 0 0],'MarkerFaceAlpha',1,'LineWidth',1.5);hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% plot one-Layer rgb (human label) %%%%%%
for kernelN = 1:size(oneLayer_corr.human.human,2)
    scatter(mean(oneLayer_corr.human.gt(:,kernelN)),mean(oneLayer_corr.human.human(:,kernelN)),40,c_oneLayer_human,'^','filled','MarkerEdgeColor',[0 0 0],'MarkerFaceAlpha',0.5,'LineWidth',lw);hold on;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% plot one-Layer rgb (gt label) %%%%%%
for kernelN = 1:size(oneLayer_corr.groundtruth.human,2)
    scatter(mean(oneLayer_corr.groundtruth.gt(:,kernelN)),mean(oneLayer_corr.groundtruth.human(:,kernelN)),40,c_oneLayer_gt,'^','filled','MarkerEdgeColor',[0 0 0],'MarkerFaceAlpha',0.5,'LineWidth',lw);hold on;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% plot two-Layers rgb (human label) %%%%%%
for kernelN = 1:size(threeLayers_corr.human.human,2)
    scatter(mean(threeLayers_corr.human.gt(:,kernelN)),mean(threeLayers_corr.human.human(:,kernelN)),50,c_oneLayer_human,'s','filled','MarkerEdgeColor',[0 0 0],'MarkerFaceAlpha',0.5,'LineWidth',lw);hold on;
    scatter(mean(threeLayers_corr.groundtruth.gt(:,kernelN)),mean(threeLayers_corr.groundtruth.human(:,kernelN)),50,c_oneLayer_gt,'s','filled','MarkerEdgeColor',[0 0 0],'MarkerFaceAlpha',0.5,'LineWidth',lw);hold on;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% plot ResNet18 (no additional) %%%%%%
scatter(ResNet18.gt,ResNet18.human,40,c_ResNet18,'>','filled','MarkerEdgeColor',[0 0 0],'MarkerFaceAlpha',0.5,'LineWidth',lw);hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% plot networks trained on additional images %%%%%%
scatter(ResNet18_additional100000.gt,ResNet18_additional100000.human,40,c_additionalImage,'>','filled','MarkerEdgeColor',[0 0 0],'MarkerFaceAlpha',0.5,'LineWidth',lw);hold on;
scatter(ResNet18_additional500000.gt,ResNet18_additional500000.human,40,c_additionalImage,'>','filled','MarkerEdgeColor',[0 0 0],'MarkerFaceAlpha',0.5,'LineWidth',lw);hold on;

hold on;

%scatter(twoArea_additional10000.gt,twoArea_additional10000.human,50,c_oneLayer_human,'s','filled','MarkerEdgeColor',[0 0 0],'MarkerFaceAlpha',0.5,'LineWidth',lw);hold on;
scatter(twoArea_additional100000.gt,twoArea_additional100000.human,50,c_additionalImage,'s','filled','MarkerEdgeColor',[0 0 0],'MarkerFaceAlpha',0.5,'LineWidth',lw);hold on;
scatter(twoArea_additional500000.gt,twoArea_additional500000.human,50,c_additionalImage,'s','filled','MarkerEdgeColor',[0 0 0],'MarkerFaceAlpha',0.5,'LineWidth',lw);hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% plot image statistics models %%%%%%
for N = 1:length(corrCoeff_imgStats.label)
    x = corrCoeff_imgStats.gtvsmodel(N,:);
    y = corrCoeff_imgStats.humanvsmodel(N,:);
    scatter(mean(x),mean(y),40,[.8 .8 .8],'d','filled','MarkerEdgeColor',[0 0 0],'MarkerFaceAlpha',1,'LineWidth',lw);hold on;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% plot image statistics Multiple regression models %%%%%%
x = corrCoeff_imgStats_multiRegression.gtvsmodel;
y = corrCoeff_imgStats_multiRegression.humanvsmodel;
scatter(mean(x),mean(y),60,[.8 .8 .8],'p','filled','MarkerEdgeColor',[0 0 0],'MarkerFaceAlpha',1,'LineWidth',lw);hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% plot image statistics Multiple regression models %%%%%%
load(fullfile('data','specularMetrics_multiRegression_corrCoeff.mat'))
x = corrCoeff_specularMetrics_multiRegression.gtvsmodel;
y = corrCoeff_specularMetrics_multiRegression.humanvsmodel;
scatter(mean(x),mean(y),60,c_specularStr,'p','filled','MarkerEdgeColor',[0 0 0],'MarkerFaceAlpha',0.5,'LineWidth',lw);hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% plot specular reflection models %%%%%%
for modellabel = {'contrast','coverage','sharpness'}
    x = abs(corrCoeff_specularMetrics.gt.(modellabel{1}));
    y = abs(corrCoeff_specularMetrics.human.(modellabel{1}));

    scatter(mean(x),mean(y),40,c_specularStr,'v','filled','MarkerEdgeColor',[0 0 0],'MarkerFaceAlpha',0.5,'LineWidth',lw);hold on;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xlim([0 1]);ylim([0 1]);axis square
xlabel('Correlation to ground-truth','FontWeight', 'Bold');ylabel('Correlation to human','FontWeight', 'Bold');

fig.Units           = 'centimeters';
fig.Position = [10,10,figp.twocolumn/2,figp.twocolumn/2];
fig.Color           = 'w';
fig.InvertHardcopy  = 'off';

xticks(0:0.25:1)
yticks(0:0.25:1)

ax.XTickLabel = {'0.00','0.25','0.50','0.75','1.00'};
ax.YTickLabel = {'0.00','0.25','0.50','0.75','1.00'};

ax.FontName = 'Arial';
ax.Color = ones(1,3);
ax.FontSize = figp.fontsize;
ax.XColor = 'k';ax.YColor = 'k';

ax.LineWidth = 0.5;
ax.Units = 'centimeters';
ax.Position = [0.95 0.85 7.6 7.6];
ticklengthcm(ax,0.0)
grid minor
box off
exportgraphics(fig,fullfile('figs','fig3(model_comparison).pdf'),'ContentType','vector')

disp('Done.')
close all
