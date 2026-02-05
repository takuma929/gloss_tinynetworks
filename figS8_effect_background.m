%--------------------------------------------------------------------------
% This script generates Supplementary Figure S5 for the manuscript.
% It compares model performance with and without background cues by:
% - Loading correlation coefficients for networks trained with/without image backgrounds,
% - Plotting paired observer results and means for each network type (one-layer and three-layer),
% - Formatting figures showing correlation to human judgments as a function of kernel number,
%
% Author: TM, 2025
%--------------------------------------------------------------------------

%% Clean workspace and load data
clearvars; close all;

disp('Generating figure S8...')

% Load figure parameters and correlation data
load(fullfile('data', 'fig_parameters'))
load(fullfile('data', 'corr_bg_nobg'), 'corr_nobg', 'corr_bg')

% Visualization parameters
gap  = 0.2;
navy = [0.2855, 0.3294, 0.3890];

%% Loop through network types and plot
for area = {'onelayer', 'threelayer'}
    % Extract correlation data
    bg_corr   = corr_bg.(area{1}).human;
    nobg_corr = corr_nobg.(area{1}).human;

    % Create figure
    fig = figure;
    hold on

    % Plot paired lines and scatter means
    for N = 1:size(bg_corr, 2)
        for n = 1:24
            line([N-gap, N+gap], [nobg_corr(n, N), bg_corr(n, N)], ...
                'LineWidth', 0.3, 'Color', [0.7 0.7 0.7]);
        end
        scatter(N-gap, mean(nobg_corr(:, N)), 25, navy, ...
            'o', 'filled', 'MarkerEdgeColor', 'none', 'LineWidth', 0.2);
        scatter(N+gap, mean(bg_corr(:, N)), 25, navy, ...
            'd', 'filled', 'MarkerEdgeColor', 'none', 'LineWidth', 0.2);
    end

    % Axes formatting
    ax = gca;
    ax.Units       = 'centimeters';
    ax.Position    = [1.0 0.9 2.4 3.2];
    ax.FontName    = 'Arial';
    ax.FontSize    = figp.fontsize;
    ax.Color       = [0.97 0.97 0.97];
    ax.XColor      = 'k';
    ax.YColor      = 'k';
    ax.TickDir     = 'out';
    ax.LineWidth   = 0.5;
    ax.XTick       = 1:4;
    ax.YTick       = 0.5:0.25:1;
    ax.YTickLabel  = {'0.50','0.75','1.00'};
    xlim([0.4 4.6]);
    ylim([0.5 1]);

    % Custom X-axis labels depending on area type
    if strcmp(area{1}, 'onelayer')
        ax.XTickLabel = {'1','2','4','9'};
    elseif strcmp(area{1}, 'threelayer')
        ax.XTickLabel = {'9','16','32','64'};
    end

    % Labels
    xlabel('kernel N', 'FontWeight', 'Bold');
    ylabel('Correlation to human', 'FontWeight', 'Bold');

    % Figure formatting
    fig.Units          = 'centimeters';
    fig.Position       = [10, 10, figp.twocolumn/5, figp.twocolumn/4];
    fig.Color          = 'w';
    fig.InvertHardcopy = 'off';

    ticklengthcm(ax, 0.1);
    grid off
    box off

    % Export to PDF
    exportgraphics(fig, fullfile('figs', ['figS8(bgvsnobg_', area{1}, ').pdf']), ...
        'ContentType', 'vector');
end


%% Statistical comparison: background vs. no-background
% Report:
% - two-sided paired t-test: t(df)=..., p=...
% - Cohen's dz for paired differences: dz = mean(diff)/std(diff)
% - 95% CI for Cohen's dz (via noncentral t inversion)
%
% Notes:
% - CI reported below is for effect size dz (NOT for mean difference).
% - If you also want CI of mean diff, keep ci_* from ttest as well.

alpha = 0.05;   % 95% CI

% -----------------------------
% Three-layer
% -----------------------------
bg_vals_two   = mean(corr_bg.threelayer.human, 2);
nobg_vals_two = mean(corr_nobg.threelayer.human, 2);
diff_two = bg_vals_two - nobg_vals_two;

% Paired t-test
[~, p_two, ~, stats_two] = ttest(diff_two, 0, 'Alpha', alpha, 'Tail', 'both');

n_two  = numel(diff_two);
df_two = stats_two.df;
t_two  = stats_two.tstat;

% Cohen's dz
dz_two = mean(diff_two) / std(diff_two);   % same as t/sqrt(n) up to rounding

% 95% CI for Cohen's dz via noncentral t inversion
[deltaL_two, deltaU_two] = local_nct_delta_ci(t_two, df_two, alpha);
dzL_two = deltaL_two / sqrt(n_two);
dzU_two = deltaU_two / sqrt(n_two);

% -----------------------------
% One-layer
% -----------------------------
bg_vals_one   = mean(corr_bg.onelayer.human, 2);
nobg_vals_one = mean(corr_nobg.onelayer.human, 2);
diff_one = bg_vals_one - nobg_vals_one;

% Paired t-test
[~, p_one, ~, stats_one] = ttest(diff_one, 0, 'Alpha', alpha, 'Tail', 'both');

n_one  = numel(diff_one);
df_one = stats_one.df;
t_one  = stats_one.tstat;

% Cohen's dz
dz_one = mean(diff_one) / std(diff_one);

% 95% CI for Cohen's dz via noncentral t inversion
[deltaL_one, deltaU_one] = local_nct_delta_ci(t_one, df_one, alpha);
dzL_one = deltaL_one / sqrt(n_one);
dzU_one = deltaU_one / sqrt(n_one);

% -----------------------------
% Display results (Nature-ish formatting)
% -----------------------------
fprintf('--- Three-layer stats ---\n');
fprintf('t(%d) = %.2f, p = %.3f, Cohen''s dz = %.3f, 95%% CI [%.3f, %.3f]\n', ...
    df_two, t_two, p_two, dz_two, dzL_two, dzU_two);

fprintf('--- One-layer stats ---\n');
fprintf('t(%d) = %.2f, p = %.3f, Cohen''s dz = %.3f, 95%% CI [%.3f, %.3f]\n', ...
    df_one, t_one, p_one, dz_one, dzL_one, dzU_one);

fprintf('Done.\n');
close all


% ========================================================================
% Local function (must be at end of script file in MATLAB)
% ========================================================================
function [deltaL, deltaU] = local_nct_delta_ci(tObs, df, alpha)
% Returns 100*(1-alpha)% CI for the noncentrality parameter (delta)
% of a noncentral t distribution, inverted from the observed t statistic.
%
% We solve:
%   nctcdf(tObs; df, deltaL) = 1 - alpha/2
%   nctcdf(tObs; df, deltaU) = alpha/2
%
% This matches the usual inversion for a two-sided CI.

    if exist('nctcdf','file') ~= 2
        error('nctcdf not found. Requires Statistics and Machine Learning Toolbox.');
    end

    targetL = 1 - alpha/2;  % 0.975
    targetU = alpha/2;      % 0.025

    % Monotonicity: nctcdf(tObs, df, delta) decreases as delta increases.
    % So deltaL < deltaU typically, but for safety we'll bracket separately.

    % Solve for deltaL where CDF = 0.975 (smaller delta)
    deltaL = local_solve_delta(tObs, df, targetL);

    % Solve for deltaU where CDF = 0.025 (larger delta)
    deltaU = local_solve_delta(tObs, df, targetU);

    % Ensure ordering
    if deltaL > deltaU
        tmp = deltaL; deltaL = deltaU; deltaU = tmp;
    end
end

function delta = local_solve_delta(tObs, df, target)
% Solve nctcdf(tObs; df, delta) = target for delta, using bracketing + fzero.

    fun = @(d) nctcdf(tObs, df, d) - target;

    % Start bracket around observed delta ~ tObs (roughly)
    lo = -1e4;
    hi =  1e4;

    flo = fun(lo);
    fhi = fun(hi);

    % Expand bracket if needed (rare, but safe)
    iter = 0;
    while sign(flo) == sign(fhi) && iter < 50
        lo = lo * 2;
        hi = hi * 2;
        flo = fun(lo);
        fhi = fun(hi);
        iter = iter + 1;
    end

    if sign(flo) == sign(fhi)
        % Best-effort fallback if bracketing fails
        delta = tObs;
        return;
    end

    delta = fzero(fun, [lo, hi]);
end
