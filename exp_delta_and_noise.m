%% exp_delta_and_noise.m  (final version for paper figures)
% Extra experiments for QET:
%  (i) Effect of quantization step size Delta on performance
% (ii) Sensitivity w.r.t. measurement noise variance R
%
% 输出：一张 2x2 的综合图：
% (a) RMSE vs Delta
% (b) Avg. comm vs Delta
% (c) RMSE vs R
% (d) Avg. comm vs R

clear; clc; close all;

%% ---------- Plot style ----------
set(groot, 'defaultAxesFontSize', 11, ...
           'defaultAxesLabelFontSizeMultiplier', 1.0, ...
           'defaultLineLineWidth', 1.3, ...
           'defaultTextInterpreter', 'latex', ...
           'defaultAxesTickLabelInterpreter', 'latex', ...
           'defaultLegendInterpreter', 'latex');

rng(1);   % 便于复现

%% ---------- System (same as main_experiments) ----------
T  = 0.1;

A  = [0.95 T    0    0;
      0    0.95 0    0;
      0    0    0.95 T;
      0    0    0    0.95];

Q  = 1e-3 * eye(4);
x0_true = [0; 1; 0; 0.5];

C{1} = [1 0 0 0];
C{2} = [0 0 1 0];
C{3} = [1 0 1 0] / sqrt(2);
m     = numel(C);

K  = 300;
MC = 50;

gamma_QET = 1.0;   % same as main_experiments
L_levels  = 2001;  % same as main_experiments

%% ========================================================================
%% Experiment 4: Delta sweep
%% ========================================================================
Delta_list = [0.05 0.10 0.20 0.30 0.40];

rmse_ss_Delta  = zeros(size(Delta_list));
comm_avg_Delta = zeros(size(Delta_list));

R0 = 0.1;   % baseline measurement noise variance

for id = 1:numel(Delta_list)
    Delta = Delta_list(id);

    % build R cell for current case
    R = cell(1,m);
    for i = 1:m
        R{i} = R0;
    end

    [rmse_ss_Delta(id), comm_avg_Delta(id)] = ...
        sweep_threshold(A, Q, C, R, x0_true, T, ...
                        gamma_QET, Delta, L_levels, K, MC);

    fprintf('[Delta sweep] %d/%d: Delta = %.3f, RMSE = %.4f, Comm = %.4f\n', ...
        id, numel(Delta_list), Delta, rmse_ss_Delta(id), comm_avg_Delta(id));
end

%% ========================================================================
%% Experiment 5: R sweep
%% ========================================================================
R_list = [0.05 0.10 0.20 0.40 0.80];
Delta_fixed = 0.20;   % use a representative quantization step size

rmse_ss_R  = zeros(size(R_list));
comm_avg_R = zeros(size(R_list));

for ir = 1:numel(R_list)
    Rval = R_list(ir);

    R = cell(1,m);
    for i = 1:m
        R{i} = Rval;
    end

    [rmse_ss_R(ir), comm_avg_R(ir)] = ...
        sweep_threshold(A, Q, C, R, x0_true, T, ...
                        gamma_QET, Delta_fixed, L_levels, K, MC);

    fprintf('[R sweep] %d/%d: R = %.3f, RMSE = %.4f, Comm = %.4f\n', ...
        ir, numel(R_list), Rval, rmse_ss_R(ir), comm_avg_R(ir));
end

%% ========================================================================
%% Combined 2x2 figure
%% ========================================================================
figure;
set(gcf, 'Units','centimeters', 'Position',[2 2 18 14]);  % 图整体尺寸
tiledlayout(2,2, 'Padding','compact', 'TileSpacing','compact');

%% ------- (a) RMSE vs Delta -------
nexttile;
plot(Delta_list, rmse_ss_Delta, '-o', 'MarkerFaceColor','w');
ylabel('Steady-state RMSE');
grid on; box on;
text(0.03, 0.92, '(a)  RMSE vs. $\Delta$', ...
     'Units','normalized','FontWeight','bold');
yl = ylim; ylim([0.95*yl(1), 1.05*yl(2)]);

%% ------- (b) Comm vs Delta -------
nexttile;
plot(Delta_list, comm_avg_Delta, '-s', 'MarkerFaceColor','w');
ylabel('Avg. comm. per sensor');
grid on; box on;
text(0.03, 0.92, '(b)  Comm vs. $\Delta$', ...
     'Units','normalized','FontWeight','bold');
yl = ylim; ylim([0.95*yl(1), 1.05*yl(2)]);

%% ------- (c) RMSE vs R -------
nexttile;
plot(R_list, rmse_ss_R, '-o', 'MarkerFaceColor','w');
xlabel('Measurement noise variance $R$');
ylabel('Steady-state RMSE');
grid on; box on;
text(0.03, 0.92, '(c)  RMSE vs. $R$', ...
     'Units','normalized','FontWeight','bold');
yl = ylim; ylim([0.95*yl(1), 1.05*yl(2)]);

%% ------- (d) Comm vs R -------
nexttile;
plot(R_list, comm_avg_R, '-s', 'MarkerFaceColor','w');
xlabel('Measurement noise variance $R$');
ylabel('Avg. comm. per sensor');
grid on; box on;
text(0.03, 0.92, '(d)  Comm vs. $R$', ...
     'Units','normalized','FontWeight','bold');
yl = ylim; ylim([0.95*yl(1), 1.05*yl(2)]);

% 导出为矢量图
print(gcf, '-depsc2', 'fig_Delta_Rnoise_Combined.eps');
print(gcf, '-dpdf',   'fig_Delta_Rnoise_Combined.pdf');
