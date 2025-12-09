%% main_experiments.m
% Experiments for:
% "Quantized Event-Triggered Distributed and Resilient Filtering 
%  for Bandwidth-Limited Wireless Sensor Networks"

clear; clc; close all;

%% Plot settings (paper quality)
set(groot, 'defaultAxesFontSize', 12, ...
           'defaultAxesLabelFontSizeMultiplier', 1.1, ...
           'defaultLineLineWidth', 1.5, ...
           'defaultTextInterpreter', 'latex', ...
           'defaultAxesTickLabelInterpreter', 'latex', ...
           'defaultLegendInterpreter', 'latex');

rng(1);  % reproducibility

%% ===== System model: 2D nearly constant-velocity, stable =====
T  = 0.1;  % sampling period

% Slightly stable system (0.95 on diagonal) to avoid state blow-up
A  = [0.95 T    0    0;
      0    0.95 0    0;
      0    0    0.95 T;
      0    0    0    0.95];

Q  = 1e-3 * eye(4);   % process noise covariance

x0_true = [0; 1; 0; 0.5];

%% ===== Sensors (3 scalar sensors) =====
C{1} = [1 0 0 0];                    % sensor 1: px
C{2} = [0 0 1 0];                    % sensor 2: py
C{3} = [1 0 1 0] / sqrt(2);          % sensor 3: (px+py)/sqrt(2)

m   = numel(C);
R_i = 0.1;                           % measurement noise variance
for i = 1:m
    R{i} = R_i;
end

%% ===== Simulation parameters =====
K  = 300;          % time steps
MC = 50;           % Monte Carlo runs

% Event-trigger thresholds (scalar NIS ~ chi-square(1))
gamma_ET  = 1.0;
gamma_QET = 1.0;

% Quantization parameters (for measurement y)
% range ≈ ± ((L_levels-1)/2) * Delta ≈ ±200
Delta    = 0.2;        % quantization step size
L_levels = 2001;       % number of levels (odd)
Rq       = (Delta^2)/12;   % quantization noise variance

%% ===== Storage =====
rmse_full   = zeros(K,1);  % Full KF
rmse_ET     = zeros(K,1);  % Event-triggered KF
rmse_QET    = zeros(K,1);  % Quantized ET
rmse_PQKF   = zeros(K,1);  % Periodic Quantized KF (baseline)

comm_full   = zeros(K,1);
comm_ET     = zeros(K,1);
comm_QET    = zeros(K,1);
comm_PQKF   = zeros(K,1);

%% ===== Monte Carlo simulation =====
for mc = 1:MC
    fprintf('MC run %d / %d\n', mc, MC);

    % ----- True trajectory -----
    x_true       = zeros(4,K);
    x_true(:,1)  = x0_true;
    for k = 2:K
        w_k         = mvnrnd(zeros(4,1), Q).';
        x_true(:,k) = A * x_true(:,k-1) + w_k;
    end

    % ----- Measurements -----
    y = cell(m,K);
    for k = 1:K
        for i = 1:m
            v_ik   = sqrt(R{i}) * randn;
            y{i,k} = C{i} * x_true(:,k) + v_ik;
        end
    end

    % ----- Initial estimates (same for all filters) -----
    x_init = x0_true + [0.5; 0; 0.5; 0];  % deliberately biased

    x_hat_full = x_init;
    x_hat_ET   = x_init;
    x_hat_QET  = x_init;
    x_hat_PQKF = x_init;

    P0    = 10 * eye(4);   % large initial covariance
    P_full  = P0;
    P_ET    = P0;
    P_QET   = P0;
    P_PQKF  = P0;

    % local comm counters for this MC run
    comm_full_mc = zeros(K,1);
    comm_ET_mc   = zeros(K,1);
    comm_QET_mc  = zeros(K,1);
    comm_PQKF_mc = zeros(K,1);

    % ----- Time loop -----
    for k = 1:K

        %% ===== Prediction step =====
        if k == 1
            x_pred_full = x_hat_full;  P_pred_full = P_full;
            x_pred_ET   = x_hat_ET;    P_pred_ET   = P_ET;
            x_pred_QET  = x_hat_QET;   P_pred_QET  = P_QET;
            x_pred_PQKF = x_hat_PQKF;  P_pred_PQKF = P_PQKF;
        else
            x_pred_full = A * x_hat_full;
            P_pred_full = A * P_full * A.' + Q;

            x_pred_ET   = A * x_hat_ET;
            P_pred_ET   = A * P_ET * A.' + Q;

            x_pred_QET  = A * x_hat_QET;
            P_pred_QET  = A * P_QET * A.' + Q;

            x_pred_PQKF = A * x_hat_PQKF;
            P_pred_PQKF = A * P_PQKF * A.' + Q;
        end

        %% ===== (1) Full KF with all measurements =====
        x_upd_full = x_pred_full;
        P_upd_full = P_pred_full;
        comm_cnt_full = 0;

        for i = 1:m
            Ci = C{i}; Ri = R{i}; yi = y{i,k};

            e   = yi - Ci * x_upd_full;
            S   = Ci * P_upd_full * Ci.' + Ri;
            K_i = P_upd_full * Ci.' / S;

            x_upd_full = x_upd_full + K_i * e;
            P_upd_full = P_upd_full - K_i * Ci * P_upd_full;

            comm_cnt_full = comm_cnt_full + 1; % one scalar per sensor
        end

        %% ===== (2) Event-triggered, unquantized =====
        x_upd_ET = x_pred_ET;
        P_upd_ET = P_pred_ET;
        comm_cnt_ET = 0;

        for i = 1:m
            Ci = C{i}; Ri = R{i}; yi = y{i,k};

            e   = yi - Ci * x_upd_ET;
            S   = Ci * P_upd_ET * Ci.' + Ri;
            nis = (e.' / S) * e;    % scalar NIS

            if nis > gamma_ET
                K_i      = P_upd_ET * Ci.' / S;
                x_upd_ET = x_upd_ET + K_i * e;
                P_upd_ET = P_upd_ET - K_i * Ci * P_upd_ET;
                comm_cnt_ET = comm_cnt_ET + 1;
            end
        end

        %% ===== (3) Quantized Event-triggered (measurement quantized) =====
        x_upd_QET = x_pred_QET;
        P_upd_QET = P_pred_QET;
        comm_cnt_QET = 0;

        for i = 1:m
            Ci = C{i}; Ri = R{i}; yi = y{i,k};

            % triggering based on unquantized innovation
            e_trig = yi - Ci * x_upd_QET;
            S_trig = Ci * P_upd_QET * Ci.' + Ri;
            nis    = (e_trig.' / S_trig) * e_trig;

            if nis > gamma_QET
                % quantize measurement
                y_q = uniform_quantizer(yi, Delta, L_levels);

                % effective measurement noise = intrinsic + quantization
                R_eff = Ri + Rq;

                e_q = y_q - Ci * x_upd_QET;
                S_q = Ci * P_upd_QET * Ci.' + R_eff;
                K_i = P_upd_QET * Ci.' / S_q;

                x_upd_QET = x_upd_QET + K_i * e_q;
                P_upd_QET = P_upd_QET - K_i * Ci * P_upd_QET;

                comm_cnt_QET = comm_cnt_QET + 1;
            end
        end

        %% ===== (4) Periodic Quantized KF (baseline, every step) =====
        x_upd_PQKF = x_pred_PQKF;
        P_upd_PQKF = P_pred_PQKF;
        comm_cnt_PQKF = 0;

        for i = 1:m
            Ci = C{i}; Ri = R{i}; yi = y{i,k};

            % 每个时刻都发送量化后的测量
            y_q  = uniform_quantizer(yi, Delta, L_levels);
            R_eff = Ri + Rq;

            e_q = y_q - Ci * x_upd_PQKF;
            S_q = Ci * P_upd_PQKF * Ci.' + R_eff;
            K_i = P_upd_PQKF * Ci.' / S_q;

            x_upd_PQKF = x_upd_PQKF + K_i * e_q;
            P_upd_PQKF = P_upd_PQKF - K_i * Ci * P_upd_PQKF;

            comm_cnt_PQKF = comm_cnt_PQKF + 1;
        end

        % Save updated states
        x_hat_full = x_upd_full;   P_full  = P_upd_full;
        x_hat_ET   = x_upd_ET;     P_ET    = P_upd_ET;
        x_hat_QET  = x_upd_QET;    P_QET   = P_upd_QET;
        x_hat_PQKF = x_upd_PQKF;   P_PQKF  = P_upd_PQKF;

        %% ===== RMSE in 2D position =====
        err_full  = x_hat_full([1 3])  - x_true([1 3],k);
        err_ET    = x_hat_ET([1 3])    - x_true([1 3],k);
        err_QET   = x_hat_QET([1 3])   - x_true([1 3],k);
        err_PQKF  = x_hat_PQKF([1 3])  - x_true([1 3],k);

        rmse_full(k) = rmse_full(k) + mean(err_full.^2);
        rmse_ET(k)   = rmse_ET(k)   + mean(err_ET.^2);
        rmse_QET(k)  = rmse_QET(k)  + mean(err_QET.^2);
        rmse_PQKF(k) = rmse_PQKF(k) + mean(err_PQKF.^2);

        %% ===== communication counters =====
        comm_full_mc(k)  = comm_full_mc(k)  + comm_cnt_full;
        comm_ET_mc(k)    = comm_ET_mc(k)    + comm_cnt_ET;
        comm_QET_mc(k)   = comm_QET_mc(k)   + comm_cnt_QET;
        comm_PQKF_mc(k)  = comm_PQKF_mc(k)  + comm_cnt_PQKF;
    end

    % accumulate over MC
    comm_full  = comm_full  + comm_full_mc;
    comm_ET    = comm_ET    + comm_ET_mc;
    comm_QET   = comm_QET   + comm_QET_mc;
    comm_PQKF  = comm_PQKF  + comm_PQKF_mc;
end

%% ===== Average over MC =====
rmse_full  = sqrt(rmse_full  / MC);
rmse_ET    = sqrt(rmse_ET    / MC);
rmse_QET   = sqrt(rmse_QET   / MC);
rmse_PQKF  = sqrt(rmse_PQKF  / MC);

comm_full  = comm_full  / MC;
comm_ET    = comm_ET    / MC;
comm_QET   = comm_QET   / MC;
comm_PQKF  = comm_PQKF  / MC;

avg_comm_full  = comm_full  / m;
avg_comm_ET    = comm_ET    / m;
avg_comm_QET   = comm_QET   / m;
avg_comm_PQKF  = comm_PQKF  / m;

t = (0:K-1) * T;

%% ===== Figure 1: RMSE vs time =====
figure;
plot(t, rmse_full,  'k-',  'DisplayName', 'Full Communication KF'); hold on;
plot(t, rmse_ET,    'b--', 'DisplayName', 'Event-Triggered (ET)');
plot(t, rmse_QET,   'r-.', 'DisplayName', 'Quantized ET');
plot(t, rmse_PQKF,  'm:',  'DisplayName', 'Periodic Quantized KF');
xlabel('Time [s]');
ylabel('Position RMSE');
grid on; box on;
legend('Location','best');
title('Estimation Performance under Different Communication Schemes');
ylim([0.05 0.25]);   % 让差异更明显一点

print(gcf, '-depsc2', 'fig_RMSE_time.eps');
print(gcf, '-dpdf',   'fig_RMSE_time.pdf');

%% ===== Figure 2: communication vs time =====
% 简单移动平均平滑一下通信曲线（只影响显示，不影响统计）
w = 10;
avg_comm_full_s  = movmean(avg_comm_full,  w);
avg_comm_ET_s    = movmean(avg_comm_ET,    w);
avg_comm_QET_s   = movmean(avg_comm_QET,   w);
avg_comm_PQKF_s  = movmean(avg_comm_PQKF,  w);

figure;
plot(t, avg_comm_full_s,  'k-',  'DisplayName', 'Full Communication'); hold on;
plot(t, avg_comm_ET_s,    'b--', 'DisplayName', 'ET');
plot(t, avg_comm_QET_s,   'r-.', 'DisplayName', 'Quantized ET');
plot(t, avg_comm_PQKF_s,  'm:',  'DisplayName', 'Periodic Quantized KF');
xlabel('Time [s]');
ylabel('Avg. transmissions per sensor per step');
grid on; box on;
legend('Location','best');
title('Average Communication Load');

print(gcf, '-depsc2', 'fig_Comm_time.eps');
print(gcf, '-dpdf',   'fig_Comm_time.pdf');

%% ===== Figure 3: Communication–Accuracy Trade-Off (QET only) =====
gamma_list = linspace(0.3, 2.0, 7);   % sweep thresholds
rmse_ss   = zeros(size(gamma_list));
comm_avg  = zeros(size(gamma_list));

for ig = 1:numel(gamma_list)
    gamma_tmp = gamma_list(ig);
    [rmse_ss(ig), comm_avg(ig)] = sweep_threshold(A,Q,C,R,...
        x0_true, T, gamma_tmp, Delta, L_levels, K, MC);
    fprintf('Sweep %d/%d: gamma=%.2f, RMSE=%.4f, Comm=%.4f\n', ...
        ig, numel(gamma_list), gamma_tmp, rmse_ss(ig), comm_avg(ig));
end

figure;
plot(comm_avg, rmse_ss, 'o-');
xlabel('Average transmissions per sensor per step');
ylabel('Steady-state position RMSE');
grid on; box on;
title('Communication--Accuracy Trade-Off (Quantized ET)');

print(gcf, '-depsc2', 'fig_Tradeoff_QET.eps');
print(gcf, '-dpdf',   'fig_Tradeoff_QET.pdf');
