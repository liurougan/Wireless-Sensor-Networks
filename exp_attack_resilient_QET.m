%% exp_attack_resilient_QET.m  (final optimized version for Figure 5)
% Sensor attack / fault experiment:
%   - Baseline: QET without attack
%   - Under attack: QET (fails)
%   - Under attack: Resilient QET (R-QET, robust innovation check)
%
% Figure:
%   RMSE vs. time, with three curves + vertical line "Attack starts"

clear; clc; close all;

%% ---------- Plot style ----------
set(groot, 'defaultAxesFontSize', 12, ...
           'defaultAxesLabelFontSizeMultiplier', 1.1, ...
           'defaultLineLineWidth', 1.5, ...
           'defaultTextInterpreter', 'latex', ...
           'defaultAxesTickLabelInterpreter', 'latex', ...
           'defaultLegendInterpreter', 'latex');

rng(1);   % reproducibility

%% ---------- System (same as main_experiments) ----------
T  = 0.1;  % sampling period

A  = [0.95 T    0    0;
      0    0.95 0    0;
      0    0    0.95 T;
      0    0    0    0.95];

Q  = 1e-3 * eye(4);
x0_true = [0; 1; 0; 0.5];

% 3 scalar sensors
C{1} = [1 0 0 0];                    % sensor 1: px
C{2} = [0 0 1 0];                    % sensor 2: py
C{3} = [1 0 1 0] / sqrt(2);          % sensor 3: (px+py)/sqrt(2)

m   = numel(C);
R_i = 0.1;
for i = 1:m
    R{i} = R_i;
end

%% ---------- QET parameters ----------
Delta    = 0.2;         % same as main_experiments
L_levels = 2001;
Rq       = (Delta^2)/12;

gamma_QET = 1.0;        % NIS-based event trigger threshold
beta_max  = 4.0;        % robust bound for normalized innovation in R-QET

%% ---------- Simulation parameters ----------
K  = 300;               % time steps
MC = 50;                % Monte Carlo runs

% Attack settings
k_attack   = 150;        % attack start time index (1...K)
t_attack   = (k_attack-1) * T;
i_attacked = 1;          % index of attacked sensor
attack_bias = 3.0;       % constant bias added to attacked sensor after k_attack

%% ---------- Storage for RMSE ----------
rmse_QET_no   = zeros(K,1);   % QET without attack
rmse_QET_att  = zeros(K,1);   % QET with attack
rmse_RQET_att = zeros(K,1);   % Resilient QET with attack

%% ---------- Monte Carlo simulation ----------
for mc = 1:MC
    fprintf('MC run %d / %d\n', mc, MC);

    % ---- True state trajectory ----
    x_true      = zeros(4,K);
    x_true(:,1) = x0_true;
    for k = 2:K
        w_k         = mvnrnd(zeros(4,1),Q).';
        x_true(:,k) = A * x_true(:,k-1) + w_k;
    end

    % ---- Measurements: no-attack & attack versions ----
    y_no  = cell(m,K);
    y_att = cell(m,K);

    for k = 1:K
        for i = 1:m
            v_ik = sqrt(R{i}) * randn;
            yi   = C{i} * x_true(:,k) + v_ik;

            % no-attack measurements
            y_no{i,k} = yi;

            % attack measurements: add bias to one sensor after k_attack
            if (i == i_attacked) && (k >= k_attack)
                yi = yi + attack_bias;
            end
            y_att{i,k} = yi;
        end
    end

    % ---- Initial estimates (same for all filters) ----
    x_init = x0_true + [0.5; 0; 0.5; 0];   % biased initial guess

    x_hat_QET_no   = x_init;
    x_hat_QET_att  = x_init;
    x_hat_RQET_att = x_init;

    P0 = 10 * eye(4);
    P_QET_no   = P0;
    P_QET_att  = P0;
    P_RQET_att = P0;

    % ---- Time loop ----
    for k = 1:K

        % prediction
        if k == 1
            x_pred_QET_no   = x_hat_QET_no;   P_pred_QET_no   = P_QET_no;
            x_pred_QET_att  = x_hat_QET_att;  P_pred_QET_att  = P_QET_att;
            x_pred_RQET_att = x_hat_RQET_att; P_pred_RQET_att = P_RQET_att;
        else
            x_pred_QET_no   = A * x_hat_QET_no;
            P_pred_QET_no   = A * P_QET_no   * A.' + Q;

            x_pred_QET_att  = A * x_hat_QET_att;
            P_pred_QET_att  = A * P_QET_att  * A.' + Q;

            x_pred_RQET_att = A * x_hat_RQET_att;
            P_pred_RQET_att = A * P_RQET_att * A.' + Q;
        end

        %% ----- (1) QET without attack -----
        x_upd_QET_no = x_pred_QET_no;
        P_upd_QET_no = P_pred_QET_no;

        for i = 1:m
            Ci = C{i}; Ri = R{i}; yi = y_no{i,k};

            % triggering based on unquantized innovation
            e_trig = yi - Ci * x_upd_QET_no;
            S_trig = Ci * P_upd_QET_no * Ci.' + Ri;
            nis    = (e_trig.' / S_trig) * e_trig;

            if nis > gamma_QET
                y_q  = uniform_quantizer(yi, Delta, L_levels);
                R_eff = Ri + Rq;

                e_q = y_q - Ci * x_upd_QET_no;
                S_q = Ci * P_upd_QET_no * Ci.' + R_eff;
                K_i = P_upd_QET_no * Ci.' / S_q;

                x_upd_QET_no = x_upd_QET_no + K_i * e_q;
                P_upd_QET_no = P_upd_QET_no - K_i * Ci * P_upd_QET_no;
            end
        end

        %% ----- (2) QET with attack -----
        x_upd_QET_att = x_pred_QET_att;
        P_upd_QET_att = P_pred_QET_att;

        for i = 1:m
            Ci = C{i}; Ri = R{i}; yi = y_att{i,k};

            e_trig = yi - Ci * x_upd_QET_att;
            S_trig = Ci * P_upd_QET_att * Ci.' + Ri;
            nis    = (e_trig.' / S_trig) * e_trig;

            if nis > gamma_QET
                y_q  = uniform_quantizer(yi, Delta, L_levels);
                R_eff = Ri + Rq;

                e_q = y_q - Ci * x_upd_QET_att;
                S_q = Ci * P_upd_QET_att * Ci.' + R_eff;
                K_i = P_upd_QET_att * Ci.' / S_q;

                x_upd_QET_att = x_upd_QET_att + K_i * e_q;
                P_upd_QET_att = P_upd_QET_att - K_i * Ci * P_upd_QET_att;
            end
        end

        %% ----- (3) Resilient QET with attack -----
        x_upd_RQET_att = x_pred_RQET_att;
        P_upd_RQET_att = P_pred_RQET_att;

        for i = 1:m
            Ci = C{i}; Ri = R{i}; yi = y_att{i,k};

            % first-level trigger: NIS (same as QET)
            e_trig = yi - Ci * x_upd_RQET_att;
            S_trig = Ci * P_upd_RQET_att * Ci.' + Ri;
            nis    = (e_trig.' / S_trig) * e_trig;

            if nis > gamma_QET
                % robust check: normalized innovation
                norm_innov = abs(e_trig) / sqrt(S_trig);

                if norm_innov <= beta_max
                    % treat as normal, use quantized measurement
                    y_q  = uniform_quantizer(yi, Delta, L_levels);
                    R_eff = Ri + Rq;

                    e_q = y_q - Ci * x_upd_RQET_att;
                    S_q = Ci * P_upd_RQET_att * Ci.' + R_eff;
                    K_i = P_upd_RQET_att * Ci.' / S_q;

                    x_upd_RQET_att = x_upd_RQET_att + K_i * e_q;
                    P_upd_RQET_att = P_upd_RQET_att - K_i * Ci * P_upd_RQET_att;
                else
                    % suspiciously large innovation -> possible attack
                    % discard this measurement (no update)
                end
            end
        end

        % save updated states
        x_hat_QET_no   = x_upd_QET_no;    P_QET_no   = P_upd_QET_no;
        x_hat_QET_att  = x_upd_QET_att;   P_QET_att  = P_upd_QET_att;
        x_hat_RQET_att = x_upd_RQET_att;  P_RQET_att = P_upd_RQET_att;

        %% ----- RMSE in 2D position -----
        idx_pos = [1 3];   % x,y positions

        err_QET_no   = x_hat_QET_no(idx_pos)   - x_true(idx_pos,k);
        err_QET_att  = x_hat_QET_att(idx_pos)  - x_true(idx_pos,k);
        err_RQET_att = x_hat_RQET_att(idx_pos) - x_true(idx_pos,k);

        rmse_QET_no(k)   = rmse_QET_no(k)   + mean(err_QET_no.^2);
        rmse_QET_att(k)  = rmse_QET_att(k)  + mean(err_QET_att.^2);
        rmse_RQET_att(k) = rmse_RQET_att(k) + mean(err_RQET_att.^2);
    end
end

%% ---------- Average over MC ----------
rmse_QET_no   = sqrt(rmse_QET_no   / MC);
rmse_QET_att  = sqrt(rmse_QET_att  / MC);
rmse_RQET_att = sqrt(rmse_RQET_att / MC);

t = (0:K-1) * T;

%% ---------- Figure: RMSE vs time with attack ----------
figure;
plot(t, rmse_QET_no,   'k-',  'DisplayName','QET (no attack)'); hold on;
plot(t, rmse_QET_att,  'b--', 'DisplayName','QET (under attack)');
plot(t, rmse_RQET_att, 'r-',  'DisplayName','Resilient QET (under attack)');

% vertical line showing attack starting time
xline(t_attack, 'k:', 'Attack starts', ...
      'Interpreter','latex', ...
      'LabelHorizontalAlignment','left', ...
      'LabelVerticalAlignment','bottom');

xlabel('Time [s]');
ylabel('Position RMSE');
grid on; box on;
legend('Location','best');
title('Estimation Performance under Sensor Attack');

% slightly tighten y-limits to make differences clearer
yl = ylim;
ylim([0.9*yl(1), 1.05*yl(2)]);

% export vector figures for paper
print(gcf, '-depsc2', 'fig_Attack_RQET.eps');
print(gcf, '-dpdf',   'fig_Attack_RQET.pdf');
