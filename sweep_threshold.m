function [rmse_ss, comm_avg] = sweep_threshold(A,Q,C,R,...
    x0_true, T, gamma_QET, Delta, L_levels, K, MC)
%SWEEP_THRESHOLD  通信–精度折中仿真（Quantized ET）
%
%   返回:
%   rmse_ss  : 稳态 RMSE（最后 1/3 时间平均）
%   comm_avg : 每个传感器每步的平均通信次数

    m  = numel(C);
    Rq = (Delta^2)/12;

    rmse_time = zeros(K,1);
    comm_time = zeros(K,1);

    for mc = 1:MC

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

        % ----- Initial estimate -----
        x_hat = x0_true + [0.5; 0; 0.5; 0];
        P_hat = 10 * eye(4);

        for k = 1:K

            % prediction
            if k == 1
                x_pred = x_hat;
                P_pred = P_hat;
            else
                x_pred = A * x_hat;
                P_pred = A * P_hat * A.' + Q;
            end

            x_upd   = x_pred;
            P_upd   = P_pred;
            comm_sc = 0;

            for i = 1:m
                Ci = C{i}; Ri = R{i}; yi = y{i,k};

                % trigger based on unquantized innovation
                e_trig = yi - Ci * x_upd;
                S_trig = Ci * P_upd * Ci.' + Ri;
                nis    = (e_trig.' / S_trig) * e_trig;

                if nis > gamma_QET
                    % quantize measurement
                    y_q  = uniform_quantizer(yi, Delta, L_levels);
                    R_eff = Ri + Rq;

                    e_q = y_q - Ci * x_upd;
                    S_q = Ci * P_upd * Ci.' + R_eff;
                    K_i = P_upd * Ci.' / S_q;

                    x_upd = x_upd + K_i * e_q;
                    P_upd = P_upd - K_i * Ci * P_upd;

                    comm_sc = comm_sc + 1;
                end
            end

            x_hat = x_upd;
            P_hat = P_upd;

            err = x_hat([1 3]) - x_true([1 3],k);
            rmse_time(k) = rmse_time(k) + mean(err.^2);
            comm_time(k) = comm_time(k) + comm_sc;
        end
    end

    % average over MC
    rmse_time = sqrt(rmse_time / MC);
    comm_time = comm_time / MC;

    % last 1/3 horizon as steady state
    idx_ss = round(2*K/3):K;
    rmse_ss  = mean(rmse_time(idx_ss));
    comm_avg = mean(comm_time(idx_ss)) / m;
end
