function [hat_x_k1_k1, P_k1_k1] = calk_ukf(x, P, y)
    global n kappa dt
    hat_x = x;
    y_k1 = y;

    % M:=(n + kappa)P(k|k))とした時のM=NN^\topを満たす行列Nのi番目の列ベクトルが\sqrt(M)
    M = (n + kappa) * P;
    [U, S, V] = svd(M);
    N = U * sqrt(S);

    X = [];
    W = [];

    %% Step 1
    X(:, 1) = hat_x;
    W(1) = kappa / (n + kappa);

    % MATLABのindexの関係で1ずらす
    for i = 1:n
        X(:, i + 1) = hat_x + N(:, i);
        W(i + 1) = 1 / (2 * (n + kappa));
        X(:, i + 1 + n) = hat_x - N(:, i);
        W(i + 1 + n) = 1 / (2 * (n + kappa));
    end

    %% Step 2
    X_k1_k = [];
    hat_x_k1_k = zeros(n, 1);
    P_k1_k = zeros(n);

    for i = 1:2 * n + 1
        X_k1_k(:, i) = X(:, i) + dt * f_c(X(:, i));
        hat_x_k1_k = hat_x_k1_k + W(i) * X_k1_k(:, i);
    end

    for i = 1:2 * n + 1
        P_k1_k = P_k1_k + W(i) * ((X_k1_k(:, i) - hat_x_k1_k) * (X_k1_k(:, i) - hat_x_k1_k).');
    end

    %% Step 3
    eta_k1_k = [];

    global a1 a2

    for i = 1:2 * n + 1
        x = X_k1_k(:, i);
        eta_k1_k(i) = a1 * sin(x(1)) + a2 * sin(x(1) + x(2));
    end

    hat_y_k1_k = sum(W .* eta_k1_k);
    P_yy = 0;

    for i = 1:2 * n + 1
        P_yy = P_yy + W(i) * (eta_k1_k(i) - hat_y_k1_k) * (eta_k1_k(i) - hat_y_k1_k)';
    end

    %% Step 4
    % これらはk+1|k
    P_vv = 0.01 ^ 2 + P_yy;
    P_xv = zeros(n, 1);

    for i = 1:2 * n + 1
        P_xv = P_xv + W(i) * (X(:, i) - hat_x_k1_k) * (eta_k1_k(i) - hat_y_k1_k).';
    end

    %% Step 5
    W_k1 = P_xv * inv(P_vv);
    v_k1 = y_k1 - hat_y_k1_k; % 予測出力誤差

    hat_x_k1_k1 = hat_x_k1_k + W_k1 * v_k1;
    P_k1_k1 = P_k1_k - W_k1 * P_vv * W_k1.';
end
