function [hat_x_k1, P_k1] = calk_ekf(x, P, y)
    global a1 a2 dt
    hat_x = x;
    y_k1 = y;
    P_k = P;

    syms theta1 theta2 dtheta1 dtheta2
    x_sym = [theta1; theta2; dtheta1; dtheta2];
    f = x_sym + dt * f_c(x_sym);
    h = h_c(x_sym);
    A_sym = jacobian(f, x_sym);
    C_sym = jacobian(h, x_sym);

    % ヤコビアン計算
    A_k = double(subs(A_sym, x_sym, hat_x));
    C_k = double(subs(C_sym, x_sym, hat_x));

    W_k = A_k * P_k * C_k.' / (C_k * P_k * C_k.' + 0.01 ^ 2);
    hat_x_k1 = (hat_x + dt * f_c(hat_x)) + W_k * (y_k1 - h_c(hat_x));
    P_k1 = A_k * P_k * A_k' - W_k * (C_k * P_k * C_k' + 0.01 ^ 2) * W_k.';
end

function y = h_c(x)
    global a1 a2
    theta1 = x(1);
    theta2 = x(2);

    y = a1 * sin(theta1) + a2 * sin(theta1 + theta2);
end
