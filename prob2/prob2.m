clc
clear
close all

% パラメータ設定
a = 1 % alpha
B = 4;
% 区間[0, 10]をN=1000を用いてdt:=10/1000として離散化
dt = 0.01;
t = 0:dt:10

%{
dP(t) = r P(t)dt + a P(t) dB(t)
      = (rdt + a dB(t)) P(t)
%}

%% r = 1.0, 2.0, 3.0についてそれぞれ計算
for r = 1:3
    figure('Position', [100, 100, 1300, 1000])

    sgtitle("r = " + num2str(r), 'FontSize', 24)

    % 10回計算
    for i = 1:10
        P(1) = 10;
        subplot(5, 2, i)

        % Euler-Maruyama methods
        for j = 2:length(t)
            f_j_1 = r * P(j - 1);
            g_j_1 = a * P(j - 1);
            P(j) = P(j - 1) + f_j_1 * dt + g_j_1 * normrnd(0, sqrt(B * dt));
        end

        plot(t, P, 'LineWidth', 3)
        grid on
        set(gca, 'FontSize', 18)
    end

end

% logスケールで一枚
for r = 1:3
    figure('Position', [100, 100, 1300, 1000])
    hold on
    title("r = " + num2str(r), 'FontSize', 24)

    % 10回計算
    for i = 1:10
        P(1) = 10;

        % Euler-Maruyama methods
        for j = 2:length(t)
            f_j_1 = r * P(j - 1);
            g_j_1 = a * P(j - 1);
            P(j) = P(j - 1) + f_j_1 * dt + g_j_1 * normrnd(0, sqrt(B * dt));
        end

        % 常用対数
        % 点線見づらいから薄くすると良いかも？
        plot(t, log10(P), 'LineStyle', '--', 'LineWidth', 3)
        grid on
        set(gca, 'FontSize', 18)
    end

    % Pの平均
end
