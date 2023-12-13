clc
clear
close all

% パラメータ設定
global l1 l2 a1 a2 m1 m2 d1 d2 g phi1 phi2 phi3
l1 = 0.5; a1 = 1.0;
l2 = 0.25; a2 = 0.5;
m1 = 10; m2 = 8;
I1 = 5; I2 = 0.5;
d1 = 0; d2 = 0;
g = 9.8;

phi1 = m1 * l1 ^ 2 + m2 * a1 ^ 2 + I1;
phi2 = m2 * l2 ^ 2 + I2;
phi3 = m2 * a1 * l2;

global dt
u = zeros(2, 1);
dt = 0.01;
x0 = [0; pi / 4; 0; 0];

% ode45を使う場合
% [t, x] = ode45(@(t, x) f_c(x, u), (0:dt:10), x0);
% x = x';

% 前進差分近似の場合
t = 0:dt:10;
x = [x0];

for i = 2:length(t)
    dx = f_c(x(:, end), u);
    x = [x, x(:, end) + dx * dt];
end

% ノイズ無しの場合の手先軌道
figure
hold on
x_ee = a1 * cos(x(1, :)) + a2 * cos(x(1, :) + x(2, :));
y_ee = a1 * sin(x(1, :)) + a2 * sin(x(1, :) + x(2, :));
plot(x_ee, y_ee, 'LineWidth', 3, 'DisplayName', 'ee')
xlabel('x')
ylabel('y')
grid on
box on
axis equal
set(gca, 'FontSize', 24)

% 白色ノイズを付加した場合
y = y_ee + wgn(1, length(t), -40);

%% (1)
global n kappa
n = 4;
kappa = 3 - n;
xest = x0 - ones(4, 1)/ 4;
% Sigmaポイントの初期共分散
P = 1e-2 * eye(4);

for i = 1:length(t) - 1
    [hat_x_k1_k1, P_k1_k1] = calc_ukf(xest(:, end), P, y(i)');
    xest = [xest, hat_x_k1_k1];
    P = P_k1_k1;
end

figure('Position', [100, 100, 1300, 900])

subplot(2, 2, 1)
hold on
plot(t, x(1, :), 'LineWidth', 3, 'DisplayName', 'True')
plot(t, xest(1, :), 'LineStyle', '--', 'LineWidth', 3, 'DisplayName', 'UKF')
xlabel('Time [s]')
ylabel('$\theta_1$', 'Interpreter', 'latex')
ylim([-6, 1])
grid on
box on
legend('Location', 'southwest')
set(gca, 'FontSize', 24)

subplot(2, 2, 2)
hold on
plot(t, x(2, :), 'LineWidth', 3, 'DisplayName', 'True')
plot(t, xest(2, :), 'LineStyle', '--', 'LineWidth', 3, 'DisplayName', 'UKF')
xlabel('Time [s]')
ylabel('$\theta_2$', 'Interpreter', 'latex')
grid on
box on
legend('Location', 'northwest')
set(gca, 'FontSize', 24)

subplot(2, 2, 3)
hold on
plot(t, x(3, :), 'LineWidth', 3, 'DisplayName', 'True')
plot(t, xest(3, :), 'LineStyle', '--', 'LineWidth', 3, 'DisplayName', 'UKF')
xlabel('Time [s]')
ylabel('$\dot{\theta}_1$', 'Interpreter', 'latex')
ylim([-10, 6])
grid on
box on
legend('Location', 'southwest')
set(gca, 'FontSize', 24)

subplot(2, 2, 4)
hold on
plot(t, x(4, :), 'LineWidth', 3, 'DisplayName', 'True')
plot(t, xest(4, :), 'LineStyle', '--', 'LineWidth', 3, 'DisplayName', 'UKF')
xlabel('Time [s]')
ylabel('$\dot{\theta}_2$', 'Interpreter', 'latex')
ylim([-30, 21])
grid on
box on
legend('Location', 'southwest')
set(gca, 'FontSize', 24)
