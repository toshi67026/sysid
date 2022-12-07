clear all;
close all;

l1 = 0.5; a1 = 1.0;
l2 = 0.25; a2 = 0.5;
m1 = 10; m2 = 5;
I1 = 5; I2 = 0.5;
d1 = 0; d2 = 0;
g = 9.8;

phi1 = @(x) m1 * l1 ^ 2 + x(5) * a1 ^ 2 + I1;
phi2 = @(x) x(5) * l2 ^ 2 + I2;
phi3 = @(x) x(5) * a1 * l2;

T = 0.01;
EndTime = 10;
t = 0:T:EndTime;
N = EndTime / T + 1;
n = 5; %óÔÏÌ
R = 1e-4; %ÏªG¹ÌªU

D = @(x) [phi1(x) + phi2(x) + 2 * phi3(x) * cos(x(2)) phi2(x) + phi3(x) * cos(x(2));
    phi2(x) + phi3(x) * cos(x(2)) phi2(x)];

H = @(x) (-phi3(x) * sin(x(2)) * [x(4) x(3) + x(4); -x(3) 0] + [d1 0; 0 d2]) * [x(3); x(4)];

G = @(x) g * [m1 * l1 * cos(x(1)) + x(5) * (a1 * cos(x(1)) + l2 * cos(x(1) + x(2)));
        x(5) * l2 * cos(x(1) + x(2))];

w = randn(1, N) * sqrtm(R);
x = zeros(n, N); y = zeros(1, N);
x(1, 1) = 0.0;
x(2, 1) = pi / 4;
x(3, 1) = 0;
x(4, 1) = 0;
x(5, 1) = m2;

% data generation
for i = 1:N - 1
    x(1, i + 1) = x(1, i) + x(3, i) * T;
    x(2, i + 1) = x(2, i) + x(4, i) * T;
    x(3:4, i + 1) = x(3:4, i) + inv(D(x(:, i))) * (-H(x(:, i)) - G(x(:, i))) * T;
    x(5, i + 1) = x(5, i);
    y(1, i) = l1 * cos(x(1, i)) + l2 * cos(x(1, i) + x(2, i)) + w(1, i);
end

y(1, N) = l1 * cos(x(1, N)) + l2 * cos(x(1, N) + x(2, N)) + w(1, N);

%------------------------------------------------
kappa = 0; % UT parameter
nn = 2 * n + 1; % total number of sigma points
nl = n + kappa;
xs = zeros(n, nn); % initializing of sigma points
xst = zeros(n, nn); %
% UKF
P0 = 0.1 * eye(n, n);
%P0(5,5)=1;
Pep = P0;
x0 = [0.0 pi / 4 0 0 m2 * 0.5]';
xep = x0;

% sigma points
[U, S, V] = svd(Pep);
SP = U * sqrt(S);
%SP=chol(Pep)';
ws = zeros(1, nn);

for j = 1:n
    xst(:, j) = xep + sqrt(nl) * SP(:, j);
    xst(:, j + n) = xep - sqrt(nl) * SP(:, j);
    ws(j) = 0.5 / nl;
    ws(j + n) = ws(j);
end

xst(:, nn) = xep;
ws(nn) = kappa / nl;
%
for i = 1:N
    % nonlinear transformation of sigma points
    yst = zeros(1, nn);

    for j = 1:nn
        yst(1, j) = l1 * cos(xst(1, j)) + l2 * cos(xst(1, j) + xst(2, j));
    end

    % mean and covariance of output y
    ys = zeros(1, 1);

    for j = 1:nn
        ys = ys + ws(j) * yst(1, j);
    end

    Pnu = zeros(1, 1);

    for j = 1:nn
        Pnu = Pnu + ws(j) * (yst(1, j) - ys) * (yst(1, j) - ys)';
    end

    %
    Pxnu = zeros(n, 1); % cross-covariance of x and nu

    for j = 1:nn
        Pxnu = Pxnu + ws(j) * (xst(:, j) - xep) * (yst(1, j) - ys)';
    end

    Pnu = Pnu + R; % covariance of innovation
    Kt = Pxnu * inv(Pnu); % UKF Kalman gain
    xef = xep + Kt * (y(1, i) - ys); % filtered estimate
    Pef = Pep - Kt * Pnu * Kt'; % filtered covariance matrix
    %
    xeff(:, i) = xef; % for figure
    Peff(:, :, i) = Pef; %
    Kff(:, i) = Kt; %

    %VO}|CgÌvZ
    % sigma points
    [U, S, V] = svd(Pef);
    SP = U * sqrt(S);
    %SP=chol(Pef)';
    for j = 1:n
        xs(:, j) = xef + sqrt(nl) * SP(:, j);
        xs(:, j + n) = xef - sqrt(nl) * SP(:, j);
    end

    xs(:, nn) = xef;

    %VO}|CgXV
    for j = 1:nn
        xst(1, j) = xs(1, j) + xs(3, j) * T;
        xst(2, j) = xs(2, j) + xs(4, j) * T;
        xst(3:4, j) = xs(3:4, j) + inv(D(xs(:, j))) * (-H(xs(:, j)) - G(xs(:, j))) * T;
        xst(5, j) = xs(5, j);
    end

    % prediction
    xt = zeros(n, 1);
    %OóÔèl
    for j = 1:nn
        xt = xt + ws(j) * xst(:, j); % predicted estimate
    end

    xep = xt;
    Ptt = zeros(n, n);
    %Oë·€ªUsñ
    for j = 1:nn
        Ptt = Ptt + ws(j) * (xst(:, j) - xep) * (xst(:, j) - xep)'; %prediction of covariance
    end

    %
    Pep = Ptt;
end % of i

%--------------------------------------------------

for i = 1:N
    RMSE(1, i) = sqrt((x(1, i) - xeff(1, i)) ^ 2);
    RMSE(2, i) = sqrt((x(2, i) - xeff(2, i)) ^ 2);
    RMSE(3, i) = sqrt((x(3, i) - xeff(3, i)) ^ 2);
    RMSE(4, i) = sqrt((x(4, i) - xeff(4, i)) ^ 2);
    RMSE(5, i) = sqrt((x(5, i) - xeff(5, i)) ^ 2);
end

%
RMSEukf1 = sum(RMSE(1, :)) / N
RMSEukf2 = sum(RMSE(2, :)) / N
RMSEukf3 = sum(RMSE(3, :)) / N
RMSEukf4 = sum(RMSE(4, :)) / N
RMSEukf5 = sum(RMSE(5, :)) / N

figure(1)
plot(t, x(1, :), 'r-', t, xeff(1, :), 'b-', 'LineWidth', 1.5)
grid
legend('x1', 'x1hat UKF')
xlabel('Time [sec]')
ylabel('Angle [rad]')
saveas(figure(1), 'x1_2.png', 'png');
%
figure(2)
plot(t, x(2, :), 'r-', t, xeff(2, :), 'b-', 'LineWidth', 1.5)
grid
legend('x2', 'x2hat UKF')
xlabel('Time [sec]')
ylabel('Angle [rad]')
saveas(figure(2), 'x2_2.png', 'png');
%
figure(3)
plot(t, x(3, :), 'r-', t, xeff(3, :), 'b-', 'LineWidth', 1.5)
grid
legend('x3', 'x3hat UKF')
xlabel('Time [sec]')
ylabel('Angular Velocity [rad/sec]')
saveas(figure(3), 'x3_2.png', 'png');
%
figure(4)
plot(t, x(4, :), 'r-', t, xeff(4, :), 'b-', 'LineWidth', 1.5)
grid
legend('x4', 'x4hat UKF')
xlabel('Time [sec]')
ylabel('Angular Velocity [rad/sec]')
saveas(figure(4), 'x4_2.png', 'png');
%
figure(5)
plot(t, x(5, :), 'r-', t, xeff(5, :), 'b-', 'LineWidth', 1.5)
grid
legend('x5', 'x5hat UKF')
xlabel('Time [sec]')
ylabel('Mass [kg]')
saveas(figure(5), 'x5_2.png', 'png');
%
figure(6)
plot(t, RMSE(1, :), 'k-', 'LineWidth', 1.5)
hold on
plot(t, RMSE(2, :), 'r-', 'LineWidth', 1.5)
plot(t, RMSE(3, :), 'g-', 'LineWidth', 1.5)
plot(t, RMSE(4, :), 'b-', 'LineWidth', 1.5)
plot(t, RMSE(5, :), 'y-', 'LineWidth', 1.5)
grid
xlabel('Time [sec]')
ylabel('E_t')
legend('x1 error', 'x2 error', 'x3 error', 'x4 error', 'x5 error')
saveas(figure(6), 'error_2.png', 'png');
