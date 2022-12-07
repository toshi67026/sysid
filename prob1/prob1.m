clc
clear
close all

mu = [1, 2];
Sigma = [1, -0.5; -0.5, 2];

%% plot pdf
figure
[X, Y] = meshgrid(linspace(-3, 5, 100), linspace(-1, 5, 100));
p = mvnpdf([X(:), Y(:)], mu, Sigma);
Z = reshape(p, 100, 100);
surf(X, Y, Z);
grid on
box on
xlabel('x')
ylabel('y')
zlabel('z')
set(gca, 'FontSize', 24)

num_samples = 10;
samples = mvnrnd(mu, Sigma, num_samples);
x = samples(:, 1);
y = samples(:, 2);
bar_x = mean(x);
bar_y = mean(y);

%% x,y
figure
hold on
plot(x, y, 'LineStyle', 'none', 'LineWidth', 3, 'Marker', 'o', 'MarkerSize', 10, 'DisplayName', 'samples')
contour(X, Y, Z, 'LineWidth', 3, 'DisplayName', 'pdf')
legend
axis equal
grid on
box on
set(gca, 'FontSize', 24)

%% 2
[V, D] = eig(Sigma);
hat_x_2 = V(1, 2) / V(2, 2) * (y - bar_y) + bar_x;

figure
hold on
plot(x, y, 'LineStyle', 'none', 'LineWidth', 3, 'Marker', 'o', 'MarkerSize', 10, 'DisplayName', 'samples')
plot(hat_x_2, y, 'LineStyle', 'none', 'LineWidth', 3, 'Marker', 'o', 'MarkerSize', 10, 'DisplayName', 'est')
contour(X, Y, Z, 'LineWidth', 3, 'DisplayName', 'pdf')
legend
axis equal
grid on
box on
set(gca, 'FontSize', 24)

%% 3
% 条件付き確率分布のhatx
A0 = Sigma(1, 2) / Sigma(2, 2);
hat_x_3 = bar_x + A0 * (y - bar_y);

figure
hold on
plot(x, y, 'LineStyle', 'none', 'LineWidth', 3, 'Marker', 'o', 'MarkerSize', 10, 'DisplayName', 'samples')
plot(hat_x_3, y, 'LineStyle', 'none', 'LineWidth', 3, 'Marker', 'o', 'MarkerSize', 10, 'DisplayName', 'est')
contour(X, Y, Z, 'LineWidth', 3, 'DisplayName', 'pdf')
legend
axis equal
grid on
box on
set(gca, 'FontSize', 24)

%% 4
E_e_2 = sum(x - hat_x_2) / num_samples
V_e_2 = sum((x - hat_x_2) .^ 2) / (num_samples - 1)

E_e_3 = sum(x - hat_x_3) / num_samples
V_e_3 = sum((x - hat_x_3) .^ 2) / (num_samples - 1)

figure
hold on
plot(x, 'LineStyle', 'none', 'LineWidth', 3, 'Marker', 'o', 'MarkerSize', 10, 'DisplayName', 'True')
plot(hat_x_2, 'LineStyle', 'none', 'LineWidth', 3, 'Marker', 'o', 'MarkerSize', 10, 'DisplayName', '2')
plot(hat_x_3, 'LineStyle', 'none', 'LineWidth', 3, 'Marker', 'o', 'MarkerSize', 10, 'DisplayName', '3')
yline(bar_x, 'LineWidth', 3, 'DisplayName', 'bar\_x')
xlim([0, 10])
legend
grid on
box on
set(gca, 'FontSize', 24)
