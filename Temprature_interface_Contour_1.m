clear all; clc; close all;

pi = 3.14;
G = 1; B = 0.05; beta = B;
q = -1.025;

h = @(x) 1 + 0.4 * sin(2 * pi * x);
q1 = @(mu, q) (-(0.8)^8 * (mu-1) * G*B/120*mu + (0.8)^6 * (mu-1) * G*B/20*mu - (0.8)^5 * (mu-5/4) * G*B/30*mu - (0.8)^4 * (mu-1) * (1+G*B/24*mu) + (0.8)^3 * (mu-3/2) * (q + 1+ G*B/30*mu) + (0.8) * (G*B/(120*mu) + 3*q/2 + 1/2))/((mu-1)*(0.8)^3 + 1);

% Grid settings
Nx = 1010;
Ny = 1500;
x_vals = linspace(0, 1, Nx);
y_vals = linspace(0, 1.6, Ny);
[X, Y] = meshgrid(x_vals, y_vals);

% --- Get h_all(x) ---
mu = 0.1;
q1_val = q1(mu, q);

x_all = []; h_all = [];

for i = 1:Nx
    x = x_vals(i);
    h_val = h(x);

    A = (mu - 1) * G*B/(120*mu);
    B1 = -(mu - 1) * G*B*h_val^2/(20*mu);
    C = (mu - 5/4) * G*B*h_val^3/(30*mu);
    D = (mu -1) * (1 + G*B*h_val^4/(24*mu));
    E = - ((mu - 3/2) * (q + h_val + G*B*h_val^5/(30*mu)) - q1_val*(mu - 1));
    F = -((q + h_val) * 3*h_val^2/2 + G*B*h_val^7/120 - h_val^3);
    G1 = q1_val*h_val^3;

    roots_array = roots([A 0 B1 C D E 0 F G1]);
    real_roots = roots_array(imag(roots_array) == 0 & real(roots_array) < h_val & real(roots_array) > 0);

    if ~isempty(real_roots)
        h_all = [h_all, real_roots(1)];
    else
        h_all = [h_all, 0];
    end
    x_all = [x_all, x];
end

% Interpolated h_all values across grid
H = arrayfun(h, x_vals);
H_all = interp1(x_vals, h_all, X(1,:), 'linear', 'extrap'); % 1st row = x
H_all_grid = repmat(H_all, Ny, 1);
H_grid = h(X);

% --- Temperature Calculations ---
T_lower = nan(size(X));
T_upper = nan(size(X));

for i = 1:Nx
    for j = 1:Ny
        if Y(j,i) <= H_all_grid(j,i)
            T_lower(j,i) = 1 - (beta/2) * (Y(j,i)^2 - H_grid(j,i)^2);
        elseif Y(j,i) <= H_grid(j,i)
            T_upper(j,i) = 1 - (beta/2) * (Y(j,i)^2 - H_grid(j,i)^2);
        end
    end
end

% --- Plotting both contour regions ---
h_vals = arrayfun(h, x_vals);  % <-- add this to define h_vals

% --- Plotting both contour regions ---
figure;
hold on;

contourf(X, Y, T_lower, 50, 'LineColor', 'none');
contourf(X, Y, T_upper, 50, 'LineColor', 'none');

colormap jet;
colorbar;
plot(x_vals, h_vals, 'k', 'LineWidth', 2, 'DisplayName', 'Wall h(x)');
plot(x_vals, h_all, '--r', 'LineWidth', 1.5, 'DisplayName', 'Interface h_{all}(x)');

xlabel('x');
ylabel('y');
title('Temperature Contours from 0 to h_{all}(x) and h_{all}(x) to h(x)');
legend('Location', 'northeast');
axis([0 1 0 1.6]);
