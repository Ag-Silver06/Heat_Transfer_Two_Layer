clear all;
close all;
clc;

% Constants and parameters
pi_val = 3.14;  
G = 0;       
B = 0;      
mu = 0.1;
beta = B;
QT=0.3;
  a0=0.5;       
q = QT-a0;

% Function definitions
h = @(x) 1 + a0 * sin(2 * pi_val * x);
q1 = @(mu, q) (-(0.8)^8 * (mu - 1) * G * B / (120 * mu) + ...
               (0.8)^6 * (mu - 1) * G * B / (20 * mu) - ...
               (0.8)^5 * (mu - 5/4) * G * B / (30 * mu) - ...
               (0.8)^4 * (mu - 1) * (1 + G * B / (24 * mu)) + ...
               (0.8)^3 * (mu - 3/2) * (q + 1 + G * B / (30 * mu)) + ...
               (0.8) * (G * B / (120 * mu) + 3 * q / 2 + 1 / 2)) / ...
               ((mu - 1) * (0.8)^3 + 1);

% X values
x_values = linspace(0,1, 81);
h_all = zeros(size(x_values));
q1_val = q1(mu, q);

% Compute h1(x)
for i = 1:length(x_values)
    x = x_values(i);
    h_val = h(x);

    A = (mu - 1) * G * B / (120 * mu);
    B1 = -(mu - 1) * G * B * h_val^2 / (20 * mu);
    C = (mu - 5/4) * G * B * h_val^3 / (30 * mu);
    D = (mu - 1) * (1 + G * B * h_val^4 / (24 * mu));
    E = -((mu - 3/2) * (q + h_val + G * B * h_val^5 / (30 * mu)) - q1_val * (mu - 1));
    F = -((q + h_val) * 3 * h_val^2 / 2 + G * B * h_val^7 / 120 - h_val^3);
    G1 = q1_val * h_val^3;

    coeffs = [A, 0, B1, C, D, E, 0, F, G1];
    rts = roots(coeffs);

    real_rts = rts(imag(rts) == 0 & real(rts) > 0 & real(rts) < h_val);

    if isempty(real_rts)
        h_all(i) = 0;
    else
        h_all(i) = real_rts(1);
    end
end

% Number of vertical points for each layer
Ny = 100;

% Initialize surfaces
E_combined = [];
Y_combined = [];
X_combined = [];

% Loop to construct seamless layers
for i = 1:length(x_values)
    x = x_values(i);
    h1_val = h_all(i);
    h_val = h(x);
    
    y_core = linspace(0, h1_val, Ny);
    y_periph = linspace(h1_val, h_val, Ny);
    
     % Core layer
            term1 = -y_core;
            term2 = y_core./2 .* (q + h_val + G.*beta./(30.*mu) .* ((mu-1).*h1_val.^5 + h_val^5)) .* ...
                    ((3.*(mu-1).*h1_val.^2 + (3.*h_val.^2 - mu.*y_core.^2)) ./ ((mu-1).*h1_val.^3 + h_val.^3));
            term3 = -G.*beta.*y_core./(120.*mu) .* (5*(mu-1).*h1_val.^4 + (5.*h_val.^4 - mu.*y_core.^4));
            E1 = term1 + term2 + term3;
    
%     E1 = (beta .* y_core).^2 + ...
%      (Br ./ T) .* ( ...
%      -3 .* mu .* y_core .* ...
%      (q + h_val + (G * beta / (30 * mu)) .* ...
%      ((mu - 1) .* h1_val.^5 + h_val.^5)) ./ ...
%      ((mu - 1) .* h1_val.^3 + h_val.^3) ).^2;




% Peripheral layer
            term1 = -y_periph;
            term2 = (q + h_val + G.*beta./(30*mu) .* ((mu-1).*h1_val.^5 + h_val.^5)) .* ...
                    ((mu-1).*h1_val.^3 + y_periph./2.*(3.*h_val.^2 - y_periph.^2)) ./ ((mu-1).*h1_val.^3 + h_val.^3);
            term3 = -G.*beta./(30.*mu) .* ((mu-1).*h1_val.^5 + y_periph./4 .* (5.*h_val.^4 -y_periph.^4));
             E2 = term1 + term2 + term3;
            
            
            
%     E2 = (beta .* y_periph).^2 + ...
%     (Br ./ T) .* ( ...
%     -3 .* mu .* y_periph .* ...
%     (q + h_val + (G * beta / (30 * mu)) .* ...
%     ((mu - 1) .* h1_val.^5 + h_val.^5)) ./ ...
%     ((mu - 1) .* h1_val.^3 + h_val.^3) + (G * beta / (6 * mu)) .* y_periph.^3 ).^2;
    
    % Stack values
    Y_combined = [Y_combined; y_core'; y_periph(2:end)'];
    X_combined = [X_combined; x * ones(Ny,1); x * ones(Ny-1,1)];
    E_combined = [E_combined; E1'; E2(2:end)'];
end

% Reshape for contour plot
points_per_column = Ny * 2 - 1;
Xgrid = reshape(X_combined, points_per_column, []);
Ygrid = reshape(Y_combined, points_per_column, []);
Egrid = reshape(E_combined, points_per_column, []);

% === Contour Plot ===
figure(3);
contourf(Xgrid, Ygrid, Egrid, 20);
xlabel('x');
ylabel('y');
% title('Contour plot of E_G(x,y)');
colorbar;
% colormap jet;
% axis tight;
box on;

% Overlay h1(x) and h(x) on contour
hold on;
plot(x_values, h_all, 'k--', 'LineWidth', 2);  % h1 interface
plot(x_values, arrayfun(h, x_values), 'k-', 'LineWidth', 2);  % wall profile

% legend('E_G Contours', 'Interface h_1(x)', 'Wall h(x)', 'Location', 'northeast');