clear all;
close all;
clc;

% Constants and parameters
pi_val = 3.14;
figure;
  B= 0.01;
for   mu= [0.01 0.05 1]
    G = 0.5;
    beta = 0.01;
    QT = 2;
    a0 = 0.5;
    q = -1.025;

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
    x_values = linspace(-0.2, 1.75, 81);
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

    % Vertical points per layer
    Ny = 100;

    % Initialize streamfunction grid data
    E_combined = [];
    Y_combined = [];
    X_combined = [];

    % Loop to build streamfunction
    for i = 1:length(x_values)
        x = x_values(i);
        h1_val = h_all(i);
        h_val = h(x);

        y_core = linspace(0, h1_val, Ny);
        y_periph = linspace(h1_val, h_val, Ny);

        % Core layer
        term1 = -y_core;
        term2 = y_core./2 .* (q + h_val + G.*beta./(30.*mu) .* ((mu-1).*h1_val.^5 + h_val.^5)) .* ...
                ((3.*(mu-1).*h1_val.^2 + (3.*h_val.^2 - mu.*y_core.^2)) ./ ((mu-1).*h1_val.^3 + h_val.^3));
        term3 = -G.*beta.*y_core./(120.*mu) .* (5*(mu-1).*h1_val.^4 + (5.*h_val.^4 - mu.*y_core.^4));
        E1 = term1 + term2 + term3;

        % Peripheral layer
        term1 = -y_periph;
        term2 = (q + h_val + G.*beta./(30*mu) .* ((mu-1).*h1_val.^5 + h_val.^5)) .* ...
                ((mu-1).*h1_val.^3 + y_periph./2.*(3.*h_val.^2 - y_periph.^2)) ./ ((mu-1).*h1_val.^3 + h_val.^3);
        term3 = -G.*beta./(30.*mu) .* ((mu-1).*h1_val.^5 + y_periph./4 .* (5.*h_val.^4 - y_periph.^4));
        E2 = term1 + term2 + term3;

        % Stack data
        Y_combined = [Y_combined; y_core'; y_periph(2:end)'];
        X_combined = [X_combined; x * ones(Ny,1); x * ones(Ny-1,1)];
        E_combined = [E_combined; E1'; E2(2:end)'];
    end

    % Create regular meshgrid
    points_per_column = Ny * 2 - 1;
    x_vals = x_values;
    y_vals = linspace(0, max(Y_combined), points_per_column);
    [Xq, Yq] = meshgrid(x_vals, y_vals);

    % Interpolate streamfunction to meshgrid
    Egrid_interp = griddata(X_combined, Y_combined, E_combined, Xq, Yq, 'linear');
    Egrid_interp(isnan(Egrid_interp)) = 0;

    % Compute velocities from streamfunction
    dx_val = mean(diff(x_vals));
    dy_val = mean(diff(y_vals));
    [dy, dx] = gradient(Egrid_interp, dy_val, dx_val);
    u = dy;
    v = -dx;

    % === Multiple Particle Trajectories ===
    tf = 5;  % Increase simulation time
    time = linspace(0, tf, 201);  % More time steps
    dt = time(2) - time(1);

    % Define multiple initial positions
    initial_positions = [0.01, 0.01, 0,01];

    hold on;

    for k = 1:size(initial_positions, 1)
        X_pos = zeros(size(time));
        Y_pos = zeros(size(time));
        X_pos(1) = initial_positions(k,1);
        Y_pos(1) = initial_positions(k,2);

        for i = 1:length(time)-1
            u_interp = interp2(Xq, Yq, u, X_pos(i), Y_pos(i), 'linear', 0);
            v_interp = interp2(Xq, Yq, v, X_pos(i), Y_pos(i), 'linear', 0);

            X_pos(i+1) = X_pos(i) + dt * u_interp;
            Y_pos(i+1) = Y_pos(i) + dt * v_interp;
        end

        plot(X_pos, Y_pos, 'LineWidth', 1.5);
    end

    % Add streamlines to visualize the flow
%     contour(Xq, Yq, Egrid_interp, 30, 'k');  % Show streamlines

    % Add velocity arrows (quiver)
%     quiver(Xq, Yq, u, v, 2, 'r');  % Arrow scale 2

    xlabel('x');
    ylabel('y');
    grid off;
end