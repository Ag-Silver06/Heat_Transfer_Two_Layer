clear all;

pi = 3.14; 
 G = 1; % Greshof number
  B = 0.01; % Heat Source Parameter
  T = 1;

h = @(x) 1+0.5 * sin(2*pi*x);
%1 + 0.5 * sin(2 * pi * x);
       
q =-1.025;
q1 = @(mu, q) (-(0.8)^8 * (mu-1) * G*B/120*mu + (0.8)^6 * (mu-1) * G*B/20*mu - (0.8)^5 * (mu-5/4) * G*B/30*mu - (0.8)^4 * (mu-1) * (1+G*B/24*mu) + (0.8)^3 * (mu-3/2) * (q + 1+ G*B/30*mu) + (0.8) * (G*B/(120*mu) + 3*q/2 + 1/2))/((mu-1)*(0.8)^3 + 1);




% Values of mu to consider
mu_values = [0.01, 0.1, 1, 10]; % Values of mu to consider
colors = {'-r', '-g', '-b', '-k'};

% Pre-allocate arrays to store the results
del_p_values = zeros(length(mu_values), 10); % 10 points for Q
Q_values = linspace(-2, 1.2, 100);

% Loop over different values of mu
for m = 1:length(mu_values)
     mu = mu_values(m);
    h1 = zeros(1, 20);  % Pre-allocate the h1 array
    fx = @(h, mu, h1) 1 /(h^3 + (mu - 1) * (h1^3));
    gx = @(h, mu, h1) (h^5) /(h^3 + (mu - 1) * (h1^3));
    lx = @(h, mu, h1) h /(h^3 + (mu - 1) * (h1^3));
    ix = @(h, mu, h1) h^2 / 2;
    x_values = linspace(0, 1, 10);

    for i = 1:length(x_values)
        x = x_values(i);
        h_value = h(x);
        q1_value = q1(mu, q);
        A = @(mu) (mu - 1) * G*B/(120*mu);
        B1 = @(mu, h) -(mu - 1) * G*B*h^2/(20*mu);
        C = @(mu, h) (mu - 5/4) * G*B*h^3/(30*mu);
        D = @(mu, h) (mu - 1) * (1 + G*B*h^4/(24*mu));
        E = @(mu, h, q1) -((mu - 3/2) * (q + h + G*B*h^5/(30*mu)) - q1*(mu - 1));
        F = @(h) -((q + h) * 3*h^2/2 + G*B*h^7/120 - h^3);
        G1 = @(q1, h) q1*h^3;

        a = A(mu);
        b1 = B1(mu, h_value);
        c = C(mu, h_value);
        d = D(mu, h_value);
        e = E(mu, h_value, q1_value);
        f = F(h_value);
        g1 = G1(q1_value, h_value);

        roots_array = roots([a 0 b1 c d e 0 f g1]);
        real_roots = roots_array(imag(roots_array) == 0 & real(roots_array) < h_value & real(roots_array) > 0);
        if ~isempty(real_roots)
            h1(i) = real_roots(1);
        end
    end

    % compute the integral parts
    for j = 1:length(x_values)
        x = x_values(j);
        h1_value = h1(j);
        h_value = h(x);
        f(j) = fx(h_value, mu, h1_value);
        g(j) = gx(h_value, mu, h1_value);
        l(j) = lx(h_value, mu, h1_value);
        i_val(j) = ix(h_value, mu, h1_value);
    end

    p = 1:(length(x_values)-1);
    I1 = (x_values(2) - x_values(1)) * sum(f(p) + f(p + 1)) / 2;
    I2 = (x_values(2) - x_values(1)) * sum(g(p) + g(p + 1)) / 2;
    I3 = (x_values(2) - x_values(1)) * sum(l(p) + l(p + 1)) / 2;
    I4 = (x_values(2) - x_values(1)) * sum(i_val(p) + i_val(p + 1)) / 2;

    % use average of non-zero h1 values
    h1_avg = mean(h1(h1 > 0));

    % corrected pressure function
    pp = @(mu, Q, h1avg, i1, i2, i3, i4) ...
         -3*mu*(Q - 0.5)*i1 ...
         - (G*B*(mu - 1)*(h1avg^5))*i1/10 ...
         - 3*mu*i3 ...
         - (G*B*i2/10) ...
         + G*(T + B*i4);

    for s = 1:length(Q_values)
        Q = Q_values(s);
        del_p = pp(mu, Q, h1_avg, I1, I2, I3, I4);
        del_p_values(m, s) = del_p;
    end
end

% Plotting
figure;
hold on;
for m = 1:length(mu_values)
    plot(Q_values, del_p_values(m, :), '-', 'Linewidth', 1.6, 'DisplayName', ['\mu = ' num2str(mu_values(m))]);
end
hold off;
xlabel('Q?');
ylabel('\Delta p');
title('\Delta p vs Q for Different \mu Values');
legend('Location', 'best');
grid on;
