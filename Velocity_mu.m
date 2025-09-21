clear all
q = -1.025;
% mu = 0.1, 1, 10; G = 1, B = 0.01; x_values = 11 for x = 0.3

G = 0.5;
B = 0.01;

h = 1.47568;
h1 = 1.27367;%0.1
h2 = 1.26063;%1
h3 = 1.23123;%10

mu1 = 0.1;
ww1 = @(q, h, mu, r, h1) - 1 + (((q + h) + (G*B/30*mu)*((mu - 1)*h1^5 + h^5)) / 2 ) * (3*(mu - 1)*h1^2 + (3*h^2 - 3*mu*r^2)) / ((h^3) + (mu - 1)*(h1^3)) - (G*B/24*mu)*((mu - 1)*h1^4 + (h^4 - mu*r^4));
ww = @(q, h, mu, r, h1) - 1 + (((q + h) + (G*B/30*mu)*((mu - 1)*h1^5 + h^5)) / 2 ) * (3*h^2 - 3*r^2) / ((h^3) + (mu - 1)*(h1^3)) - (G*B/24*mu)*(h^4 - r^4);
r1_values = linspace(0, h1, 20);
for i = 1:length(r1_values)
    r = r1_values(i);
    w1 = ww1(q, h, mu1, r, h1);
    
    % fprintf('%4.5f ', w1);
end
r_values1 = linspace(h1, h, 20);
for i = 1:length(r_values1)
    r = r_values1(i);
    w = ww(q, h, mu1, r, h1);
    
    % fprintf('%4.5f ', w);
end

mu2 = 1;
r2_values = linspace(0, h2, 20);
for i = 1:length(r2_values)
    r = r2_values(i);
    w2 = ww1(q, h, mu2, r, h2);
    
    % fprintf('%4.5f ', w2);
end
r_values2 = linspace(h2, h, 20);
for i = 1:length(r_values2)
    r = r_values2(i);
    w = ww(q, h, mu2, r, h2);
    
     % fprintf('%4.5f ', w);
end
mu3 = 10;
r3_values = linspace(0, h3, 20);
for i = 1:length(r3_values)
    r = r3_values(i);
    w3 = ww1(q, h, mu3, r, h3);
    
     % fprintf('%4.5f ', w3);
end
r_values3 = linspace(h3, h, 20);
for i = 1:length(r_values2)
    r = r_values3(i);
    w = ww(q, h, mu3, r, h3);
    
     fprintf('%4.5f ', w);
end
%Result
%0.001
A1 = [-0.64168 -0.64191 -0.64258 -0.64370 -0.64527 -0.64729 -0.64976 -0.65268 -0.65605 -0.65986 -0.66412 -0.66884 -0.67400 -0.67961 -0.68567 -0.69217 -0.69913 -0.70653 -0.71439 -0.72269];
A  = [-0.72269 -0.73627 -0.74996 -0.76377 -0.77769 -0.79172 -0.80586 -0.82012 -0.83449 -0.84897 -0.86356 -0.87827 -0.89309 -0.90803 -0.92307 -0.93823 -0.95351 -0.96889 -0.98439 -1.00000];
%0.002
B1 = [-0.54169 -0.54262 -0.54541 -0.55005 -0.55655 -0.56490 -0.57512 -0.58718 -0.60110 -0.61688 -0.63450 -0.65398 -0.67531 -0.69849 -0.72352 -0.75040 -0.77912 -0.80969 -0.84210 -0.87635];
B  = [-0.87635 -0.88238 -0.88846 -0.89459 -0.90078 -0.90702 -0.91331 -0.91966 -0.92606 -0.93252 -0.93902 -0.94559 -0.95220 -0.95887 -0.96559 -0.97237 -0.97919 -0.98608 -0.99301 -1.00000];
%0.005
C1 = [-0.45435 -0.45594 -0.46070 -0.46862 -0.47968 -0.49384 -0.51107 -0.53132 -0.55453 -0.58064 -0.60956 -0.64122 -0.67553 -0.71238 -0.75166 -0.79325 -0.83702 -0.88285 -0.93058 -0.98005];
C  = [-0.98005 -0.98106 -0.98206 -0.98308 -0.98410 -0.98512 -0.98615 -0.98719 -0.98823 -0.98928 -0.99033 -0.99139 -0.99245 -0.99352 -0.99459 -0.99566 -0.99674 -0.99782 -0.99891 -1.00000];




% Marker indices for each curve
idx1 = round(linspace(1, length(r1_values), 20));
idx2 = round(linspace(1, length(r2_values), 20));
idx3 = round(linspace(1, length(r3_values), 20));

idx4 = round(linspace(1, length(r_values1), 12));
idx5 = round(linspace(1, length(r_values2), 8));
idx6 = round(linspace(1, length(r_values3), 4));

figure(2)
hold on

% Plot with markers at selected points only
plot(A1, r1_values, 'o-', 'MarkerIndices', idx1);
plot(B1, r2_values, 's-', 'MarkerIndices', idx2);
plot(C1, r3_values, 'd-', 'MarkerIndices', idx3);

plot(A, r_values1, 'o-', 'MarkerIndices', idx4);
plot(B, r_values2, 's-', 'MarkerIndices', idx5);
plot(C, r_values3, 'd-', 'MarkerIndices', idx6);

xlabel('u')
ylabel('y')
legend('\mu = 0.1', '\mu = 1', '\mu = 10', 'Location', 'best')


