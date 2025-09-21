clear all
q = -1.025;
% G = 1, 2, 3; mu = 0.1, B = 0.1; x_values = 11 for x = 0.3

% G = 1;
B = 0.1;

h = 1.47568;
h1 = 1.13881;%1
h2 = 1.06926;%2
h3 = 1.01338;%3
h4 = 1.23758;%0

mu = 0.1;
ww1 = @(q, h, G, r, h1) - 1 + (((q + h) + (G*B/30*mu)*((mu - 1)*h1^5 + h^5)) / 2 ) * (3*(mu - 1)*h1^2 + (3*h^2 - 3*mu*r^2)) / ((h^3) + (mu - 1)*(h1^3)) - (G*B/24*mu)*((mu - 1)*h1^4 + (h^4 - mu*r^4));
ww = @(q, h, G, r, h1) - 1 + (((q + h) + (G*B/30*mu)*((mu - 1)*h1^5 + h^5)) / 2 ) * (3*h^2 - 3*r^2) / ((h^3) + (mu - 1)*(h1^3)) - (G*B/24*mu)*(h^4 - r^4);

G1 = 1;
r1_values = linspace(0, h1, 20);
for i = 1:length(r1_values)
    r = r1_values(i);
    w1 = ww1(q, h, G1, r, h1);
    
    % fprintf('%4.5f ', w1);
end
r_values1 = linspace(h1, h, 20);
for i = 1:length(r_values1)
    r = r_values1(i);
    w = ww(q, h, G1, r, h1);
    
    % fprintf('%4.5f ', w);
end


G2 = 2;
r2_values = linspace(0, h2, 20);
for i = 1:length(r2_values)
    r = r2_values(i);
    w2 = ww1(q, h, G2, r, h2);
    
    % fprintf('%4.5f ', w2);
end
r_values2 = linspace(h2, h, 20);
for i = 1:length(r_values2)
    r = r_values2(i);
    w = ww(q, h, G2, r, h2);
    
     % fprintf('%4.5f ', w);
end


G3 = 3;
r3_values = linspace(0, h3, 20);
for i = 1:length(r3_values)
    r = r3_values(i);
    w3 = ww1(q, h, G3, r, h3);
    
     % fprintf('%4.5f ', w3);
end
r_values3 = linspace(h3, h, 20);
for i = 1:length(r_values3)
    r = r_values3(i);
    w = ww(q, h, G3, r, h3);
    
     % fprintf('%4.5f ', w);
end


G4=0;
r4_values = linspace(0, h4, 20);
for i = 1:length(r4_values)
    r = r4_values(i);
    w4 = ww1(q, h, G4, r, h4);
    
     fprintf('%4.5f ', w4);
end
r_values4 = linspace(h4, h, 20);
for i = 1:length(r_values4)
    r = r_values4(i);
    w = ww(q, h, G4, r, h4);
    
     fprintf('%4.5f ', w);
end


%Result
%0.001
A1 = [-0.63742 -0.63755 -0.63794 -0.63858 -0.63949 -0.64065 -0.64208 -0.64376 -0.64570 -0.64790 -0.65035 -0.65307 -0.65604 -0.65927 -0.66276 -0.66650 -0.67051 -0.67477 -0.67928 -0.68406];
A  = [-0.68406 -0.69867 -0.71351 -0.72857 -0.74386 -0.75937 -0.77510 -0.79106 -0.80724 -0.82365 -0.84028 -0.85713 -0.87421 -0.89151 -0.90903 -0.92678 -0.94475 -0.96294 -0.98136 -1.00000];
%0.002
B1 = [-0.63240 -0.63250 -0.63281 -0.63332 -0.63404 -0.63496 -0.63608 -0.63741 -0.63894 -0.64067 -0.64261 -0.64475 -0.64710 -0.64965 -0.65240 -0.65535 -0.65850 -0.66186 -0.66542 -0.66918];
B  = [-0.66918 -0.68399 -0.69910 -0.71449 -0.73018 -0.74615 -0.76241 -0.77896 -0.79580 -0.81293 -0.83035 -0.84805 -0.86604 -0.88432 -0.90288 -0.92173 -0.94087 -0.96030 -0.98000 -1.00000];
%0.005
C1 = [-0.62762 -0.62770 -0.62796 -0.62839 -0.62899 -0.62976 -0.63070 -0.63181 -0.63309 -0.63454 -0.63617 -0.63796 -0.63992 -0.64205 -0.64436 -0.64683 -0.64946 -0.65227 -0.65524 -0.65839];
C  = [-0.65839 -0.67327 -0.68850 -0.70408 -0.72000 -0.73627 -0.75288 -0.76984 -0.78714 -0.80479 -0.82278 -0.84111 -0.85978 -0.87879 -0.89815 -0.91784 -0.93787 -0.95824 -0.97895 -1.00000];
%0
D1 = [-0.64163 -0.64182 -0.64239 -0.64334 -0.64467 -0.64638 -0.64848 -0.65095 -0.65380 -0.65704 -0.66065 -0.66465 -0.66902 -0.67378 -0.67892 -0.68443 -0.69033 -0.69661 -0.70327 -0.71031];
D  = [-0.71031 -0.72429 -0.73841 -0.75267 -0.76707 -0.78161 -0.79630 -0.81112 -0.82609 -0.84119 -0.85644 -0.87183 -0.88736 -0.90303 -0.91884 -0.93479 -0.95088 -0.96711 -0.98349 -1.00000];
% Marker indices for each curve
idx1 = round(linspace(1, length(r1_values), 14));
idx2 = round(linspace(1, length(r2_values), 14));
idx3 = round(linspace(1, length(r3_values), 14));
idx7 = round(linspace(1, length(r4_values), 14));

idx4 = round(linspace(1, length(r_values1), 20));
idx5 = round(linspace(1, length(r_values2), 20));
idx6 = round(linspace(1, length(r_values3), 20));
idx8 = round(linspace(1, length(r_values4), 20));

figure(2)
hold on

% Plot with markers at selected points only
plot(D1, r4_values, '-', 'MarkerIndices', idx7);
plot(A1, r1_values, '-', 'MarkerIndices', idx1);
plot(B1, r2_values, '-', 'MarkerIndices', idx2);
plot(C1, r3_values, '-', 'MarkerIndices', idx3);


plot(D, r_values4, '-', 'MarkerIndices', idx8);
plot(A, r_values1, '-', 'MarkerIndices', idx4);
plot(B, r_values2, '-', 'MarkerIndices', idx5);
plot(C, r_values3, '-', 'MarkerIndices', idx6);

xlabel('u')
ylabel('y')
legend( 'G = 0','G = 1', 'G = 2', 'G = 3')
