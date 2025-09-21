clear all
q = -1.025;
% B = 1, 2, 3; mu = 0.1, G = 1; x_values = 11 for x = 0.3

G = 1;
% B = 0.01;

h = 1.47568;
h1 = 0.774801;%0.01
h2 = 0.597004;%0.02
h3 = 0.480137;%0.03

mu = 0.1;
ww1 = @(q, h, B, r, h1) - 1 + (((q + h) + (G*B/30*mu)*((mu - 1)*h1^5 + h^5)) / 2 ) * (3*(mu - 1)*h1^2 + (3*h^2 - 3*mu*r^2)) / ((h^3) + (mu - 1)*(h1^3)) - (G*B/24*mu)*((mu - 1)*h1^4 + (h^4 - mu*r^4));
ww = @(q, h, B, r, h1) - 1 + (((q + h) + (G*B/30*mu)*((mu - 1)*h1^5 + h^5)) / 2 ) * (3*h^2 - 3*r^2) / ((h^3) + (mu - 1)*(h1^3)) - (G*B/24*mu)*(h^4 - r^4);

B1 = 1;
r1_values = linspace(0, h1, 20);
for i = 1:length(r1_values)
    r = r1_values(i);
    w1 = ww1(q, h, B1, r, h1);
    
    % fprintf('%4.5f ', w1);
end
r_values1 = linspace(h1, h, 20);
for i = 1:length(r_values1)
    r = r_values1(i);
    w = ww(q, h, B1, r, h1);
    
    % fprintf('%4.5f ', w);
end

B2 = 2;
r2_values = linspace(0, h2, 20);
for i = 1:length(r2_values)
    r = r2_values(i);
    w2 = ww1(q, h, B2, r, h2);
    
    % fprintf('%4.5f ', w2);
end
r_values2 = linspace(h2, h, 20);
for i = 1:length(r_values2)
    r = r_values2(i);
    w = ww(q, h, B2, r, h2);
    
     % fprintf('%4.5f ', w);
end
B3 = 3;
r3_values = linspace(0, h3, 20);
for i = 1:length(r3_values)
    r = r3_values(i);
    w3 = ww1(q, h, B3, r, h3);
    
     fprintf('%4.5f ', w3);
end
r_values3 = linspace(h3, h, 20);
for i = 1:length(r_values2)
    r = r_values3(i);
    w = ww(q, h, B3, r, h3);
    
     fprintf('%4.5f ', w);
end
%Result
%0.001
A1 = [-0.60261 -0.60265 -0.60278 -0.60299 -0.60328 -0.60366 -0.60413 -0.60467 -0.60530 -0.60602 -0.60682 -0.60770 -0.60866 -0.60971 -0.61084 -0.61205 -0.61334 -0.61471 -0.61617 -0.61770];
A  = [-0.61770 -0.63226 -0.64746 -0.66330 -0.67978 -0.69690 -0.71464 -0.73300 -0.75199 -0.77159 -0.79179 -0.81260 -0.83400 -0.85600 -0.87858 -0.90173 -0.92546 -0.94975 -0.97460 -1.00000];
%0.002
B1 = [-0.58061 -0.58063 -0.58070 -0.58082 -0.58099 -0.58121 -0.58148 -0.58180 -0.58216 -0.58257 -0.58303 -0.58354 -0.58409 -0.58470 -0.58535 -0.58604 -0.58679 -0.58757 -0.58841 -0.58929];
B  = [-0.58929 -0.60307 -0.61781 -0.63351 -0.65015 -0.66771 -0.68618 -0.70554 -0.72577 -0.74686 -0.76878 -0.79152 -0.81505 -0.83936 -0.86441 -0.89018 -0.91665 -0.94380 -0.97159 -1.00000];
%0.005
C1 = [-0.56476 -0.56478 -0.56483 -0.56491 -0.56502 -0.56516 -0.56534 -0.56555 -0.56579 -0.56606 -0.56636 -0.56669 -0.56706 -0.56745 -0.56788 -0.56834 -0.56883 -0.56935 -0.56990 -0.57048];
C  = [-0.57048 -0.58344 -0.59766 -0.61312 -0.62978 -0.64763 -0.66662 -0.68673 -0.70793 -0.73016 -0.75340 -0.77760 -0.80272 -0.82871 -0.85552 -0.88309 -0.91139 -0.94034 -0.96990 -1.00000];

% Marker indices for each curve
idx1 = round(linspace(1, length(r1_values), 10));
idx2 = round(linspace(1, length(r2_values), 10));
idx3 = round(linspace(1, length(r3_values), 10));

idx4 = round(linspace(1, length(r_values1), 20));
idx5 = round(linspace(1, length(r_values2), 20));
idx6 = round(linspace(1, length(r_values3), 20));

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
legend('\beta = 1', '\beta = 2', '\beta = 3')
