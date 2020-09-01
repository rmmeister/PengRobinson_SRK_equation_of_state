%% INITIALIZE MATLAB
clear all;
clc;
close all;
format long;

%% DASHBOARD

% problem parameters
T = 150 + 273.15; % Kelvin
P = [10, 20, 40, 60, 100, 300, 400, 500];
Z_xp = [0.985, 0.970, 0.942, 0.913, 0.869, 0.762, 0.824, 0.910]; 
R = 8.31446261815324e-2; % L.bar/K/mol 

% CO2 properties
Pc = 73.8; % bar
Tc = 30.98 + 273.15; % Kelvin
omega = .225;
% ref: https://www.engineeringtoolbox.com/CO2-carbon-dioxide-properties-d_2017.html
%%%%%%%%%%%%%
%% SRK EoS %%
%%%%%%%%%%%%%
m = .480 + 1.574*omega - .176*omega^2;
a = .42748*(R^2)*(Tc^2)/Pc*(1+ m*(1-sqrt(T/Tc)))^2;
b = .08664*R*Tc/Pc;

% Solving the Equation for Molar Volume (V, ltrs/mole)
for i = 1:length(P)
    syms v
    eqn = P(i)*v^3 - R*T*v^2 - v*(R*T*b + P(i)*b^2 - a) - a*b == 0;
    S = solve(eqn, v, 'Real', true);
    V(i) = double(S); % Ltrs/mole
    Z_SRK(i) = P(i)*V(i)/R/T;
end
SRK_Z_err = sqrt((sum((Z_SRK-Z_xp).^2))/length(Z_xp));

%%%%%%%%%%%%%%%
%% PR-76 EoS %%
%%%%%%%%%%%%%%%
m1 = .375 + 1.574*omega - .267*omega^2;
a1 = .45724*(R^2)*(Tc^2)/Pc*(1+ m1*(1-sqrt(T/Tc)))^2;
b1 = .07780*R*Tc/Pc;

% Solving the Equation for Molar Volume (V, ltrs/mole)
for i = 1:length(P)
    syms x
    eqn1 = -P(i)*x^3 + (R*T-P(i)*b1)*x^2 + (2*b1*R*T - a1 + 3*P(i)*b1^2)*x - P(i)*b1^3 - R*T*b1^2 + a1*b1 == 0;
    S1 = solve(eqn1, x, 'Real', true);
    X(i) = double(S1); % Ltrs/mole
    Z_PR(i) = P(i)*X(i)/R/T;
end

PR_Z_err = sqrt((sum((Z_PR-Z_xp).^2))/length(Z_xp));
%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting the Results %%
%%%%%%%%%%%%%%%%%%%%%%%%%%

% % % SRK % % %

% Z vs. P
figure
plot(P, Z_SRK, 'g--', 'LineWidth', 2);
hold on
plot(P, Z_xp, 'r', 'LineWidth', 2);
xlabel('Pressure, bar');
ylabel('Z-Factor');
legend('Analytical', 'Experimental');
title('Z-Factor vs. Pressure using SRK EoS');
grid on
hold off

% Z vs. exp. Z
figure 
scatter(Z_SRK, Z_xp, 'b*', 'LineWidth', 2);
title('SRK Z-Factor vs. Experimental Results');
hold on
% fitting a line through data
[xData, yData] = prepareCurveData(Z_SRK, Z_xp );
ft = fittype( 'poly1' );
[SRK_fitresult, SRK_gof] = fit( xData, yData, ft );
plot(SRK_fitresult);
hold off
% meta
xlabel('SRK Z-Factor');
ylabel('Experimental Z-Factor');
legend('Data', 'Fitted Curve');
str = ['R2 = ', num2str(SRK_gof.rsquare)];
annotation('textbox',[0.15 0.75 0.3 0.15], 'String', str,'FitBoxToText','on');

% % % PR-76 % % %

% Z vs. P
figure
plot(P, Z_PR, 'g--', 'LineWidth', 2);
hold on
plot(P, Z_xp, 'r', 'LineWidth', 2);
xlabel('Pressure, bar');
ylabel('Z-Factor');
legend('Analytical', 'Experimental');
title('Z-Factor vs. Pressure using Peng-Robinson 1976 EoS');
grid on
hold off

% Z vs. exp. Z
figure 
plot(Z_PR, Z_xp, 'b*', 'LineWidth', 2);
title('Peng Robinson Z-Factor vs. Experimental Results');
hold on
% fitting a line through data
[xD, yD] = prepareCurveData(Z_PR, Z_xp );
ft1 = fittype( 'poly1' );
[PR_fitresult, PR_gof] = fit( xD, yD, ft1 );
plot(PR_fitresult);
hold off
% meta
xlabel('Peng-Robinson Z-Factor');
ylabel('Experimental Z-Factor');
legend('Data', 'Fitted Curve');
str1 = ['R2 = ', num2str(PR_gof.rsquare)];
annotation('textbox',[0.15 0.75 0.3 0.15], 'String', str1,'FitBoxToText','on');

% % % Comparison % % %
for i = 1:length(P)
    V_xp(i) = Z_xp(i)*R*T/P(i);
end
% V vs. P
figure
plot(P, V, 'b', 'LineWidth', 2);
hold on
plot(P, X, 'm--', 'LineWidth', 2);
hold on
plot(P, V_xp, 'c:', 'LineWidth', 2);
xlabel('Pressure, bar');
ylabel('Molar Volume, L/mole');
legend('SRK', 'Peng-Robinson', 'Experimental Data');
title('Comparison of SRK and PR in determining Molar Volume vs. P');
grid on
hold off
