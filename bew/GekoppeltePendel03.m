% -------------------------------------------------------------------------
% GekoppeltePendel03.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Verwendet Programmteile von
% Lagrange Mechanics  Version 1.4.0.0 (4.48 KB) 
% von Peter van Alem
% Simulation of double and coupled pendulum
% -------------------------------------------------------------------------
% Dynamik der gekoppelten Federpendel 
% Vergleich analytische und numerische Lösung für unterschiedliche Pendel
% große Winkel
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

%%
% Generalisierte Koordinaten
% Achtung !!!
% Gegenüber der Aufgabenstellung sind die generalisierten Koordinaten hier
% zwecks einfacherer Lesbarkeit

% theta1 ---> theta
% theta2 ---> phi

%%


syms t dum_
theta = str2sym('theta(t)');
phi = str2sym('phi(t)');

% Konstanten für Pendel unterschiedlicher Masse und Länge, aber kleine
% Winkel
L1 = 6;
L2 = 4;
m1 = 5;
m2 = 3;
LF = 4;
g   = 9.81;
c1   = LF/L1;
c2   = LF/L2;
keff = 0.5; %Federkonstante
kappa1 = keff*c1^2/m1;
kappa2 = keff*c2^2/m2;

d_J1_J2 = 10; % Abstand Aufhängung 
d_0 = 10;     % Restlänge der Feder

% Positionen und Geschwindigkeiten
x1 = -d_J1_J2/2 + L1 * sin(theta);
y1 = -L1 * cos(theta);
x2 = d_J1_J2/2 + L2 * sin(phi); 
y2 = -L2 * cos(phi); 
x1_dot = diff(x1, t);
x2_dot = diff(x2, t);
y1_dot = diff(y1, t);
y2_dot = diff(y2, t);

% Kinetische Energie 
T = m1/2 * (x1_dot^2 + y1_dot^2) + m2/2 * (x2_dot^2 + y2_dot^2);

% Federbefestigung
xf1 = -d_J1_J2/2 + L1 * sin(theta)*c1;
yf1 = y1*c1;
xf2 = d_J1_J2/2 + L2 * sin(phi)*c2;
yf2 = y2*c2;

% Potentielle Energie
V = m1 * g * y1 + m2 * g * y2 +...
    1/2 * keff * (sqrt((xf2 - xf1)^2 + (yf2 - yf1)^2) - d_0)^2;

% Lagrange-Funktion
L = T - V;

% dL/d(qdot)
dL_dthetadot = subs(diff(subs(L, diff(theta, t), dum_), dum_), dum_, diff(theta, t));
dL_dphidot = subs(diff(subs(L, diff(phi, t), dum_), dum_), dum_, diff(phi, t));
% dL/dq
dL_dtheta = subs(diff(subs(L, theta, dum_), dum_), dum_, theta);
dL_dphi = subs(diff(subs(L, phi, dum_), dum_), dum_, phi);

% dFdq Reibungsfunktion
eta = 0.000001; % Reibungskonstante
F = 1/2 * eta * (x1_dot^2 + y1_dot^2 + x2_dot^2 + y2_dot^2);
dF_dthetadot = subs(diff(subs(F, diff(theta, t), dum_), dum_), dum_, diff(theta, t));
dF_dphidot = subs(diff(subs(F, diff(phi, t), dum_), dum_), dum_, diff(phi, t));

% Lagrange-Gleichungen
deq_1 = diff(dL_dthetadot, t) - dL_dtheta + dF_dthetadot;  %Pendel1
deq_2 = diff(dL_dphidot, t) - dL_dphi + dF_dphidot;        %Pendel2

% Variablen
variables = {theta, phi, diff(theta, t), diff(phi, t), diff(theta, t, 2), diff(phi, t, 2)};    
variables_short = arrayfun(@str2sym, {'x(1)', 'x(2)', 'x(3)', 'x(4)', 'thetaddot', 'phiddot'});
deq_1 = subs(deq_1, variables, variables_short);
deq_2 = subs(deq_2, variables, variables_short);

% Numerische Lösung für thetaddot, phiddot
solution = solve(deq_1, deq_2, str2sym('thetaddot'), str2sym('phiddot'));
THETADDOT = solution.thetaddot;
PHIDDOT = solution.phiddot;

% ODE DGL-System
NPts = 1000;
tend = 50;
time = linspace(0, tend, NPts);
% AB [theta, phi, thetadot, phidot]
theta10 = deg2rad(-60);
x_0  = [theta10 0  0 0];

str = ['x_dot = @(t, x)[x(3); x(4);', char(THETADDOT), ';', char(PHIDDOT), '];'];
eval(str);
[t, q] = ode45(x_dot, time, x_0);

% Position  Massen kartesische Koordinaten
X1 = -d_J1_J2/2 + L1 * sin(q(:, 1));
Y1 = -L1 * cos(q(:, 1));
X2 = d_J1_J2/2 + L2 * sin(q(:, 2));
Y2 = -L2 * cos(q(:, 2));

% Feder in kartesischen Korodinaten
XF1 = -d_J1_J2/2 + L1 * sin(q(:, 1))*c1;
YF1 = Y1*c1;
XF2 = d_J1_J2/2 + L2 * sin(q(:, 2))*c2;
YF2 = Y2*c2;


%% Graphik

% Simulation
set(gcf, 'color', 'w')
set(gcf, 'position', [100, 100, 750, 750])
h = plot([]);
hold on
box on
grid on;
axis equal
for ie = 1 : numel(time)/5
    i = 5*ie;
    if ~ishghandle(h)
        break
    end
    cla
    plot([-d_J1_J2/2, X1(i)],[0, Y1(i)],'Color',Colors(2,:),'Linewidth',2);
    plot(X1(i), Y1(i), 'o', 'MarkerFaceColor', Colors(2,:), ...
    'MarkerEdgeColor', Colors(2,:), 'MarkerSize', 5 * m1);
    plot([d_J1_J2/2, X2(i)], [0, Y2(i)],'Color',Colors(4,:),'Linewidth',2);
    plot(X2(i), Y2(i), 'o', 'MarkerFaceColor',Colors(4,:) ,...
    'MarkerEdgeColor', Colors(4,:), 'MarkerSize', 5 * m2);
    axis([-12, 12, -10, 5]);
    line([-d_J1_J2, +d_J1_J2],[0,0],'LineWidth',2);
    h = draw_spring_2D([XF1(i); YF1(i)], [XF2(i); YF2(i)], 20, 0.6);
    drawnow
end

% Trajektorien numerisch
figure()
hold on
plot(t,rad2deg(q(:, 1)),'Color', Colors(2,:),'LineWidth',2,'LineStyle',Style(1));
plot(t,rad2deg(q(:, 2)),'Color', Colors(4,:),'LineWidth',2,'LineStyle',Style(2));
grid on;
legend('m_1','m_2');
legend box off
ylabel('\theta_1, \theta_2','FontSize',13), xlabel('\it t \rm in s','FontSize',13)
legend box off
axis([0, tend, -rad2deg(abs(x_0(1))+abs(x_0(2))),...
      rad2deg(abs(x_0(1))+abs(x_0(2)))]); 
set(gca,'FontSize',14)

% Trajektorien analytisch kleine Winkel
omegaI  = sqrt(g/L1);
omegaII = sqrt(g/L2);
omega1 = sqrt((omegaI^2+omegaII^2+kappa1+kappa2)/2-0.5*...
               sqrt((omegaI^2-omegaII^2)^2+(kappa1+kappa2)^2+ ...
                 2*(omegaI^2-omegaII^2)*(kappa1-kappa2)));
omega2 = sqrt((omegaI^2+omegaII^2+kappa1+kappa2)/2+0.5*...
               sqrt((omegaI^2-omegaII^2)^2+(kappa1+kappa2)^2+ ...
                 2*(omegaI^2-omegaII^2)*(kappa1-kappa2)));
C1 = (omegaII^2-omega1^2+kappa2)/kappa2;
C2 = (omegaII^2-omega2^2+kappa2)/kappa2;

theta1 = rad2deg(theta10*(C1*cos(omega1*time)-C2*cos(omega2*time))/(C1-C2));
theta2 = rad2deg(theta10*(cos(omega1*time)-cos(omega2*time))/(C1-C2));
figure()
hold on
plot(time,theta1,'Color', Colors(2,:),'LineWidth',2,'LineStyle',Style(1));
plot(time,theta2,'Color', Colors(4,:),'LineWidth',2,'LineStyle',Style(2));
grid on;
legend('m_1','m_2');
legend box off
ylabel('\theta_1, \theta_2','FontSize',13), xlabel('\it t \rm in s','FontSize',13)
legend box off
axis([0, tend, -rad2deg(abs(x_0(1))+abs(x_0(2))),...
      rad2deg(abs(x_0(1))+abs(x_0(2)))]); 
set(gca,'FontSize',14)
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


%%
% Federsimulation 
% Benützt to Lagrange Mechanics  Version 1.4.0.0 (4.48 KB) 
% von Peter van Alem
% Simulation of double and coupled pendulum

function h = draw_spring_2D(A, B, number_of_coils, y_amplitude)    
    persistent t
    
    normalvector_AB = (B - A) / norm(B - A);
    offset_A = A + 1.25 * normalvector_AB;
    offset_B = B - 1.25 * normalvector_AB;
    distance_between_offsets = norm(offset_B - offset_A);
    
    t = linspace(-pi, number_of_coils * 2 * pi, 500);
    x_coordinate_between_offsets = distance_between_offsets * linspace(0, 1, numel(t));
    
    ratio_X_div_Y = 0.5;
    
    x = x_coordinate_between_offsets + ratio_X_div_Y * y_amplitude * cos(t);
    y = y_amplitude * sin(t);
    
    coil_positions = [x; y];

    rotation_matrix = [normalvector_AB, null(normalvector_AB')];
    rotated_coil_positions = rotation_matrix * coil_positions;
    h = plot([A(1), offset_A(1) + rotated_coil_positions(1,:), B(1)], ...
         [A(2), offset_A(2) + rotated_coil_positions(2,:), B(2)], 'k');
end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------
