% -------------------------------------------------------------------------
% FallendeKette01.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet Dynamik der fallenden Kette auf Basis der
% Lagrange-Gleichungen  
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];


%% Numerische Lösung

% Parameter

L       = 10.0;             % Länge Kette m
mu      = 5.0;              % Massendichte in kg/m 
g       = 9.81;             % g in m/s 
z0      = 0;                % Anfangsgposition der Kette
dz0     = 0;                % Anfangsgeschwindigkeit der Kette
tmax    = sqrt(L/g);        % Fallzeit freier Fall

fprintf('\n ');
fprintf('\n L   = %8.2f m', L);
fprintf('\n mu  = %8.2f kg/m', mu);
fprintf('\n ');
fprintf('\n z0  = %8.2f m', z0);
fprintf('\n dz0 = %8.2f m/s', dz0);
fprintf('\n g   = %8.2f m/s^2', g);
fprintf('\n tmax= %8.2f m/s', tmax);
fprintf('\n ');


%% Berechnungen

% Anfangswerte
AB=[z0;dz0]; % AB für ode45

% Numerische Lösung Volle Lagrangegleichung
opt=odeset('AbsTol',1.e-7,'RelTol',1.e-6);
% Numerische Lösung LGL
[t,Y]=ode45(@dgl_FallingChain,[0.0,tmax],AB,opt,g,mu,L); 

z  = Y(:,1);
dz = Y(:,2);

zf = -0.5*g*t.^2;
dzf= -g*t;

%% Graphische Ausgabe

% Lösung Lagrangegleichung z(x)
figure();
subplot(1,3,1)
hold on
plot(t,abs(z),'Color',Colors(3,:), 'LineWidth',2);
plot(t,abs(zf),'Color',Colors(4,:), 'LineWidth',2,'LineStyle',Style(3));
axis([0 tmax 0 1.2*L/2]);
grid on
ylabel('z in m','FontSize',14)
xlabel('t in s','FontSize',14)
h=title('Numerische Lösung Lagrange-Gl.');
set(h,'FontSize',14,'FontWeight','normal'); 
set(gca,'FontSize',16);

subplot(1,3,2)
hold on
plot(t,abs(dz),'Color',Colors(3,:), 'LineWidth',2);
plot(t,abs(dzf),'Color',Colors(4,:),'LineWidth',2,'LineStyle', Style(3));
axis([0 tmax 0 1.2*max(abs(dz))]);
axis([0 tmax 0 12]);

grid on
ylabel('Geschwindigkeit in m/s','FontSize',14)
xlabel('t in s','FontSize',14)
set(gca,'FontSize',16);

subplot(1,3,3)
hold on
plot(abs(z),abs(dz),'Color',Colors(3,:), 'LineWidth',2);
plot(abs(zf),abs(dzf),'Color',Colors(4,:),'LineWidth',2,'LineStyle', Style(3));
axis([0 L/2 0 1.2*max(abs(dz))]);
axis([0 L/2 0 12]);
grid on
ylabel('Geschwindigkeit in m/s','FontSize',14)
xlabel('z in m','FontSize',14)
set(gca,'FontSize',16);
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------




%% Funktioen

% Lagrange-Gleichung
function dY = dgl_FallingChain(t,Y,g,mu,L)
    % Y(1)- Position z(t), 
    % Y(2)- Geschwindigkeit dz(t)
    dY    = [Y(2);...    
             -g-0.5*Y(2).^2./(L-Y(1))];
end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------


