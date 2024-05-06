% -------------------------------------------------------------------------
% VibDaempfung01.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet Dynamik eines Stoßdämpfers 
% auf Basis Lagrange-Gleichungen 
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];


%% Numerik
% Parameter
MF      = 2000;            % Masse Auto in kg
mD      =   10;            % Masse Stoßdämpfer in kg
kF      = 3000;            % Federkonstante Auto
kD      = 8000;            % Federkonstante Stoßdämpfer
etaD    = 5000;            % Dämpfungskoeffizient
tmax    = 10;

fprintf('\n ');
fprintf('\n ');
fprintf('\n Masse Auto :                 MF   = %8.2f kg', MF);
fprintf('\n Masse Stoßdämpfer :          mD   = %8.2f kg', mD);
fprintf('\n Federkonstante Auto :        kF   = %8.2f kg/s² (N/m)', kF);
fprintf('\n Federkonstante Stoßdämpfer : kD   = %8.2f kg/s² (N/m)', kD);
fprintf('\n Dämpfung Stoßdämpfer :       etaD = %8.2f kg/s  (N/(m/s)', etaD);
fprintf('\n ');

% Parameter in LGL
% Anfangswerte
zF0  = 0.1;
dzF0 = -0.1;
zD0  = 0;
dzD0 = 0;

%% Berechnungen ODE45

% Anfangswerte
AB = [zF0;dzF0;zD0;dzD0]; % AB für ode45

P1.MF    = MF;
P1.mD    = mD;
P1.kD    = kD;
P1.kF    = kF;
P1.etaD  = etaD;


% Numerische Lösung Lagrangegleichung
opt=odeset('AbsTol',1.e-7,'RelTol',1.e-7,'events',@MyEvent);
% Numerische Lösung LGL mit Reibung
[t, Y, TE, YE, IE] = ode45(@dgl_VibDamp01,[0.0,tmax],AB,opt,P1); 
zF  = Y(:,1);
dzF = Y(:,2); 
zD  = Y(:,3);
kend = length(zF);


%% 
% Graphische Ausgabe

% Zeitentwicklung 
figure();
hold on
p(1) = plot(t, zF,'Color',Colors(2,:), 'LineWidth',2);
p(2) = plot(t, zD,'Color',Colors(3,:), 'LineWidth',2);
% p(3) = plot(t, dzF,'Color',Colors(4,:), 'LineWidth',2);
grid on
ylabel('zF in m','FontSize',14)
xlabel('t in s','FontSize',14)
h=title('Stossdämpfer ');
set(h,'FontSize',12,'FontWeight','normal'); 
legend(p,'z_F ','z_D', 'location','northeast',...
        'NumColumns',1);
legend box off
set(gca,'FontSize',16);
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------



%% Funktionen

% Lagrangegleichung
function dY = dgl_VibDamp01(t,Y,P1)
    % Y(1)- Position zF(t), 
    % Y(2)- Geschwindigkeit dotzF(t)
    % Y(3)- Position zD(t), 
    % Y(4)- Geschwindigkeit dotzD(t)
    zF  = Y(1);
    dzF = Y(2);
    zD  = Y(3);
    dzD = Y(4);
    MF    = P1.MF;
    mD    = P1.mD;
    kD    = P1.kD;
    kF    = P1.kF;
    etaD  = P1.etaD; 
    dY = [dzF;...
      -kF/MF*zF - kD/MF*(zF-zD);...
      dzD;...
      -etaD/mD*dzD + kD/mD*(zF-zD)];
end

function [value,isterminal,direction] = MyEvent(t,Y,P1)
    value = (0.0-(Y(2)*Y(2)+Y(4)*Y(4))); % kinetische Energie = 0
    isterminal = 1;                      % stop the integration
    direction = 0;                       % negative direction
end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------


