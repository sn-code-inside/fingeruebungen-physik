% -------------------------------------------------------------------------
% Doppelkegel03.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Aufwärtsrollender Doppelkegel
% 
% Programm berechnet Lösungen den Phasenraum aus den
% Lagrange-Gleichungen des aufwärtsrollenden Doppelkegels und die
% Energiekonversion
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

H       = 6.5;              % Länge Schaukel cm
R       = 3.0;              % Länge Oberkörper cm 
phi     = deg2rad(15.3);    % Öffnungswinkel Schienen
theta   = deg2rad(06.5);    % Anstiegswinkel Schienen
psi     = atan(R/H);        % Öffnungswinkel Kegel
M       = 122.5;            % Masse Doppelkegel in g
J       = 3*M*R^2/10;       % Trägheitsmoment Doppelkegel in g cm²
g       = 981;              % g in cm/s

alpha   = asin(tan(psi)*tan(phi));  %Richtungswinkel alpha
q0      = 0;             % Anfangsgposition gen. KO des Doppelkegels
dq0     = 0;             % Anfangsgeschwindigkeit gen. KO des Doppelkegels
ys0     = H*tan(psi)*sin(alpha);
                         % Anfangsgposition ys im KOS Sigma


fprintf('\n ');
fprintf('\n phi   = %4.2f°', rad2deg(phi));
fprintf('\n psi   = %4.2f°', rad2deg(psi));
fprintf('\n theta = % 4.2f°', rad2deg(theta));
fprintf('\n alpha = % 4.2f°', rad2deg(alpha));
fprintf('\n ');
fprintf('\n H  = %8.2f cm', H);
fprintf('\n R  = %8.2f cm', R);
fprintf('\n M  = %8.2f g', M);
fprintf('\n J  = %8.2f gcm²', J);
fprintf('\n ');
fprintf('\n ys0= %8.2f cm', ys0);
fprintf('\n g  = %8.2f cm/s', g);
fprintf('\n ');
fprintf('\n Rollbedingung erfüllt ??? ');
b1 = tan(theta);
b2 = tan(phi)*tan(psi)/sqrt(1-(tan(phi)*tan(psi))^2);
if b1 < b2
    fprintf('\n %s  < %s ! Rollbedingung erfüllt!!!',num2str(b1,4),...
        num2str(b2,4));
    fprintf('\n ');
else
    fprintf('\n %s  > %s ! Rollbedingung nicht(!) erfüllt! ',...
        num2str(b1,4), num2str(b2,4));
    fprintf('\n ');
end
fprintf('\n ');



%% Berechnungen

% Anfangswerte
AB=[q0;dq0]; % AB für ode45
tmax  = 25;  % maximale Berechnungszeit in s

% Parameterset für ODE45
P1.H =H;
P1.R =H;
P1.alpha = alpha;
P1.psi   = psi;
P1.phi   = phi;
P1.C1 = 3*R*R/10/(tan(alpha))^2;
P1.C2 = -g*sin(alpha-theta);
P1.C3 = H*tan(psi)/tan(alpha);


%% Berechnungen

qmax  = 0.9999*H*tan(psi)/tan(alpha);
q = linspace(0,qmax,100);

E0      = 0;
% dq = sqrt(2*(E0+M*g*q*sin(alpha-theta))./...
%      (M+(J./(H*tan(psi)-q*tan(alpha)).^2)));
dq = sqrt(-2*P1.C2*q./(1+P1.C1./(P1.C3-q).^2));

gamma   = H*tan(psi)/tan(alpha);
U       = E0+M*g*q*sin(alpha-theta);
T_trans = M*(tan(alpha))^2.*(gamma-q).^2.*(E0+M*g*q*sin(alpha-theta))./...
          (M*(tan(alpha))^2.*(gamma-q).^2 +J);
T_rot   = J*(E0+M*g*q*sin(alpha-theta))./(M*(tan(alpha))^2.*(gamma-q).^2+J);

%% Graphische Ausgabe
% Lösung Phasenraum und
% Vergleich mit numerischer Berechnung aus Doppelkegel02.m (LGL)

figure();
hold on
plot(q,dq,'Color',Colors(3,:), 'LineWidth',2);
line([qmax qmax],[0, 1.2*max(dq)],'Color',Colors(4,:), 'LineWidth',1);
axis([0 1.2*qmax 0 1.2*max(dq)]);
grid on
xlabel('q in cm','FontSize',14)
ylabel('dq in cm/s','FontSize',14)
h=title('Phasenraum');
set(h,'FontSize',12,'FontWeight','normal'); 
set(gca,'FontSize',16);

% Energiegleichung
figure();
hold on
plot(q,T_trans/1e5,'Color',Colors(3,:), 'LineWidth',2);
plot(q,T_rot/1e5,'Color',Colors(4,:), 'LineWidth',2);
plot(q,U/1e5,'Color',Colors(2,:), 'LineWidth',2);
line([qmax qmax],[0, max(E0+U)/1e5],'Color',Colors(4,:), 'LineWidth',1);
axis([0 1.2*qmax 0 max(E0+U)/1e5]);
grid on
xlabel('q in cm','FontSize',14)
ylabel('E(q) in Ws ','FontSize',14)
legend('T_{trans}', 'T_{rot}', 'U', 'location','northwest');
legend box off
h=title('Energieverteilung');
set(h,'FontSize',12,'FontWeight','normal');
set(gca,'FontSize',16);
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


