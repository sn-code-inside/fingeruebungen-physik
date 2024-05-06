% ------------------------------------------------------------------------
% AchterbahnKontBeschl.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik des Kontinuums" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet die Bewegung einer Achterbahn in einem
% kontinuierlichen Modell eines elastischen Bandes. Die dabei wirkenden
% Beschleunigungen auf einen Fahrgast werden berechnet.
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

%% Parameter
H      = 8;                           % Höhe Berg
B      = 10;                          % Breite Berg
a      = sqrt(log(2))/B;              % Gauss-Parameter 
g      = 9.81;                        % Schwerebeschleunigung
 
L       = 10.;                        % Länge des Zuges in m
rho     = 500;                        % Massenlängendichte in kg/m
Phi     = 200000;                     % Längen-Elastizitätsmodul

NG = 1001;                            % materielle Gitterpunkte
dl = L/(NG-1);                        % Abstand zweier mat. Gitterpunkte
l(:,1) = (0:NG-1)*dl;
                
Tmax = 8*sqrt(2*H/g);                 % maximale Zeit
NT = 100001;                          % Schritte in der Zeit
dt = Tmax/(NT-1);                     % Zeitschritt
t(:)= (0:NT-1)*dt;                    % Zeitgitter

%% Bereite Gitter vor
% Positionen des Zuges
x = zeros(NT,NG);
% Beschleunigungen der Personen im Zug in x-Richtung
apx = zeros(NT,NG);
% z-Richtung
apz = zeros(NT,NG);
% Normalkomponenten
apn = zeros(NT,NG);

% Anfangsbedingung (Zug ist entspannt)
x(1,(1:NG)) = l(:,1)-4/a;
% 2. Zeitschritt für DGL 2. Ordnung notwendig, Bewegung mit 
% geschwindigkeit analog zu diskretem Modell
for i = 1:NG
  dhdx = hp(x(1,i),a,H);
  x(2,i) = x(1,i)+0.98*sqrt(2*g*(H-gauss(x(1,i),a,H))/(1+dhdx^2))*dt;
end

%% Lösen der DGL mit der Finite-Differenzen-Methode
for j = 3:NT
    for i = 2:NG-1
        dhdx   = hp(x(j-1,i),a,H);                 % dh/dx
        d2hdx2 = h2p(x(j-1,i),a,H);                % d^2h/dx^2
        b      = dhdx/(1+dhdx^2);
        c      = d2hdx2*b;

        dxdt   = x(j-1,i)-x(j-2,i);                % dx/dt *dt
        dxdl   = x(j-1,i)-x(j-1,i-1);              % dx/dl *dl
        d2xdl2 = x(j-1,i+1)+x(j-1,i-1)-2*x(j-1,i); % d^2x/dl^2 *dl^2
        x(j,i) = 2*x(j-1,i)-x(j-2,i) - c*dxdt^2 - g*b*dt^2 ...
                 + Phi/rho*(dt/dl)^2*( c*dxdl^2 + d2xdl2 );

        % Schreiben der Beschleunigungen, die eine Person erfährt,
        % so lange sie nicht gegen die Halterungen des Wagens stößt
        
        % Normalkomponente
        apn(j,i) = (d2hdx2*(dxdt/dt)^2 - g)/sqrt(1+dhdx^2);
        
        % x- und z-Komponenten
        apx(j,i) = -c*(dxdt/dt)^2;
        apz(j,i) = -apx(j,i)/dhdx-g;
    end

    % Anfang und Ende bewegen sich mit und werden relativ zu den
    % benachbarten Punkten im Abstand dl gesetzt
    x(j,1)  = x(j,2)-dl;
    x(j,NG) = x(j,NG-1)+dl;
end

%% Darstellung der Ergebnisse
figure
subplot(3,1,1);
plot(x(:,998),apn(:,998)/g,'linewidth',4, 'color', Colors(4,:), ...
    'LineStyle','--');
hold on;
plot(x(:,998),apx(:,998)/g,'linewidth',2, 'color', Colors(2,:));
hold on;
plot(x(:,998),apz(:,998)/g,'linewidth',2, 'color', Colors(3,:));
hold off;
grid on;
axis([-3/a 3/a -3 0.5]);
legend('Normalkomponenten', 'a_x', 'a_z', 'location','best', ...
    'numcolumns',1);
xlabel('\itx');
ylabel('\ita/g');
titletext = "Fahrgast am vorderen Ende des Zuges";
texto = text(-34,-2.5,titletext);
set(texto, 'FontSize',14, 'FontWeight' ,'normal','FontName','Times');
set(gca,'FontSize',16,'FontName','Times');

subplot(3,1,2);
plot(x(:,501),apn(:,501)/g,'linewidth',4, 'color', Colors(4,:), ...
    'LineStyle','--');
hold on;
plot(x(:,501),apx(:,501)/g,'linewidth',2, 'color', Colors(2,:));
hold on;
plot(x(:,501),apz(:,501)/g,'linewidth',2, 'color', Colors(3,:));
hold off;
grid on;
axis([-3/a 3/a -3 0.5]);
xlabel('\itx');
ylabel('\ita/g');
titletext = "Fahrgast in der Mitte des Zuges";
texto = text(-34,-2.5,titletext);
set(texto, 'FontSize',14, 'FontWeight' ,'normal','FontName','Times');
set(gca,'FontSize',16,'FontName','Times');

subplot(3,1,3);
plot(x(:,3),apn(:,3)/g,'linewidth',4, 'color', Colors(4,:), ...
    'LineStyle','--');
hold on;
plot(x(:,3),apx(:,3)/g,'linewidth',2, 'color', Colors(2,:));
hold on;
plot(x(:,3),apz(:,3)/g,'linewidth',2, 'color', Colors(3,:));
hold off;
grid on;
axis([-3/a 3/a -3 0.5]);
xlabel('\itx');
ylabel('\ita/g');
titletext = "Fahrgast am hinteren Ende des Zuges";
texto = text(-34,-2.5,titletext);
set(texto, 'FontSize',14, 'FontWeight' ,'normal','FontName','Times');
set(gca,'FontSize',16,'FontName','Times');
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


%%
% Funktionen

function zg = gauss(x,a,H)
    zg = H*exp(-a^2*x.^2);
end

function hprime = hp(x,a,H)
    hprime = -2*a^2*x.*gauss(x,a,H);
end

function h2prime = h2p(x,a,H)
    h2prime = 2*a^2*(2*a^2*x.^2-1).*gauss(x,a,H);
end

% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------
