% ------------------------------------------------------------------------
% AchterbahnKont.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik des Kontinuums" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet die Bewegung einer Achterbahn in einem
% kontinuierlichen Modell eines elastischen Bandes.
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
Phi     = 100000;                     % Längen-Elastizitätsmodul

NG = 1001;                            % materielle Gitterpunkte
dl = L/(NG-1);                        % Abstand zweier mat. Gitterpunkte
l(:,1) = (0:NG-1)*dl;
                
Tmax = 6*sqrt(2*H/g);                 % maximale Zeit
NT = 100001;                          % Schritte in der Zeit
dt = Tmax/(NT-1);                     % Zeitschritt
t(:)= (0:NT-1)*dt;                    % Zeitgitter

%% Bereite Gitter vor
x = zeros(NT,NG);

% Anfangsbedingung (Zug ist entspannt)
x(1,(1:NG)) = l(:,1)-4/a;
% 2. Zeitschritt für DGL 2. Ordnung notwendig, Bewegung mit 
% geschwindigkeit analog zu diskretem Modell
for i = 1:NG
  dhdx = hp(x(1,i),a,H);
  x(2,i) = x(1,i)+1.018*sqrt(2*g*(H-gauss(x(1,i),a,H))/(1+dhdx^2))*dt;
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
    end

    % Anfang und Ende bewegen sich mit und werden relativ zu den
    % benachbarten Punkten im Abstand dl gesetzt
    x(j,1)  = x(j,2)-dl;
    x(j,NG) = x(j,NG-1)+dl;
end

%% Darstellung der Ergebnisse
xplot = -4/a:8/a/500:4/a;
k     = [10000 20000 30000 50000 60000 90000];  % Zeiten für Plots
pos   = [1 251 501 751 1001];                   % markierte Positionen
figure
for i = 1:6
  subplot(3,2,i);
  plot(xplot(:),gauss(xplot(:),a,H),'linewidth',2, 'color', Colors(3,:));
  hold on
  plot(x(k(i),:),gauss(x(k(i),:),a,H),'linewidth',5, 'color', Colors(1,:));
  for j = 1:5
    plot(x(k(i),pos(j)),gauss(x(k(i),pos(j)),a,H),'o', ...
        'MarkerFaceColor', Colors(2,:), 'MarkerSize', 8);
  end
  hold off
  grid on;
  axis([-4/a 4/a -0.5 1.5*H]);
  xlabel('{\itx}({\itt})');
  ylabel('{\itz}({\itt})');
  titletext = sprintf('t = %.2f s',t(k(i)));
  texto = text(-3.8/a,1.4*H,titletext);
  set(texto, 'FontSize',14,'FontWeight','normal','FontName','Times');
  set(gca,'FontSize',16,'FontName','Times');
end
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