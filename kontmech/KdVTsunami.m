% -------------------------------------------------------------------------
% KdVTsunami.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik des Kontinuums" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm propagiert ein Soliton der Korteweg-De-Fries-Gleichung auf eine
% ansteigende Küste hin.
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

%% Parameter
L = 500;                              % Länge der Wasserlinie

NG = 201;                             % Gitterpunkte im Ort
dx = L/(NG-1);                        % Abstand zweier Gitterpunkte
x(:,1) = (0:NG-1)*dx;
                
Tmax = 220.;                          % maximale Zeit
NT = 10001;                           % Schritte in der Zeit
dt = Tmax/(NT-1);                     % Zeitschritt
t(:)= (0:NT-1)*dt;                    % Zeitgitter

g = 9.91;                             % Schwerebeschleunigung in m/s^2
s = 0.07495;                          % Oberflächenspannung Wasser, 5°C
                                      % in N/m
rho = 1000.;                          % Dichte Wasser, 5°C in kg/m^3

%% Bereite Gitter vor
y = zeros(NT,NG);

% Wassertiefe
d(1:(NG-1)/2) = 10;
d((NG-1)/2+1:NG) = 10*(1-(0:(NG-1)/2)/((NG-1)/2));
% KdV-Paramter
sigma(:) = d(:).^3/3.-s.*d(:)./rho./g;

% Anfangsbedingung: Soliton mit Geschwindigkeit c
c = 1.;
y(1,(1:NG)) = 2*c*sqrt(d(1)/g) ...
   *sech(0.5*sqrt(2*c*sqrt(d(1)/g)/sigma(1))*(x(:,1)-80)).^2;
% 2. Zeitschritt für symmetrische Form der diskreten Ableitung notwendig, 
% Soliton ist mit Geschwindigkeit c um dt vorangeschritten
y(2,(1:NG)) = 2*c*sqrt(d(1)/g) ...
   *sech(0.5*sqrt(2*c*sqrt(d(1)/g)/sigma(1))*(x(:,1)-80-c*dt)).^2;

%% Lösen der DGL mit der Finite-Differenzen-Methode
for j = 3:NT
    for i = 3:NG-2
        yx  = y(j-1,i+1)-y(j-1,i-1);
        yx3 = y(j-1,i+2)-2*y(j-1,i+1)+2*y(j-1,i-1)-y(j-1,i-2);
        y(j,i) = y(j-2,i)-(dt/dx)*0.5*sqrt(g/d(i))*(3*y(j-1,i)*yx ...
            + sigma(i)*yx3/dx^2); 
    end                   
end

%% Darstellung der Ergebnisse
k = [1 7000 8000 9000 9500 10001];
figure
for i = 1:6
  subplot(3,2,i);
  plot(x,y(k(i),:),'linewidth',3, 'color', Colors(3,:));
  hold on
  plot(x,-d,'linewidth',2, 'color', Colors(5,:));
  grid on;
  axis([min(x) max(x) -12 6]);
  xlabel('{\itx} in m');
  ylabel('{\ity}({\itt}) in m');
  titletext = sprintf('t = %.1f s',t(k(i)));
  texto = text(L/4,4,titletext);
  set(texto,'FontSize',16,'FontWeight','normal','FontName','Times');
  set(gca,'FontSize',16,'FontName','Times');
end

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------

