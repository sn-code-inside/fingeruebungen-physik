% -------------------------------------------------------------------------
% SchwingendeSaite03.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik des Kontinuums" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet die Schwingung einer kontinuierlichen Saite, die an
% einem Ende zu Schwingungen angeregt wird.
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

%% Parameter
L       = 0.1;                        % Länge der Saite in m
rho     = 1150;                       % Massen(längen)dichte in kg/m
Phi     = 3e9;                        % Elastizitäts(längen)modul

NG = 1001;                            % Gitterpunkte im Ort
dx = L/(NG-1);                        % Abstand zweier Gitterpunkte
x(:) = (0:NG-1)*dx;
mNG = (NG+1)/2;                       % Mitte der Ortsachse                
                
Tmax = 0.0001;                        % maximale Zeit
NT = 100001;                          % Schritte in der Zeit
dt = Tmax/(NT-1);                     % Zeitschritt
t(:)= (0:NT-1)*dt;                    % Zeitgitter

%% Bereite Gitter vor
y = zeros(NT,NG);

% Anregung
y((1:NT),1) = sin(10*pi/Tmax*t(:));

%% Lösen der DGL durch Methode der finiten Differenzen
for j = 3:NT
    for i = 2:NG-1
        yx2 = y(j-1,i-1)-2*y(j-1,i)+y(j-1,i+1);
        y(j,i) = 2*y(j-1,i)-y(j-2,i)+Phi/rho*(dt/dx)^2*yx2;    
    end                   
end

%% Darstellung der Ergebnisse
k = [1 20000 60000 65000 80000 100001];
figure
for i = 1:6
  subplot(3,2,i);
  plot(x,y(k(i),:),'linewidth',2, 'color', Colors(1,:));
  grid on;
  axis([min(x) max(x) -2 2]);
  xlabel('\itx');
  ylabel('{\ity}({\itt})');
  titletext = sprintf('t = %.2f ms',t(k(i))*1000);
  texto = text(L/4,1.5,titletext);
  set(texto, 'FontSize',14,'FontWeight','normal','FontName','Times');
  set(gca,'FontSize',16,'FontName','Times');
end
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------

