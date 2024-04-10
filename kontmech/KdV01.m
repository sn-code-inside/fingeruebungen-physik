% -------------------------------------------------------------------------
% KdV01.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik des Kontinuums" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm propagiert ein Soliton der Korteweg-De-Fries-Gleichung.
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

%% Parameter
L = 40;                               % Länge der Wasserlinie

NG = 401;                             % Gitterpunkte im Ort
dx = L/(NG-1);                        % Abstand zweier Gitterpunkte
x(:,1) = (0:NG-1)*dx;
                
Tmax = 20.;                           % maximale Zeit
NT = 100001;                          % Schritte in der Zeit
dt = Tmax/(NT-1);                     % Zeitschritt
t(:)= (0:NT-1)*dt;                    % Zeitgitter

%% Bereite Gitter vor
y = zeros(NT,NG);

% Anfangsbedingung: Soliton mit Geschwindigkeit c
c = 1.;
y(1,(1:NG)) = KdVSoliton(x,0,c,10,NG);
% 2. Zeitschritt für symmetrische Form der diskreten Ableitung notwendig, 
% Soliton ist mit Geschwindigkeit c um dt vorangeschritten
y(2,(1:NG)) = KdVSoliton(x,dt,c,10,NG);

%% Lösen der DGL mit der Finite-Differenzen-Methode
for j = 3:NT
    for i = 3:NG-2
        yx  = y(j-1,i+1)-y(j-1,i-1);
        yx3 = y(j-1,i+2)-2*y(j-1,i+1)+2*y(j-1,i-1)-y(j-1,i-2);
        y(j,i) = y(j-2,i)-(dt/dx)*(6*y(j-1,i)*yx+yx3/dx^2);    
    end                   
end

%% Darstellung der Ergebnisse
k = [1 20000 30000 50000 70000 100001];
figure
for i = 1:6
  subplot(3,2,i);
  plot(x,KdVSoliton(x,dt*(k(i)-1),c,10,NG),'linewidth',2, 'color', ...
      Colors(3,:));
  hold on;
  plot(x,y(k(i),:),'linewidth',5, 'color', Colors(4,:),'LineStyle','--');
  grid on;
  axis([min(x) max(x) -1.2 2]);
  xlabel('{\itx}');
  ylabel('{\it\psi}({\itt})');
  titletext = sprintf('{\\itt} = %.1f',t(k(i)));
  texto = text(L/4,1.5,titletext);
  set(texto,'FontSize',16,'FontWeight','normal','FontName','Times');
  set(gca,'FontSize',16,'FontName','Times');
end
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


%% Funktion
% analytische Lösung eines KdV-Solitions zu einem festen Zeitpunkt t
function  y  =  KdVSoliton(x,t,c,x0,NG)
    y(1:NG) = 0.5*c*sech(0.5*sqrt(c)*(x(1:NG,1)-x0-c*t)).^2;
end