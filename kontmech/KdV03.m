% -------------------------------------------------------------------------
% KdV03.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik des Kontinuums" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm propagiert eine Summe aus Solitonien der Korteweg-De-Fries-
% Gleichung und zeigt exemplarisch, dass diese Summe kein Soliton
% darstellt bzw. dies näherungsweise erfüllt ist, wenn die Solitonen sich
% nicht überschneiden (Ergänzung zur analytischen Rechnung).
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

%% Parameter
L = 80;                               % Länge der Wasserlinie

NG = 801;                             % Gitterpunkte im Ort
dx = L/(NG-1);                        % Abstand zweier Gitterpunkte
x(:,1) = (0:NG-1)*dx;
                
Tmax = 20.;                           % maximale Zeit
NT = 100001;                          % Schritte in der Zeit
dt = Tmax/(NT-1);                     % Zeitschritt
t(:)= (0:NT-1)*dt;                    % Zeitgitter

%% Bereite Gitter vor
ys1 = zeros(NT,NG);
ys2 = zeros(NT,NG);

%% Propagation für eine Summe aus zwei sich überschneidenden Solitonen
% Anfangsbedingung: Solitonen mit Geschwindigkeiten c1 und c2
c1 = 1.;
c2 = 1.5;
ys1(1,(1:NG)) = 0.5*c1*sech(0.5*sqrt(c1)*(x(:,1)-10)).^2 ...
    + 0.5*c2*sech(0.5*sqrt(c2)*(x(:,1)-12)).^2;
% 2. Zeitschritt für symmetrische Form der diskreten Ableitung notwendig, 
% Soliton ist mit Geschwindigkeit c um dt vorangeschritten
ys1(2,(1:NG)) = 0.5*c1*sech(0.5*sqrt(c1)*(x(:,1)-10-c1*dt)).^2 ...
    + 0.5*c2*sech(0.5*sqrt(c2)*(x(:,1)-12-c2*dt)).^2;

% Lösen der DGL mit der Finite-Differenzen-Methode
for j = 3:NT
    for i = 3:NG-2
        yx  = ys1(j-1,i+1)-ys1(j-1,i-1);
        yx3 = ys1(j-1,i+2)-2*ys1(j-1,i+1)+2*ys1(j-1,i-1)-ys1(j-1,i-2);
        ys1(j,i) = ys1(j-2,i)-(dt/dx)*(6*ys1(j-1,i)*yx+yx3/dx^2);    
    end                   
end

% Darstellung der Ergebnisse
k = [1 20000 30000 50000 70000 100001];
figure
for i = 1:6
  subplot(3,2,i);
  plot(x,ys1(k(i),:),'linewidth',3, 'color', Colors(4,:));
  grid on;
  axis([min(x) max(x) -0.1 2]);
  xlabel('{\itx}');
  ylabel('{\it\psi}({\itt})');
  titletext = sprintf('{\\itt} = %.1f',t(k(i)));
  texto = text(L/4,1.8,titletext);
  set(texto,'FontSize',16,'FontWeight','normal','FontName','Times');
  set(gca,'FontSize',16,'FontName','Times');
end

%% Propagation für eine Summe aus zwei sich Solitonen, die sich nicht
%% überschneiden
ys2(1,(1:NG)) = 0.5*c1*sech(0.5*sqrt(c1)*(x(:,1)-10)).^2 ...
    + 0.5*c2*sech(0.5*sqrt(c2)*(x(:,1)-20)).^2;
ys2(2,(1:NG)) = 0.5*c1*sech(0.5*sqrt(c1)*(x(:,1)-10-c1*dt)).^2 ...
    + 0.5*c2*sech(0.5*sqrt(c2)*(x(:,1)-20-c2*dt)).^2;

% Lösen der DGL mit der Finite-Differenzen-Methode
for j = 3:NT
    for i = 3:NG-2
        yx  = ys2(j-1,i+1)-ys2(j-1,i-1);
        yx3 = ys2(j-1,i+2)-2*ys2(j-1,i+1)+2*ys2(j-1,i-1)-ys2(j-1,i-2);
        ys2(j,i) = ys2(j-2,i)-(dt/dx)*(6*ys2(j-1,i)*yx+yx3/dx^2);    
    end                   
end

% Darstellung der Ergebnisse
k = [1 20000 30000 50000 70000 100001];
figure
for i = 1:6
  subplot(3,2,i);
  plot(x,ys2(k(i),:),'linewidth',3, 'color', Colors(3,:));
  grid on;
  axis([min(x) max(x) -0.1 2]);
  xlabel('{\itx}');
  ylabel('{\it\psi}({\itt})');
  titletext = sprintf('{\\itt} = %.1f',t(k(i)));
  texto = text(L/4,1.8,titletext);
  set(texto,'FontSize',16,'FontWeight','normal','FontName','Times');
  set(gca,'FontSize',16,'FontName','Times');
end

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------

