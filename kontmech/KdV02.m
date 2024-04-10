% -------------------------------------------------------------------------
% KdV02.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik des Kontinuums" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm propagiert Solitonien der Korteweg-De-Fries-Gleichung und
% vergleicht verschiedene Solitonen für verschiedene Parameter c.
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
y1 = zeros(NT,NG);
y2 = zeros(NT,NG);
y3 = zeros(NT,NG);
y4 = zeros(NT,NG);

%% Lösen der DGL für die erste Geschwindigkeit
% Anfangsbedingung: Soliton mit Geschwindigkeit c = 0.5
c = 0.5;
y1(1,(1:NG)) = KdVSoliton(x,0,c,10,NG);
% 2. Zeitschritt für symmetrische Form der diskreten Ableitung notwendig, 
% Soliton ist mit Geschwindigkeit c um dt vorangeschritten
y1(2,(1:NG)) = KdVSoliton(x,dt,c,10,NG);

% Lösen der DGL mit der Finite-Differenzen-Methode
for j = 3:NT
    for i = 3:NG-2
        yx  = y1(j-1,i+1)-y1(j-1,i-1);
        yx3 = y1(j-1,i+2)-2*y1(j-1,i+1)+2*y1(j-1,i-1)-y1(j-1,i-2);
        y1(j,i) = y1(j-2,i)-(dt/dx)*(6*y1(j-1,i)*yx+yx3/dx^2);    
    end                   
end

%% Lösen der DGL für alle weiteren Geschwindigkeiten
c = 1.;
y2(1,(1:NG)) = KdVSoliton(x,0,c,10,NG);
y2(2,(1:NG)) = KdVSoliton(x,dt,c,10,NG);
for j = 3:NT
    for i = 3:NG-2
        yx  = y2(j-1,i+1)-y2(j-1,i-1);
        yx3 = y2(j-1,i+2)-2*y2(j-1,i+1)+2*y2(j-1,i-1)-y2(j-1,i-2);
        y2(j,i) = y2(j-2,i)-(dt/dx)*(6*y2(j-1,i)*yx+yx3/dx^2);    
    end                   
end

c = 2.;
y3(1,(1:NG)) = KdVSoliton(x,0,c,10,NG);
y3(2,(1:NG)) = KdVSoliton(x,dt,c,10,NG);
for j = 3:NT
    for i = 3:NG-2
        yx  = y3(j-1,i+1)-y3(j-1,i-1);
        yx3 = y3(j-1,i+2)-2*y3(j-1,i+1)+2*y3(j-1,i-1)-y3(j-1,i-2);
        y3(j,i) = y3(j-2,i)-(dt/dx)*(6*y3(j-1,i)*yx+yx3/dx^2);    
    end                   
end

c = 3.;
y4(1,(1:NG)) = KdVSoliton(x,0,c,10,NG);
y4(2,(1:NG)) = KdVSoliton(x,dt,c,10,NG);
for j = 3:NT
    for i = 3:NG-2
        yx  = y4(j-1,i+1)-y4(j-1,i-1);
        yx3 = y4(j-1,i+2)-2*y4(j-1,i+1)+2*y4(j-1,i-1)-y4(j-1,i-2);
        y4(j,i) = y4(j-2,i)-(dt/dx)*(6*y4(j-1,i)*yx+yx3/dx^2);    
    end                   
end

%% Darstellung der Ergebnisse
k = [1 20000 30000 50000 70000 100001];
figure
for i = 1:6
  subplot(3,2,i);
  val1 = plot(x,KdVSoliton(x,dt*(k(i)-1),0.5,10,NG),'linewidth',1, ...
      'color', Colors(1,:));
  hold on;
  plot(x,y1(k(i),:),'linewidth',3, 'color', Colors(1,:),'LineStyle','--');
  hold on;
  val2 = plot(x,KdVSoliton(x,dt*(k(i)-1),1,10,NG),'linewidth',1, ...
      'color', Colors(2,:));
  hold on;
  plot(x,y2(k(i),:),'linewidth',3, 'color', Colors(2,:),'LineStyle','--');
  hold on;
  val3 = plot(x,KdVSoliton(x,dt*(k(i)-1),2,10,NG),'linewidth',1, ...
      'color', Colors(3,:));
  hold on;
  plot(x,y3(k(i),:),'linewidth',3, 'color', Colors(3,:),'LineStyle','--');
  hold on;
  val4 = plot(x,KdVSoliton(x,dt*(k(i)-1),3,10,NG),'linewidth',1, ...
      'color', Colors(4,:));
  hold on;
  plot(x,y4(k(i),:),'linewidth',3, 'color', Colors(4,:),'LineStyle','--');
  grid on;
  axis([min(x) max(x) 0 2]);
  xlabel('{\itx}');
  ylabel('{\it\psi}({\itt})');
  if (i==1)
      legend([val4,val3,val2,val1],{'{\itc} = 3','{\itc} = 2', ...
          '{\itc} = 1','{\itc} = 0.5'})
  end
  titletext = sprintf('{\\itt} = %.1f',t(k(i)));
  texto = text(L/4,1.8,titletext);
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