% -------------------------------------------------------------------------
% KdV04.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik des Kontinuums" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm propagiert zwei Solitonien der Korteweg-De-Fries-Gleichung,
% die sich anfänglich nicht überschneiden, dann interagieren und
% wieder trennen. Die Phasenverschiebung im Vergleich zu einer
% ungestörten Propagation einzelner Solitionen wird extrahiert.
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
                
Tmax = 30.;                           % maximale Zeit
NT = 100001;                          % Schritte in der Zeit
dt = Tmax/(NT-1);                     % Zeitschritt
t(:)= (0:NT-1)*dt;                    % Zeitgitter

%% Bereite Gitter vor
ys = zeros(NT,NG);

%% Integration der Lösung
% Zwei anfänglich getrennte Solitonen bewegen sich mit unterschiedlichen
% Geschwindigkeiten nach rechts. Das linke Soliton ist schneller,
% sodass sich beide nach einiger Zeit Überschneiden.
c1  = 2.;
x01 = 10;
c2  = 1.;
x02 = 20;
ys(1,(1:NG)) = KdVSoliton(x,0,c1,x01,NG) + KdVSoliton(x,0,c2,x02,NG);
% 2. Zeitschritt für symmetrische Form der diskreten Ableitung notwendig, 
% Soliton ist mit Geschwindigkeit c um dt vorangeschritten
ys(2,(1:NG)) = KdVSoliton(x,dt,c1,x01,NG) + KdVSoliton(x,dt,c2,x02,NG);

% Lösen der DGL mit der Finite-Differenzen-Methode
for j = 3:NT
    for i = 3:NG-2
        yx  = ys(j-1,i+1)-ys(j-1,i-1);
        yx3 = ys(j-1,i+2)-2*ys(j-1,i+1)+2*ys(j-1,i-1)-ys(j-1,i-2);
        ys(j,i) = ys(j-2,i)-(dt/dx)*(6*ys(j-1,i)*yx+yx3/dx^2);    
    end                   
end

% Extrahieren der beiden Maxima für t(66666)
nmax = 1;
for i = 2:NG-1
    if ((ys(66666,i) > ys(66666,i-1)) && (ys(66666,i) > ys(66666,i+1)) ...
            && (ys(70000,i) > 0.1))
        x_maxima1(nmax) = x(i);
        nmax = nmax+1;
    end
end                   

% Extrahieren der beiden Maxima für t(NT)
nmax = 1;
for i = 2:NG-1
    if ((ys(NT,i) > ys(NT,i-1)) && (ys(NT,i) > ys(NT,i+1)) ...
            && (ys(NT,i) > 0.1))
        x_maxima2(nmax) = x(i);
        nmax = nmax+1;
    end
end    

%% Darstellung der Ergebnisse
k = [1 20000 30000 50000 66666 100001];
figure
for i = 1:6
  subplot(3,2,i);
  plot(x,KdVSoliton(x,dt*(k(i)-1),c1,x01,NG),'linewidth',2, 'color', ...
      Colors(1,:));
  hold on;
  plot(x,KdVSoliton(x,dt*(k(i)-1),c2,x02,NG),'linewidth',2, 'color', ...
      Colors(2,:));
  hold on;
  plot(x,ys(k(i),:),'linewidth',3, 'color', Colors(3,:));
  if (k(i) == 66666)
      xline(x_maxima1(1),'linewidth',2,'color',Colors(2,:));
      xline(x02+c2*t(k(i)),'linewidth',2,'color',Colors(2,:));
      labeltext = sprintf('{\\it\\Deltax} = %.1f', ...
          x02+c2*t(k(i))-x_maxima1(1));
      text1 = text(x_maxima1(1)-15,0.8,labeltext);
      set(text1,'FontSize',16,'FontWeight','normal','FontName','Times');
      xline(x_maxima1(2),'linewidth',2,'color',Colors(1,:));
      xline(x01+c1*t(k(i)),'linewidth',2,'color',Colors(1,:));
      labeltext = sprintf('{\\it\\Deltax} = %.1f', ...
          x01+c1*t(k(i))-x_maxima1(2));
      text2 = text(x_maxima1(2)+2,1.1,labeltext);
      set(text2,'FontSize',16,'FontWeight','normal','FontName','Times');
  end
  if (k(i) == NT)
      xline(x_maxima2(1),'linewidth',2,'color',Colors(2,:));
      xline(x02+c2*t(k(i)),'linewidth',2,'color',Colors(2,:));
      labeltext = sprintf('{\\it\\Deltax} = %.1f', ...
          x02+c2*t(k(i))-x_maxima2(1));
      text3 = text(x_maxima2(1)-15,0.8,labeltext);
      set(text3,'FontSize',16,'FontWeight','normal','FontName','Times');
      xline(x_maxima2(2),'linewidth',2,'color',Colors(1,:));
      xline(x01+c1*t(k(i)),'linewidth',2,'color',Colors(1,:));
      labeltext = sprintf('{\\it\\Deltax} = %.1f', ...
          x01+c1*t(k(i))-x_maxima2(2));
      text4 = text(x01+c1*t(k(i))-16,1.1,labeltext);
      set(text4,'FontSize',16,'FontWeight','normal','FontName','Times');
  end
  grid on;
  axis([min(x) max(x) 0 1.2]);
  xlabel('{\itx}');
  ylabel('{\it\psi}({\itt})');
  titletext = sprintf('{\\itt} = %.1f',t(k(i)));
  texto = text(L/5,1.1,titletext);
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