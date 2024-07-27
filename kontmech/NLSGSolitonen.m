% -------------------------------------------------------------------------
% NLSGSolitonen.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik des Kontinuums" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm propagiert Soliton-Lösungen der nichtlinearen
% Schrödinger-Gleichung, indem sowhohl analytische Lösungen dargestellt
% werden als auch die Differentialgleichung direkt mit der
% Finite-Differenzen-Methode gelöst wird.
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

%% Parameter
L = 80;                               % Länge des Gitters

NG = 801;                             % Gitterpunkte im Ort
dx = L/(NG-1);                        % Abstand zweier Gitterpunkte
x(:,1) = -L/4+(0:NG-1)*dx;
                
Tmax = 20.;                           % maximale Zeit
NT = 100001;                          % Schritte in der Zeit
dt = Tmax/(NT-1);                     % Zeitschritt
t(:)= (0:NT-1)*dt;                    % Zeitgitter

%% Soliton des Fundamentaltyps
% Bereite Gitter vor
y1 = zeros(NT,NG);

% Anfangsbedingung: Soliton mit Geschwindigkeit v
vr = 1.5;
gr = 0.8;

% Setzen der Parameter für die Rechnung
v  = vr;
g  = gr;
x0 = -1./v;
t0 = 1./v^2;
A0 = v/sqrt(g);
y1(1,(1:NG)) = NLSGSoliton1(x,0,A0,v,t0,x0,NG);
% 2. Zeitschritt für symmetrische Form der diskreten Ableitung notwendig, 
% Soliton ist mit Geschwindigkeit c um dt vorangeschritten
y1(2,(1:NG)) = NLSGSoliton1(x,dt,A0,v,t0,x0,NG);

% Lösen der DGL mit der Finite-Differenzen-Methode
for j = 3:NT
    for i = 3:NG-2
        yx2 = y1(j-1,i+1)+y1(j-1,i-1)-2*y1(j-1,i);
        y1(j,i) = y1(j-2,i) + 1i*dt/dx^2*yx2 ...
            + 2i*dt*g*abs(y1(j-1,i))^2*y1(j-1,i);
    end                   
end

%% Darstellung der Ergebnisse 
% Abbildung mit der Variation verschiedenenr Paramter, zunächst v
figure
subplot(2,1,1)
v=0.5; g=1.;
plot(x,abs(NLSGSoliton1(x,0,v/sqrt(g),v,1./v^2,-1./v,NG)).^2, ...
    'linewidth', 3, 'color', Colors(1,:));
hold on;
v=1.;
plot(x,abs(NLSGSoliton1(x,0,v/sqrt(g),v,1./v^2,-1./v,NG)).^2, ...
    'linewidth', 3, 'color', Colors(2,:));
hold on;
v=2.;
plot(x,abs(NLSGSoliton1(x,0,v/sqrt(g),v,1./v^2,-1./v,NG)).^2, ...
    'linewidth', 3, 'color', Colors(3,:));
hold on;
v=3.;
plot(x,abs(NLSGSoliton1(x,0,v/sqrt(g),v,1./v^2,-1./v,NG)).^2, ...
    'linewidth', 3, 'color', Colors(4,:));
axis([-6 6 0 10]);
xlabel('{\itx}');
ylabel('|{\it\psi}|^2');
legend('{\itv}=0.5','{\itv}=1.','{\itv}=2.','{\itv}=3.');
set(gca,'FontSize',16,'FontName','Times');
% Variation von g
subplot(2,1,2)
v=1; g=0.5;
plot(x,abs(NLSGSoliton1(x,0,v/sqrt(g),v,1./v^2,-1./v,NG)).^2, ...
    'linewidth', 3, 'color', Colors(1,:));
hold on;
g=1.;
plot(x,abs(NLSGSoliton1(x,0,v/sqrt(g),v,1./v^2,-1./v,NG)).^2, ...
    'linewidth', 3, 'color', Colors(2,:));
hold on;
g=2.;
plot(x,abs(NLSGSoliton1(x,0,v/sqrt(g),v,1./v^2,-1./v,NG)).^2, ...
    'linewidth', 3, 'color', Colors(3,:));
hold on;
g=3.;
plot(x,abs(NLSGSoliton1(x,0,v/sqrt(g),v,1./v^2,-1./v,NG)).^2, ...
    'linewidth', 3, 'color', Colors(4,:));
axis([-6 6 0 2.5]);
xlabel('{\itx}');
ylabel('|{\it\psi}|^2');
legend('{\itg}=0.5','{\itg}=1.','{\itg}=2.','{\itg}=3.');
set(gca,'FontSize',16,'FontName','Times');

% Abbildung zur numerischen Lösung (bewegte Solitonen)
k = [1 20000 30000 50000 70000 100001];
figure
for i = 1:6
  subplot(3,2,i);
  plot(x,abs(NLSGSoliton1(x,dt*(k(i)-1),A0,vr,t0,x0,NG)).^2, ...
      'linewidth', 2, 'color', Colors(3,:));
  hold on;
  plot(x,abs(y1(k(i),:)).^2,'linewidth',5, 'color', Colors(4,:), ...
      'LineStyle','--');
  axis([min(x) max(x) 0 3]);
  xlabel('{\itx}');
  ylabel('{|\it\psi}({\itt})|^2');
  titletext = sprintf('{\\itt} = %.1f',t(k(i)));
  texto = text(-18,2.5,titletext);
  set(texto,'FontSize',16,'FontWeight','normal','FontName','Times');
  set(gca,'FontSize',16,'FontName','Times');
end
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


%% Funktionen
% analytische Lösung eines Solitions des ersten Typs (Fundamentallösung)
% in der nichtlinearen Schrödingergleichung
function  y  =  NLSGSoliton1(x,t,A0,v,t0,x0,NG)
    y(1:NG) = A0*sech((t-x/v)/t0).*exp(-1i*x/x0);
end