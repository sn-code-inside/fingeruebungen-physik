% -------------------------------------------------------------------------
% SinusGordonSolitonen.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik des Kontinuums" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm propagiert Soliton-Lösungen der Sinus-Gordon-Gleichung, indem
% sowhohl analytische Lösungen dargestellt werden als auch die
% Differentialgleichung direkt mit der Finite-Differenzen-Methode gelöst
% wird.
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

%% Numerische Lösung für Soliton vom Typ 1 (kink)
% Bereite Gitter vor
y1 = zeros(NT,NG);

% Anfangsbedingung: Soliton mit Geschwindigkeit v bei
% Maximalgeschwindigkeit c
c = 10.;
v = 2.;
y1(1,(1:NG)) = SGSoliton1(x,0,v,c,NG);
% 2. Zeitschritt für symmetrische Form der diskreten Ableitung notwendig, 
% Soliton ist mit Geschwindigkeit c um dt vorangeschritten
y1(2,(1:NG)) = SGSoliton1(x,dt,v,c,NG);
% Setzen der Ränder des Gitters auf die erwarteten Werte
Soliton10       = SGSoliton1(x,0,v,c,NG);
y1((1:NT),1)    = Soliton10(1);
y1((1:NT),2)    = Soliton10(2);
y1((1:NT),NG-1) = Soliton10(NG-1);
y1((1:NT),NG)   = Soliton10(NG);

% Lösen der DGL mit der Finite-Differenzen-Methode
for j = 3:NT
    for i = 3:NG-2
        yx2 = y1(j-1,i+1)+y1(j-1,i-1)-2*y1(j-1,i);
        y1(j,i) = 2*y1(j-1,i)-y1(j-2,i)+c^2*dt^2*(yx2/dx^2 ...
            -sin(y1(j-1,i)));    
    end                   
end

%% Numerische Lösung für Soliton vom Typ 2 (anti-kink)
% Bereite Gitter vor
y2 = zeros(NT,NG);

% Anfangsbedingung: Soliton mit Geschwindigkeit v bei
% Maximalgeschwindigkeit c
c = 10.;
v = 2.;
y2(1,(1:NG)) = SGSoliton2(x,0,v,c,NG);
% 2. Zeitschritt für symmetrische Form der diskreten Ableitung notwendig, 
% Soliton ist mit Geschwindigkeit c um dt vorangeschritten
y2(2,(1:NG)) = SGSoliton2(x,dt,v,c,NG);
% Setzen des rechten Randes des Gitters auf die Randwerte des Solitons
Soliton20       = SGSoliton2(x,0,v,c,NG);
y2((1:NT),1)    = Soliton20(1);
y2((1:NT),2)    = Soliton20(2);
y2((1:NT),NG-1) = Soliton20(NG-1);
y2((1:NT),NG)   = Soliton20(NG);

% Lösen der DGL mit der Finite-Differenzen-Methode
for j = 3:NT
    for i = 3:NG-2
        yx2 = y2(j-1,i+1)+y2(j-1,i-1)-2*y2(j-1,i);
        y2(j,i) = 2*y2(j-1,i)-y2(j-2,i)+c^2*dt^2*(yx2/dx^2 ...
            -sin(y2(j-1,i)));    
    end                   
end

%% Darstellung der Ergebnisse
% Abbildung mit der Variation verschiedenenr Paramter (Typ 1, kink)
subplot(2,1,1);
plot(x,SGSoliton1(x,0,0.01,1,NG),'linewidth',3, 'color', ...
    Colors(1,:));
hold on;
plot(x,SGSoliton1(x,0,0.2,1,NG),'linewidth',3, 'color', ...
    Colors(3,:));
hold on;
plot(x,SGSoliton1(x,0,0.8,1,NG),'linewidth',3, 'color', ...
    Colors(4,:));
hold on;
plot(x,SGSoliton1(x,0,0.999,1,NG),'linewidth',3, 'color', ...
    Colors(5,:));
grid on;
axis([-10 10 0 2*pi]);
xlabel('{\itx}');
ylabel('{\it\psi}({\itt})');
legend('{\itv/c}=0.01','{\itv/c}=0.2','{\itv/c}=0.8','{\itv/c}=0.99');
set(gca,'FontSize',16,'FontName','Times');
% Abbildung mit der Variation verschiedenenr Paramter (Typ 2, anti-kink)
subplot(2,1,2);
plot(x,SGSoliton2(x,0,0.01,1,NG),'linewidth',3, 'color', ...
    Colors(1,:));
hold on;
plot(x,SGSoliton2(x,0,0.2,1,NG),'linewidth',3, 'color', ...
    Colors(3,:));
hold on;
plot(x,SGSoliton2(x,0,0.8,1,NG),'linewidth',3, 'color', ...
    Colors(4,:));
hold on;
plot(x,SGSoliton2(x,0,0.999,1,NG),'linewidth',3, 'color', ...
    Colors(5,:));
grid on;
axis([-10 10 0 2*pi]);
xlabel('{\itx}');
ylabel('{\it\psi}({\itt})');
legend('{\itv/c}=0.01','{\itv/c}=0.2','{\itv/c}=0.8','{\itv/c}=0.99');
set(gca,'FontSize',16,'FontName','Times');


% Abbildung zur numerischen Lösung (bewegte Solitonen)
k = [1 20000 30000 50000 70000 100001];
figure
for i = 1:6
  subplot(3,2,i);
  plot(x,SGSoliton1(x,dt*(k(i)-1),v,c,NG),'linewidth',2, 'color', ...
      Colors(3,:));
  hold on;
  plot(x,SGSoliton2(x,dt*(k(i)-1),v,c,NG),'linewidth',2, 'color', ...
      Colors(3,:));
  hold on;
  plot(x,y1(k(i),:),'linewidth',5, 'color', Colors(2,:),'LineStyle','--');
  hold on;
  plot(x,y2(k(i),:),'linewidth',5, 'color', Colors(4,:),'LineStyle','--');
  grid on;
  axis([min(x) max(x) 0 2*pi]);
  xlabel('{\itx}');
  ylabel('{\it\psi}({\itt})');
  titletext = sprintf('{\\itt} = %.1f',t(k(i)));
  texto = text(-18,3,titletext);
  set(texto,'FontSize',16,'FontWeight','normal','FontName','Times');
  set(gca,'FontSize',16,'FontName','Times');
end
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


%% Funktionen
% analytische Lösung eines Solitions des ersten Typs (kink) in der
% Sinus-Gordon-Gleichung
function  y  =  SGSoliton1(x,t,v,c,NG)
    y(1:NG) = 4*atan(exp((x-v*t)/sqrt(1-(v/c)^2)));
end

% analytische Lösung eines Solitions des zweiten Typs (anti-kink) in der
% Sinus-Gordon-Gleichung
function  y  =  SGSoliton2(x,t,v,c,NG)
    y(1:NG) = 4*atan(exp(-(x-v*t)/sqrt(1-(v/c)^2)));
end
