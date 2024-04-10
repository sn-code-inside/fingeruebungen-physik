% -------------------------------------------------------------------------
% TrebuchetSim02.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet Lagrange-Funktion und 
% Lagrange-Gleichungen des Trebuchets 
% mit Naeherungen 1-4 - Lösung ohne Zwangsbedingung
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

fzero = 1e-5; lW0 = 5;  % Standardparameter
g = 9.81;               % g
h = 10;                 % Drehpunkthöhe des Trebuchets

%% Eingabe Parameter 
% Soll eine Länge LG oder lS gleich Null gesetzt werden, dann einfach 
% fzero eingeben, also z.B. LG = fzero.

% Fall A
B = 4; 
b = 12; 
MG = 10000; 
mW = 100; 
mT = 1040; 
gamma0 = deg2rad(55);
LG = 3; 
lW = 8;  

% Fall B
% B=3; b=15; MG=15000; mW= 50; mT = 2160; gamma0 = deg2rad(40);
% LG = 1.5; lW = 5; 

%% Definiert die verallgemeinerten Koordinaten und weitere Parameter 
% Parameter Trägheitsmomente etc.
JB  = mT*(B^3+b^3)/(3*(B+b));
J0  = mW*b^2 + MG*B^2 + JB;
U0  = g*(MG*B - mW*b - (1/2)*mT*(b-B));
JG  = MG*LG*LG; 
JW  = mW*lW*lW; 
JBG = MG*LG*B;  
JbW = mW*lW*b;  

% Parametersatzstruktur für ODE45 und ODE15s
P.LG =LG; P.lW=lW; P.J0 = J0; P.JG= JG; P.JBG= JBG ; 
P.JbW = JbW; P.JW = JW; P.MG = MG; P.mW = mW; P.mT = mT; 
P.g = g; P.U0 = U0; P.B = B; P.b = b; P.h = h;

% Zeitspanne in s
tmax = 2.5;
% Simulationszeit, Zeit-Vektor
NPoints = 500; 
tstep = tmax/NPoints;             
tv=(0.0:tstep:tmax); 
N=length(tv);         

% Anfangswerte
dgamma0 = 0;
beta0   = deg2rad(-90);  
dbeta0  = 0;
if lW == fzero   %Abfangen falls Schleuderlänge = 0;
    delta0 = 0;
else
    delta0  = asin((-h+b*sin(gamma0))/lW); 
end
ddelta0 = 0;

% Vektor der Anfangsbedingungen
AB=[gamma0;dgamma0;beta0;dbeta0;delta0;ddelta0];         

%% Numerische Integration
% Runge-Kutta-Verfahren-MATLAB ODE45
opts = odeset('Mass',@(t,Y) Mass(t,Y,P),'RelTol',1e-6);
[t,Y] = ode45(@(t,Y) f(t,Y,P),tv,AB,opts);
gamma  = Y(:,1); dgamma = Y(:,2); 
beta   = Y(:,3);  dbeta = Y(:,4);
delta  = Y(:,5); ddelta = Y(:,6); 

% Abbruch Berechnung bei gamma_max = -pi/2 (Balken steht dann senkrecht)
gamma_max = -pi/2; k = 1;
while gamma(k) > gamma_max &&  k < NPoints-1
    k= k+1;
end
if k > 2 kend = k-1; else kend=1; end 
for k=kend:length(t)
    gamma(k) = NaN; dgamma(k)= NaN; 
    beta(k)  = NaN; dbeta(k) = NaN;
    delta(k) = NaN; ddelta(k)= NaN; 
end
%__________________________________________________________________________


%% Koordinaten und Geschwindigkeiten

xGB  = +B*cos(gamma);         % Aufhängung Gegengewicht
zGB  = +B*sin(gamma) ;
xG   = xGB + LG*cos(beta);    % Gegengewicht
zG   = zGB + LG*sin(beta);
dxG  = - B*dgamma.*sin(gamma) - LG*dbeta.*sin(beta);
dzG  = + B*dgamma.*cos(gamma) + LG*dbeta.*cos(beta);

xWB  = -b*cos(gamma);         % Aufhängung Wurfgewicht
zWB  = -b*sin(gamma);
xW   = xWB + lW*cos(delta);   % Wurfgewicht
zW   = zWB + lW*sin(delta);
dxW = +b*sin(gamma).*dgamma - lW*sin(delta).*ddelta;  % Wurfgewicht
dzW = -b*cos(gamma).*dgamma + lW*cos(delta).*ddelta;

    
% Kinetische Energie Wurfgewicht
vW = sqrt(dxW.^2+dzW.^2);
TW = mW*vW.^2/2;

%Wurfweite
alpha = atan2(dzW,dxW);
W  = vW.*cos(alpha).*(vW.*sin(alpha) + ...
     sqrt(vW.^2.*sin(alpha).*sin(alpha)+2*g*(zW>0)))/g; W  = real(W);
[W_max,k_max] = max(W);
alpha_max     = rad2deg(alpha(k_max));  %Abwurfwinkel
t_max         = t(k_max);               %Abwurfzeit 

% Änderung der potentiellen Energie Delta U_pot
dU_pot = -U0*(sin(gamma(k_max)) - sin(gamma(1))) - ...  
         mW*lW*g*(sin(delta(k_max)) - sin(delta(1))) + ... 
         MG*LG*g*(sin(beta(k_max)) - sin(beta(1)));


     
%% Ausgabe Trebuchet- und Effizienzparameter     

% Effizienzparameter
W0    = (2*dU_pot/mW/g)*sqrt(2*zW(k_max)*mW*g/dU_pot + 1);
epsW  = W_max/W0;
epsT  = TW(k_max)/dU_pot;
Ausgabewerte(t_max,alpha_max,vW(k_max),W_max,epsW,epsT,P);

%%  Graphische Ausgabe
% Bild 1 Abwurfzeit, Winkel und Geschwindigkeit
figure();
yyaxis left
p(1)=plot(t,real(W),'Color',Colors(3,:));
p(2)=line([t_max t_max],[0 W_max*1.2],'Color',Colors(3,:),'LineStyle',Style(4));
xlabel('\it{t} \rm in s','FontSize',14);
ylabel('Reichweite \it{W} \rm in m','FontSize',14)
axis([0,tmax,0,W_max*1.2])
yyaxis right
p(3)=plot(t,vW,'Color',Colors(4,:),'LineStyle',Style(1));
hold on
p(4)=plot(t,rad2deg(alpha),'Color',Colors(4,:),'LineStyle',Style(2));
ylabel('\alpha in °, \it{v_w} \rm in m/s ','FontSize',14)
axis([0,tmax,0,max(vW)*2])
h2=legend(p(1:4),'Reichweite','Abwurfzeit(opt.)',...
         'Abwurfgeschwindigkeit','Abwurfwinkel','location','west'); 
set(h2,'FontSize',14)
grid on;
legend box off

% Bild 2 Simulation
figure(); %Simulation
xAmin = -b-lW0; xAmax = +b; yAmin = -b-lW0; yAmax = +b+lW0;
hold on
axis([xAmin,xAmax,yAmin,yAmax]);
axis equal;
PlotTrebuchet(h,1,xW(1),xAmax,Colors,Style);
h2(1)= line([xAmin xW(1)],[-h-lW/2 -h-lW/2]);
h2(2)= line([xW(1) xW(1)],[-h-lW/2 -h]);
set(h2,'Color', Colors(8,:),'LineStyle',Style(2), 'LineWidth',2); 
for k=1:20:k_max
     plot(xW(k),zW(k),'o',...
          'MarkerSize',3,...
          'MarkerEdgeColor',Colors(2,:),...
          'MarkerFaceColor',Colors(2,:));             %mw
     plot(xG(k),zG(k),'o',...
          'MarkerSize',6,...
          'MarkerEdgeColor',Colors(4,:),...
          'MarkerFaceColor',Colors(4,:));              %MG
     line([xGB(k) xWB(k)],[zGB(k) zWB(k)],'Color', Colors(3,:),...
         'LineStyle',Style(1));
     line([xGB(k) xG(k)],[zGB(k) zG(k)],'Color', Colors(4,:),...
         'LineStyle',Style(1));
     line([xWB(k) xW(k)],[zWB(k) zW(k)],'Color', Colors(2,:),...
          'LineStyle',Style(1));
     pause(.01)
end
plot(xWB,zWB,'Color',Colors(2,:));
hold on
plot(xGB,zGB,'Color',Colors(4,:));
line([xGB(k_max) xG(k_max)],[zGB(k_max) zG(k_max)],'Color', Colors(4,:),...
      'LineStyle',Style(1), 'LineWidth',2);
line([xWB(k_max) xW(k_max)],[zWB(k_max) zW(k_max)],'Color', Colors(2,:),...
      'LineStyle',Style(1), 'LineWidth',2);
ha(1)=plot(xW(k_max),zW(k_max),'o','MarkerSize',4,'MarkerEdgeColor',Colors(2,:),...
    'MarkerFaceColor',Colors(2,:)); %mW
ha(2)=plot(xG(k_max),zG(k_max),'o','MarkerSize',7,'MarkerEdgeColor',Colors(4,:),...
    'MarkerFaceColor',Colors(4,:)); %MG
ha(3)= line([xGB(k_max) xWB(k_max)],[zGB(k_max) zWB(k_max)],'Color', Colors(3,:),...
     'LineStyle',Style(1), 'LineWidth',2);
grid on;
xlabel('\it{x} \rm in m','FontSize',14);
ylabel('\it{z} \rm in m','FontSize',14)
hb=legend(ha,{'Wurfgewicht','Gegengewicht', 'Balken'},'location','northeast'); 
legend box off
set(hb,'FontSize',14)
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------




function M = Mass(t,q,P)
    % Extract parameters
    JG  = P.JG;
    JW  = P.JW;
    JbW = P.JbW;
    JBG = P.JBG;
    J0  = P.J0;

    % Mass matrix function
    M = zeros(6,6);
    M(1,1) = 1;
    M(2,2) = J0;
    M(3,3) = 1;
    M(4,4) = JG;
    M(5,5) = 1;
    M(6,6) = JW;
    M(2,4) = +JBG*cos(q(1)-q(3));
    M(4,2) = M(2,4);
    M(2,6) = -JbW*cos(q(1)-q(5));
    M(6,2) = M(2,6);
end

function dYdt = f(t,Y,P)
    % Extract parameters
    LG  = P.LG;  
    lW  = P.lW;
    MG  = P.MG;  
    mW  = P.mW;
    JBG = P.JBG; 
    JbW = P.JbW; 
    g   = P.g;   
    U0  = P.U0;

    % DGL
    dYdt = [Y(2)
        -JBG*Y(4)^2*sin(Y(1)-Y(3))+JbW*Y(6)^2*sin(Y(1)-Y(5))-U0*cos(Y(1));
        Y(4)
        +JBG*Y(2)^2*sin(Y(1)-Y(3)) - MG*LG*g*cos(Y(3))
        Y(6)
        -JbW*Y(2)^2*sin(Y(1)-Y(5)) - mW*lW*g*cos(Y(5))];
end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------


