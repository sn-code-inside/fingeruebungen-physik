% -------------------------------------------------------------------------
% TrebuchetSim01.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Floating Arm Trebuchet
% 
% Programm berechnet Lagrange-Funktion und 
% Lagrange-Gleichungen des Floating Arm Trebuchets 
% 
% Benutzt die symbolische Euler/Lagrange Berechnung nach 
% Morten Veng (2021). Euler-Lagrange Solver 
% (https://www.mathworks.com/matlabcentral/fileexchange/93275-euler-lagrange-solver), 
% MATLAB Central File Exchange. Retrieved December 27, 2021.
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

fzero = 1e-5; lS0 = 5;  % Standardparameter
g=9.81;                 % g
h = 10;                 % Drehpunkthöhe des Trebuchets

%% Eingabe Parameter 
% -------------------------------------------------------------------------
%

% Fall A
B=4; b=12; MG=10000; mW=100; mT = 2040; 
lS = 8; 

% Fall B
% B=3; b=15; MG=15000; mW= 50; mT = 2160; gamma0 = deg2rad(40);
% lS = 5; 

% -------------------------------------------------------------------------

% Parametersatzstruktur für ODE45 und ODE15s
P.lS=lS;  
JB  = mT*(B^3+b^3)/(3*(B+b));
P.LG =0; P.J0 = 0; P.JG= 0; P.JBG= 0 ; P.JB = JB; P.lW = 0;
P.JbW = 0; P.JW = 0; P.MG = MG; P.mW = mW; P.mT = mT; 
P.g = g; P.U0 = 0; P.B = B; P.b = b; P.h = h;

% Zeitspanne in s
tmax    = 4; NPoints =200; tstep = tmax/NPoints;             
tv=(0.0:tstep:tmax); N=length(tv);          %Simulationszeit, Zeit-Vektor

% Anfangswerte
gamma0  = asin(h/b);  
gamma0d = 0;                                %AB Winkel in rad
delta0  = deg2rad(0);  
delta0d = 0;                                %AB Winkel in rad
AB=[gamma0;gamma0d;delta0;delta0d];         %AB Vektor

xw0 = -(b+B)*cos(gamma0) + ls*cos(delta0);
zw0 = - b*sin(gamma0)    + ls*sin(delta0);
xG0 = 0;
zG0 = B*sin(gamma0);

% Runge-Kutta-Verfahren-MATLAB ODE45
opt=odeset('AbsTol',1.e-9,'RelTol',1.e-5);           
[t,Y]=ode45(@(t,Y)DGL(t,Y,P),tv,AB,opt); 
gamma   = Y(:,1); dgamma  = Y(:,2); delta   = Y(:,3); ddelta  = Y(:,4);
    
% Abbruch Berechnung bei zG = - h
k = 1;
while  gamma(k) > -pi/2  && k<NPoints
    k=k+1;
end
if k > 2 kend = k-1; else kend=1; end 
for k=kend:length(t)
    gamma(k) = NaN; dgamma(k) = NaN; delta(k) = NaN; ddelta(k) = NaN;
end
    

% Koordinaten und Geschwindigkeiten aus Winkeln
xR   = -B*cos(gamma);       % Rolle
zR   = 0*t;   
xWB  = -(b+B)*cos(gamma);   % Aufhängung Wurfgewicht
zWB  =  -b*sin(gamma);
xG   = 0*t;    % Gegengewicht
zG   = B*sin(gamma);
xW   = xWB + lS*cos(delta);   % Wurfgewicht
zW   = zWB + lS*sin(delta);
dxW = (b+B)*sin(gamma).*dgamma -lS*sin(delta).*ddelta;  % Wurfgewicht
dzW = -b*cos(gamma).*dgamma  + lS*cos(delta).*ddelta;

% Kinetische Energie Wurfgewicht
vW = sqrt(dxW.^2+dzW.^2);
TW = mW*vW.^2/2;

% Wurfweite
alpha = atan2(dzW,dxW);
W  = vW.*cos(alpha).*(vW.*sin(alpha) + ...
     sqrt(vW.^2.*sin(alpha).*sin(alpha)+2*g*(zW>0)))/g; W  = real(W);
[W_max,k_max] = max(W);
alpha_max     = rad2deg(alpha(k_max));  %Abwurfwinkel
t_max         = t(k_max);               %Abwurfzeit 

% Änderung der potentiellen Energie Delta U_pot
dU_pot = -MG*B*g*(sin(gamma(k_max)) - sin(gamma(1)));


     
%% Ausgabe Trebuchet-Parameter und Effizienzparameter     

% Effizienzparameter
W0    = (2*dU_pot/mW/g)*sqrt(2*zW(k_max)*mW*g/dU_pot + 1);
epsW  = W_max/W0;
epsT  = TW(k_max)/dU_pot;
Ausgabewerte(t_max,alpha_max,vW(k_max),W_max,epsW,epsT,P);

%%  Graphische Ausgabe

% Bild 1 Reichweite, Winkel, Abwurfgeschwindigkeit
figure();
yyaxis left
plot(t,real(W),'Color',Colors(3,:));
line([t_max t_max],[0 W_max],'Color',Colors(3,:),'LineStyle',Style(4));
xlabel('\it{t} \rm in s','FontSize',14);
ylabel('Reichweite \it{W} \rm in m','FontSize',14)
yyaxis right
plot(t,vW,'Color',Colors(4,:),'LineStyle',Style(1));
hold on
plot(t,rad2deg(alpha),'Color',Colors(4,:),'LineStyle',Style(2));
ylabel('\alpha, {v_w}','FontSize',14)
axis([0,tmax,0,180])
grid on;
h1=legend('Reichweite','Abwurfzeit(opt.)',...
         'Abwurfgeschwindigkeit','Abwurfwinkel','location','northwest'); 
set(h1,'FontSize',14)
legend box off

% Bild 3 bzw. Simulation
xAmin = -b-B-lS;
xAmax = +b+lS;
yAmin = -b-lS;
yAmax = +b+lS;

figure(); % Simulation
line([xWB(1) xG(1)],[zWB(1) zG(1)],'Color', Colors(2,:),'LineStyle',Style(3));
line([xWB(1) xW(1)],[zWB(1) zW(1)],'Color', Colors(3,:),'LineStyle',Style(3));
axis([xAmin,xAmax,yAmin,yAmax]);
plotFloatArmTrebShape(h,0.5,xAmin,xAmax,Colors,Style);
pbaspect([1 1 1]);              % Aspektverhältnis Achsen
hold on
for k=1:10:k_max
     ha(1)=plot(xW(k),zW(k),'o',...
          'MarkerSize',4,...
          'MarkerEdgeColor',Colors(2,:),...
          'MarkerFaceColor',Colors(2,:)'); %mw
     ha(2)=plot(xG(k),zG(k),'s',...
          'MarkerSize',7,...
          'MarkerEdgeColor',Colors(4,:),...
          'MarkerFaceColor',Colors(4,:));  %MG
     ha(3)=plot(xR(k),zR(k),'o',...
          'MarkerSize',6,...
          'MarkerEdgeColor',Colors(8,:),...
          'MarkerFaceColor',Colors(8,:));  %Rolle
     line([xWB(k) xG(k)],[zWB(k) zG(k)],'Color', Colors(3,:),...
          'LineStyle',Style(1));
     line([xWB(k) xW(k)],[zWB(k) zW(k)],'Color', Colors(2,:),...
          'LineStyle',Style(1));
     pause(.01)
end
plot(xW,zW,'Color',Colors(2,:));
hold on
plot(xG,zG,'Color',Colors(4,:));
line([xWB(k_max) xG(k_max)],[zWB(k_max) zG(k_max)],'Color', Colors(3,:),...
     'LineStyle',Style(1), 'LineWidth',2);
line([xWB(k_max) xW(k_max)],[zWB(k_max) zW(k_max)],'Color', Colors(2,:),...
     'LineStyle',Style(1), 'LineWidth',2);
plot(xR(k_max),zR(k_max),'o','MarkerSize',6,'MarkerEdgeColor',Colors(8,:),...
          'MarkerFaceColor',Colors(8,:));  %Rolle
plot(xW(k_max),zW(k_max),'o','MarkerSize',4,'MarkerEdgeColor',Colors(2,:),...
    'MarkerFaceColor',Colors(2,:)'); %mw
plot(xG(k_max),zG(k_max),'s','MarkerSize',7,'MarkerEdgeColor',Colors(2,:),...
    'MarkerFaceColor',Colors(4,:)); %MG
line([xAmin xAmax],[-h -h],'Color', Colors(8,:),...
     'LineStyle',Style(4), 'LineWidth',1);
grid on;
axis([xAmin,xAmax,yAmin,yAmax]);
xlabel('\it{x} \rm','FontSize',14); ylabel('\it{z} \rm','FontSize',14)
hb=legend(ha,{'Wurfgewicht','Gegengewicht','Rolle'},'location','southeast'); 
legend box off
set(hb,'FontSize',14)
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


%% DGL System
function dY = DGL(t,Y,P)
    % Y(1)=gamma
    % Y(2)=gamma_dot
    % Y(3)=delta
    % Y(4)=delta_dot
    
    B = P.B;
    b = P.b;
    MG= P.MG;
    mw= P.mW;
    ls= P.lS;
    g = P.g;
    JB= P.JB;
    mT= P.mT;

    dY    = zeros(4,1);
    gamma = Y(1);
    dY(1) = Y(2);
    dgamma= Y(2);
    delta = Y(3);
    dY(3) = Y(4);
    ddelta= Y(4);
    ddgamma = dY(2);
    dddelta = dY(4);
    dY(2) = (B^2*MG*dgamma^2*sin(2*gamma) - B*g*mT*cos(gamma) - ...
        2*B*MG*g*cos(gamma) + b*g*mT*cos(gamma) + 2*b*g*mw*cos(gamma) -...
        B^2*dgamma^2*mw*sin(2*gamma) - 2*B*b*dgamma^2*mw*sin(2*gamma) +...
        2*b*ddelta^2*ls*mw*cos(delta)*sin(gamma) - ...
        2*b*ddelta^2*ls*mw*cos(gamma)*sin(delta) + ...
        2*B*dddelta*ls*mw*sin(delta)*sin(gamma) + ...
        2*b*dddelta*ls*mw*cos(delta)*cos(gamma) + ...
        2*b*dddelta*ls*mw*sin(delta)*sin(gamma) + ...
        2*B*ddelta^2*ls*mw*cos(delta)*sin(gamma))/...
        (2*(JB + B^2*mw + b^2*mw + 2*B*b*mw + B^2*MG*cos(gamma)^2 - ...
        B^2*mw*cos(gamma)^2 - 2*B*b*mw*cos(gamma)^2));
    dY(4) = (B*dgamma^2*cos(gamma)*sin(delta) - g*cos(delta) - ...
        b*dgamma^2*cos(delta)*sin(gamma) + b*dgamma^2*cos(gamma)*sin(delta)...
        + B*ddgamma*sin(delta)*sin(gamma) + b*ddgamma*cos(delta)*cos(gamma)...
        + b*ddgamma*sin(delta)*sin(gamma))/ls;
end

%% Funktion fuer die graphische Ausgabe

function plotFloatArmTrebShape(h,base,groundl,groundr,Colors,Style)
    line([groundl groundr],[-h -h],'Color', Colors(8,:),...
         'LineStyle',Style(4), 'LineWidth',1);
    line([-base -base],[-h h/2],'Color', Colors(8,:),...
         'LineStyle',Style(2), 'LineWidth',2); 
    line([-5*base -base],[-h 0],'Color', Colors(8,:),...
         'LineStyle',Style(2), 'LineWidth',2); 
    line([5*base base],[-h 0],'Color', Colors(8,:),...
         'LineStyle',Style(2), 'LineWidth',2); 
    line([+base base],[-h h/2],'Color', Colors(8,:),...
         'LineStyle',Style(2), 'LineWidth',2); 
    line([+base -base],[-h -h],'Color', Colors(8,:),...
         'LineStyle',Style(2), 'LineWidth',2); 
    line([+4*base -10*base],[0 0],'Color', Colors(8,:),...
         'LineStyle',Style(2), 'LineWidth',2); 
    line([+base -base],[h/2 h/2],'Color', Colors(8,:),...
         'LineStyle',Style(2), 'LineWidth',2); 
end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------
