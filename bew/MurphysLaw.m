% -------------------------------------------------------------------------
% MurphysLaw.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet Dynamik eines fallenden "Butterbrotes"
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];


%% Numerische Lösung

% Parameter
LB      =  20;              % Länge Brett in cm
H       =  75;              % Höhe Tisch in  cm
H       =  215 ;            % Höhe Tisch in  cm
g       =  981.0;           % g in cm/s^2
js      =  LB^2/12;         % Trägheitsmoment /m in cm^2
phi0    =  0;               % Anfangswinkel
dphi0   =  0.0;             % Anfangswinkelgeschwindigkeit in 1/s
rs0     =  LB/10;           % Anfangsposition
drs0    =  0.0;             % Anfangsgeschwindigkeit in m/s

tmax    =  0.5;
tspan   =  linspace(0.0,tmax,30); % t-Bereich
tspan2  =  linspace(0.0,2*tmax,60); % t-Bereich

% Parametersatz für ode45
P1.g    =  g;
P1.js   =  js;
P1.H    =  H;
% Achsenparameter
if H == 75
    yAu = -90;
    yAo = 10;
    xAu = -30;
    xAo = 30; 
else
    yAu = -220;
    yAo = 10;
    xAu = -75;
    xAo = 75; 
end
fprintf('\n ');
fprintf('\n LB   = %8.1f cm',  LB);
fprintf('\n js   = %8.1f cm^2', js);
fprintf('\n H    = %8.1f cm', H);
fprintf('\n g    = %8.1f cm/s^2', g);
fprintf('\n ');

%% Berechnungen

% Anfangswerte für ode45
AB=[phi0,dphi0,rs0,drs0]; % AB für ode45

%%
% Dynamikberechnung 

% Achtung: Winkel im Uhrzeigersinn werden positiv gerechnet!

% Abrollphase
options1 = odeset('AbsTol',1.e-9,'RelTol',1.e-7,'events',@MyEvent1);
[t, Y, TE, YE, IE] = ode45(@dgl_Brett1, tspan, AB, options1, P1); 
% interp1 equal intervals for ode output
Yeq = interp1(t,Y,tspan);
phi  = Yeq(:,1);
dphi = Yeq(:,2);
rs   = Yeq(:,3);
drs  = Yeq(:,4);
ZB   = g*cos(phi)-2*drs.*dphi;
Lend = Drehung(phi,rs,LB);
xo = Lend(1,:);
xu = Lend(2,:);
yo = Lend(3,:);
yu = Lend(4,:);
ltest = sqrt((xo-xu).^2+(yo-yu).^2);

% Freie Rotation
phiF  = YE(1);
dphiF = YE(2);
rsF   = YE(3);
drsF  = YE(4);

xsF   = +rsF*cos(phiF);
ysF   = -rsF*sin(phiF);

dxsF  = +drsF*cos(phiF);
dysF  = -drsF*sin(phiF);

xs2   = xsF + dxsF*tspan2;
ys2   = ysF + dysF*tspan2 - g*tspan2.^2/2;
phi2  = phiF + dphiF*tspan2; 
rs2   = sqrt(xs2.^2+ys2.^2);

yo2   = ys2 + LB*sin(phi2)/2;
yu2   = ys2 - LB*sin(phi2)/2;
xo2   = xs2 - LB*cos(phi2)/2;
xu2   = xs2 + LB*cos(phi2)/2;


%% 
% Graphische Ausgabe

% Abrollphase
fig = figure();
set(fig,'defaultAxesColorOrder',[Colors(2,:); Colors(3,:)]);
yyaxis left;
hold on; 
rp(1) = plot(tspan,rad2deg(phi));
set(rp(1),'Color',Colors(2,:), 'LineWidth',2,'LineStyle',Style(1));
rp(2) = plot(TE+tspan2,rad2deg(phi2));
set(rp(2),'Color',Colors(2,:), 'LineWidth',2,'LineStyle',Style(2));
rp(3) = plot(tspan,rs);
set(rp(3),'Color',Colors(4,:), 'LineWidth',2,'LineStyle',Style(1));
rp(4) = plot(TE+tspan2,rs2);
set(rp(4),'Color',Colors(4,:), 'LineWidth',2,'LineStyle',Style(2));
grid on
xlabel('Zeit t in s','FontSize',14)
ylabel('phi in °, rs in cm','FontSize',14)
set(gca,'FontSize',16);

yyaxis right;
hold on; 
rp(5) = plot(tspan,ZB);
set(rp(5),'Color',Colors(3,:), 'LineWidth',2,'LineStyle',Style(3));
ylabel('Zwangsbedingung ','FontSize',14)
xlabel('Zeit t in s','FontSize',14)
% legend(rp,'Winkel phi','Abstand rs', 'ZB',...
%        'location','west','numcolumns',1);
% legend box off
grid on
set(gca,'FontSize',16);

%%
% Abrollphase
fig = figure();
hold on
PlotButterSide(rs(1), 0, 0, LB, Colors);
for k= 1:2:length(xo)
   x = [xo(k) xu(k)];
   y = [yo(k) yu(k)];
   line(x,y,'Color',Colors(2,:), 'LineWidth',2,'LineStyle',Style(1))
end
for k= 1:1:length(xo2)
   x2 = [xo2(k) xu2(k)];
   y2 = [yo2(k) yu2(k)];
   if (yo2(k) > -H) && (yu2(k) > -H) 
       line(x2,y2,'Color',Colors(7,:), 'LineWidth',2,'LineStyle',Style(1));
       PlotButterSide(xs2(k), ys2(k), phi2(k), LB, Colors);
   end
end
pl(1) = line([xAu 0],[0 0]);
pl(2) = line([0 0],[0 -LB/2]);
pl(3) = line([-LB/2 0],[-LB/2 -LB/2]);
pl(4) = line([-LB/2 -LB/2],[-LB/2 -H]);
set(pl(1:4),'color',Colors(3,:),'LineWidth',2,'LineStyle',Style(3));
pl(5) = line([xAu xAo],[-H -H]);
set(pl(5),'color',Colors(15,:),'LineWidth',2,'LineStyle',Style(2));
axis equal
ylim([yAu yAo]);
xlim([xAu xAo]);
grid on
set(gca,'FontSize',16);
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


%% Funktionen

function Lend = Drehung(phi,rs,LB)
    Lend(1,:) = -(LB/2 - rs).*cos(phi);
    Lend(2,:) = (LB/2 + rs).*cos(phi);
    Lend(3,:) = (LB/2-rs).*sin(phi);
    Lend(4,:) = -(LB/2+rs).*sin(phi);
end

function PlotButterSide(xs, ys, phi, LB, Colors)
    P = LB/5;
    xP = [xs xs+P*sin(phi)/2];
    yP = [ys ys+P*cos(phi)/2];
    line(xP, yP,'color',Colors(10,:),'LineWidth',6,'LineStyle','-'); 
end

% Brett in Kontaktphase
function dY = dgl_Brett1(t,Y,P1)
    % Y(1)- Winkel phi(t), 
    % Y(2)- Winkelgeschwindigkeit dphi(t)/dt
    % Y(3)- Abstand rs(t), 
    % Y(4)- Geschwidnigkeit drs(t)/dt
    phi   = Y(1);
    dphi  = Y(2);
    rs    = Y(3);
    drs   = Y(4);
    dY    = [dphi;...   
         rs*(P1.g*cos(phi)-2*drs*dphi)/(rs*rs+P1.js);... 
         drs;...
         -rs*dphi*dphi+P1.g*sin(phi)];
end


% Ereignisfunktion

function [value,isterminal,direction] = MyEvent1(t,Y,P1)
    %Ereignisfunktion bis Eintreten des Gleitens bei Slip Winkel thetaS bzw.
    %thetaC
    value = P1.g*cos(Y(1))-2*Y(4).*Y(2); % detect length of line
    isterminal = 1;      % stop the integration
    direction = 0;       % negative direction
end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------

