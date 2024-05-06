% -------------------------------------------------------------------------
% GravityAssistSimulation.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Finger√ºbungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm simuliert die Streuung einer Raumsonde an einem Planeten
% im heliozentrischen KOS
% 
% Beispiel Jupiter,  andere Planeten ausw‰hlbar
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];
Marker = ['o','d','o','s','+'];

%% Parameter
% Variablen, Konstanten, Parameter, hier alles in m, kg, s
% Parameter hier lles in m, kg, s
G   = 6.671e-11;      % G in in m^3/s^2/kg
AE  = 149597870700;   % in m
AEk = 149597870.700;  % in km
mP  = [0.33022; 4.8685; 5.9737; 0.64185; 1898.7; 568.51; ...
       86.849; 102.44]*1e24;  % in kg
aP  = [0.3870893; 0.72333199; 1.000000; 1.52366231; 5.20336301;...
       9.53707032; 19.191264; 30.068964]*AE;  %in m
MS  = 1.989e30;

% Planetendaten ab hier alles in km,s,kg
vPlan= sqrt(G*MS./aP)/1000; % Planetengeschwindigkeiten heliozentrisch 

%AB und Parameter Jupiter
muS  = MS*G*1e-9;      % Sonne
muJ  = mP(5)*G*1e-9;   % Jupiter in km^3/s^2
vJ    = vPlan(5);
RSOI  = 48.2e06;    % in km
aJ    = aP(5)*1e-3; % in km   
omegJ = sqrt(muS/(aJ)^3);  % in 1/s
TJ    = 2*pi/86400/365/omegJ;   % Umlaufzeit in Jahren
x0J   = aJ;
y0J   = 0;
dx0J  = -omegJ*aJ;
dy0J  = +omegJ*aJ;
P1.AB = [x0J dx0J y0J dy0J];  
P1.omegJ = omegJ;
P1.aJ    = aJ;
P1.muJ   = muJ;

%AB und Parameter Sonde
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v1    = 13.0;  % "Einfluggeschwindigkeit" in km/s
gamma1 = deg2rad(75);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0    = aJ-RSOI;
y0    = 0;
dx0   = v1*cos(gamma1);
dy0   = v1*sin(gamma1);
AB = [x0 dx0 y0 dy0, x0J dx0J y0J dy0J];       % AB Sonde & Jupiter

tv    = linspace(0,1.5*RSOI/v1,10001);         % Gesch‰tze Dauer Vorbeiflug

%% Berechnung 
% opt = odeset('AbsTol',1.e-10,'RelTol',1.e-8);
% [tout,Yout]=ode45(@dgl_GravityAssist, tv,  AB, opt, P1);
[tout,Yout] = SatellitenPfad(P1, tv, AB);

% Raumschiffbahn und Geschwindigkeit
xs = Yout(:,1);
ys = Yout(:,3);
dxs = Yout(:,2);
dys = Yout(:,4);
vs  = sqrt(dxs.^2+dys.^2);
% Jupiterbahn
xJ = Yout(:,5);
yJ = Yout(:,7);



%% Graphische Ausgabe

u=linspace(0,360,360);

% Alles normieren auf AE
xJ = xJ/AEk;
yJ = yJ/AEk;
xs = xs/AEk;
ys = ys/AEk;
xRSOI(1,:) = xJ(1) + RSOI/AEk*cosd(u);
yRSOI(1,:) = yJ(1) + RSOI/AEk*sind(u);
xRSOI(2,:) = xJ(end) + RSOI/AEk*cosd(u);
yRSOI(2,:) = yJ(end) + RSOI/AEk*sind(u);

% Bahnen

figure('name','Gravity Assist am Jupiter' )
grid on
ttl = title('Gravity Assist am Jupiter' );
set(ttl,'FontSize',14, 'FontWeight','normal')
xlabel('\it x \rm in AE','FontSize',13); 
ylabel('\it y \rm in AE','FontSize',14);
hold on
axis equal
for k=1:2500:length(xJ)
  p(1) = plot(xJ(k),yJ(k),'o', 'MarkerFaceColor', Colors(5,:),...
         'MarkerEdgeColor',Colors(5,:),'MarkerSize',10);
  p(2) = plot(xs(k),ys(k),'+', 'MarkerFaceColor', Colors(2,:),...
         'MarkerEdgeColor',Colors(2,:),'MarkerSize',8, 'Linewidth',2);
end
plot(xJ,yJ,'Color',Colors(5,:),'Linewidth',2);
plot(xs, ys,'Color',Colors(2,:),'Linewidth',2);
p(3)=plot(xRSOI(1,:),yRSOI(1,:),'Color',Colors(5,:),'Linewidth',1,...
                'LineStyle',Style(4));
plot(xRSOI(2,:),yRSOI(2,:),'Color',Colors(5,:),'Linewidth',1,...
                'LineStyle',Style(4));
grid on
h1=legend(p(1:3),'Jupiter','Raumsonde','SOI',...
                 'location','northeast');
legend box off
set(h1,'FontSize',14)
set(gca,'FontSize',14)

% Geschwindigkeit
figure('name','Geschwindigkeit der Raumsonde')
plot(tout/86400, vs,'Color',Colors(2,:),'Linewidth',2,...
                'LineStyle',Style(1));
hold on
ttl = title('Geschwindigkeit der Raumsonde' );
set(ttl,'FontSize',14, 'FontWeight','normal')
xlabel('\it t \rm in d','FontSize',13); 
ylabel('\it v \rm in km/s','FontSize',14);
grid on
set(gca,'FontSize',14)

% -------------------------------------------------------------------------
% Funktionen
% DGL
function [tout,Yout] = SatellitenPfad(P1, tv, AB)
    opt = odeset('AbsTol',1.e-9,'RelTol',1.e-7,'Events',@events);
    [tout,Yout,te,ye,ie]=ode45(@(t,Y)dgl_GravityAssist(t,Y,P1),tv,AB,opt);
    function [value,isterminal,direction] = events(t,Y)
        r1  = sqrt((Y(1)-Y(5)).^2 + (Y(3)-Y(7)).^2);
        value = (r1 - 49e06);       % Abbruch bei Abstand RSOI, Wert hier 
                                    % eingeben. Jupiter = 48.6e6 km
        isterminal = 1;             % Stop Integration
        direction  = 0;             % Eagl welche Richtung
    end
end

function dY = dgl_GravityAssist(t, Y, P1)
    % Y(1,8)=xS,vSx,yS,vSy, xJ, vxJ, yJ, vyJ
    % Es muss ein Spaltenvektor zur√ºckgegeben werden 
    dY     =  zeros(8,1); 
    omegJ = P1.omegJ;
    aJ    = P1.aJ;
    muJ   = P1.muJ;
    x  = Y(1);
    y  = Y(3);
    xJ = Y(5);
    yJ = Y(7);
    r  = sqrt((x-xJ).^2 + (y-yJ).^2);
    dY(1)  =  Y(2);
    dY(2)  = -muJ*(Y(1)-xJ)./r.^3;    
    dY(3)  =  Y(4);
    dY(4)  = -muJ*(Y(3)-yJ)./r.^3;  
    %Vereinfachte Berechnung Jupiterbahn (=Kreisbahn)
    dY(5)  =  Y(6);
    dY(6)  = -omegJ^2*xJ;    
    dY(7)  =  Y(8);
    dY(8)  = -omegJ^2*yJ;       
end
% Ende Funktionen
% -------------------------------------------------------------------------



