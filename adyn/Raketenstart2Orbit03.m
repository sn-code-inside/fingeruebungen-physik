% -------------------------------------------------------------------------
% RakentenStart2Orbit03.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Astrodynamik" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% 
% Programm berechnet Startazimut und Startfenster
% für Start eines Raumschiffs zur ISS von verschiedenen Weltraumzentren
% 
%--------------------------------------------------------------------------
%%
clc
clear 
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors=GetColorLines;
xAchse1 = [0 30 60 90 120 150 180 210 240 270 300 330 360];
xAchse2 = [0 45 90 135 180 225 270 315 360];
yAchse1 = [0 30 60 90 120 150 180];

%%
% Initialisierung
G   = 6.67430e-11;                 % G in m3 /kg /s2
ME  = 5.98e24;                     % Masse der Erde in kg
RE  = 6.378e6;                     % Erdradius in m
gE  = G*ME/RE^2;                   % gE in m/s^2
omega = 2*pi/86400;                % Erdrotation in rad/s

% Startpunkte
% 1 - Cape Canaveral, 2 - Baikonur, 3 - Kourou
lats = deg2rad([28.383;45.617;5.158]);
lons = deg2rad([-80.6;63.317;-52.65]);
OrtsNamen = ["Cape Canaveral", "Baikonur", "Kourou"];


%% Inklination als Funktion des Startazimuts

beta = deg2rad(linspace(0,360,361));
for k=1:length(lats)
    incl(k,:) = acos(cos(lats(k))*sin(beta(:)));
end

figure('name','Bahnneigung als Funktion des Startazimuts')
hold on
for k=1:length(lats)
  p(k) = plot(rad2deg(beta(:)),rad2deg(incl(k,:)), ...
         'color',Colors(k+1,:),'Linewidth',2);
end
line([0,360],[90,90],'color',Colors(8,:),'Linewidth',1,'lineStyle',':')
str= "Bahnneigung als Funktion des Startazimuts";
grid on
hp1 = title(str,'FontSize',12);
set(hp1,'FontSize',14,'FontWeight','normal'); 
hp2=legend(p,OrtsNamen,'location','northwest'); 
set(hp2,'FontSize',14,'FontWeight','normal'); 
axis([0, 360, 0 180]);
xticks(xAchse2)
yticks(yAchse1)
legend box off;
grid on
set(gca,'FontSize',16);
xlabel('Startazimut \beta in °','FontSize',14); 
ylabel('Bahnneigung in °','FontSize',14)
grid on
set(gca,'FontSize',16);


%% Berechnung Startfenster und Startazimut für ISS

% Inklination ISS
incISS  = deg2rad(51.6);
Omega   = 122.3522;
dtepoch = datetime('10.02.2006');
T  = juliandate(dtepoch);
t  = Jd2JJht(T);

beta1=asin(cos(incISS)./cos(lats));
beta2= pi-beta1;

theta1 = acos(cos(beta1)./sin(incISS));
theta2 = acos(cos(beta2)./sin(incISS));

beta1d  = rad2deg(beta1);
beta2d  = rad2deg(beta2);
theta1d = rad2deg(theta1);
theta2d = rad2deg(theta2);

LWST1 = Omega + theta1d;
LWST2 = Omega + theta2d;

% Startfenster am 10.2.2006

MJ    = MJD(T);
GMST0 = rad2deg(GMST(MJ));
LST   = GMST0 + rad2deg(lons);


figure('name','Startfenster')
hold on
for k=1:length(lons)
   p(k) = plot(rad2deg(lons(k)),(LWST1(k)-GMST0)/15,'o', ...
         'color',Colors(k+1,:),'Linewidth',3);      
   plot(rad2deg(lons(k)),(LWST2(k)-GMST0)/15,'d', ...
         'color',Colors(k+1,:),'Linewidth',3); 
   line([rad2deg(lons(k)) rad2deg(lons(k))],...
        [(LWST1(k)-GMST0)/15 (LWST2(k)-GMST0)/15],...
        'color',Colors(k+1,:),'Linewidth',1, 'linestyle',':');
end
grid on
line([-90,90],[0,0],'color',Colors(3,:),'Linewidth',1,'lineStyle',':')
str= strjoin(["Startfenster am ", datestr(dtepoch)]);
grid on
hp1 = title(str,'FontSize',12);
set(hp1,'FontSize',14,'FontWeight','normal'); 
hp2=legend(p,OrtsNamen,'location','north'); 
set(hp2,'FontSize',14,'FontWeight','normal'); 
axis([-90, +90, -4 16]);
yticks([-4 0 4 8 12 16])
legend box off;
grid on
set(gca,'FontSize',16);
xlabel('Geographische Länge \lambda_0 in °','FontSize',14); 
ylabel('Zeit vs GMST in h','FontSize',14)
grid on
set(gca,'FontSize',16);
