% -------------------------------------------------------------------------
% Mondbahn02.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet die heliozentrische Mondbahn aus den Ephemeriden
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

% Initialisierung
SunData = PerturbImport('SonnePos.csv');

%% Einfache koplanare Modellrechnung
%
AE     = 149.7;        % AE in Mio km

%Fall  aP*OmegP^2  > aM*omegM^2 (Erde-Mond)
TP     = 365.25;        % in Tagen
TM     = 29.53;
aP     = 1;             % Halbachse Planet in AE
aM     = 0.384/AE;      % Halbachse Mond in AE
ratio  = (aP/TP^2)/(aM/TM^2);   % muss > 1 sein
fprintf('\n Verhältnis (aP/TP^2)/(aM/TM^2) = %4.2f \n',ratio);

% Berechnung
OmegP  = 2*pi/TP;    % Omega Planet
omegM  = 2*pi/TM;    % omega Mond
t      = linspace(0,2*pi/min(OmegP,omegM),1000);


xP     = aP*cos(OmegP*t);
yP     = aP*sin(OmegP*t);
xM     = xP+aM*cos(omegM*t);
rP     = sqrt(xP.^2+yP.^2);
yM     = yP+aM*sin(omegM*t);
rM     = sqrt(xM.^2+yM.^2);

% Ausgabe 
figure('Name',' Mondbahn zu ZK');
subplot(1,2,1)
title('Abstand ZK über Zeit');
hold on
h(1) = plot(t,rP,'color',Colors(3,:),'LineWidth',2,'LineStyle',':');
h(2) = plot(t,rM,'color',Colors(10,:),'LineWidth',2,'LineStyle','-');
grid on;
xlabel('Zeit');
% ylim([0 inf]);
legend(h,'Planetenbahn', 'Mondbahn', 'location', 'south');
legend box off
ylim([0.95*min(rM) 1.05*max(rM)]);
ylabel('Abstand Mond-ZK in km');
set(gca,'XGrid','on', 'YGrid', 'on', 'Fontsize', 14, 'linewidth', 1);

subplot(1,2,2)
hold on
title('Bahnebene koplanare Kreisbahnen')
h(1) = plot(xP,yP,'color',Colors(3,:),'LineWidth',2,'LineStyle',':');
h(2) = plot(xM,yM,'color',Colors(10,:),'LineWidth',1,'LineStyle','-');
grid on;
axis equal
axis square
PlotCircle (0,0,0.05,Colors(10,:),3);  %große Masse
legend(h,'Planetenbahn', 'Mondbahn', 'location', 'south');
legend box off
xlim([-1.2 1.2]);
ylim([-1.2 1.2]);
xlabel('x in AE');
ylabel('y in AE');
set(gca,'XGrid','on', 'YGrid', 'on', 'Fontsize', 14, 'linewidth', 1);




%% Berechnung für Erde aus Ephemeriden-Daten

% Vorbereitung Felder
Jahr = 2000;
dt0  = datetime(Jahr,01,01);    % Startdatum ET
T0   = juliandate(dt0);
epsE = deg2rad(EpsErde(T0));
t0   = Jd2JJht(T0);             % Zeit in Jul. Jahrhunderten seit JD2000
delt = 1/36525/4;               % 6 stundenweise Berechnung
Nend = 365*4;                   % Stützpunkte eines Jahres
% Schleife t
for k=1:Nend
  t(k)  = t0 + (k-1)*delt; 
end
rMass = 'km';
for k=1:Nend
  SunPos = SonneExakt(t(k),SunData, epsE, rMass);
  xE(k)  = SunPos.xyz(1);
  yE(k)  = SunPos.xyz(2);
  zE(k)  = SunPos.xyz(3);
  rE(k)  = sqrt(xE(k).^2+yE(k).^2+zE(k).^2);
  MoonPos= MondExakt(t(k),epsE,rMass);
  xM(k)  = MoonPos.xyz(1);
  yM(k)  = MoonPos.xyz(2);
  zM(k)  = MoonPos.xyz(3);
  XM(k)  = xE(k) + xM(k);
  YM(k)  = yE(k) + yM(k);
  ZM(k)  = zE(k) + zM(k);
  RM0(k)  = sqrt(XM(k).^2+YM(k).^2+1*ZM(k).^2);
  RM1(k)  = sqrt(XM(k).^2+YM(k).^2+0*ZM(k).^2);
end
AE = 149597870; % in km

% Ausgabe 
figure('Name',' Mondbahn');
subplot(1,2,1)
title('Abstand Mohn von Sonne')
hold on
taxis = datetime(JJht2Jd(t),'ConvertFrom','juliandate');
h(1) = plot(taxis,rE,"DatetimeTickFormat","MM-dd",'color',Colors(3,:),...
     'LineWidth',2,'LineStyle',':');
h(2) = plot(taxis,RM1,"DatetimeTickFormat","MM-dd",'color',Colors(10,:),...
    'LineWidth',2,'LineStyle','-');
grid on;
legend(h,'Planetenbahn', 'Mondbahn', 'location', 'south');
legend box off
xlabel('Zeit');
% ylim([0 inf]);
ylabel('Abstand Mond-Sonne in km');
set(gca,'XGrid','on', 'YGrid', 'on', 'Fontsize', 14, 'linewidth', 1);

subplot(1,2,2)
title('Mondbahn projiziert auf Ekliptik')
hold on
h1(1) = plot(xE/AE,yE/AE,'color',Colors(3,:),'LineWidth',2,'LineStyle',':');
h1(2) = plot(XM/AE,YM/AE,'color',Colors(10,:),'LineWidth',1,'LineStyle','-');
grid on;
axis equal
axis square
PlotCircle (0,0,0.05,Colors(10,:),3);  % große Masse
legend(h1,'Planetenbahn', 'Mondbahn', 'location', 'south');
legend box off
xlim([-1.2 1.2]);
ylim([-1.2 1.2]);
xlabel('x in AE');
ylabel('y in AE');
set(gca,'XGrid','on', 'YGrid', 'on', 'Fontsize', 14, 'linewidth', 1);


% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------
