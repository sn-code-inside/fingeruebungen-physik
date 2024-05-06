% -------------------------------------------------------------------------
% LangzeitJahreszeiten.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Berechnet die Dauer der Jahreszeiten und den Energieeintrag der Sonne 
% auf die Nordhemisphäre für einen Zeitraum von 8000 Jahren. 
% Einfache Modellierung der Veränderungen von Exzentrizität, 
% Ekliptikneigung und Richtung Erdachse.
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

% Initialisierung
ups0 = 0;
pi2=2*pi;

dt1 = datetime('-2000-01-01 00:00:00');
T1 = juliandate(dt1); % Julianisches Datum
MJuDa1 = juliandate(dt1,'modifiedjuliandate');% Modif. Julianisches Datum
DatumStr1 =  string(dt1,'dd.MM.yyyy');

dt2 = datetime('6000-01-01 00:00:00');
T2 = juliandate(dt2); % Julianisches Datum  % Bedeckungsdatum
MJuDa2 = juliandate(dt2,'modifiedjuliandate');% Modif. Julianisches Datum
DatumStr2 =  string(dt2,'dd.MM.yyyy');

T=linspace(T1,T2,10000);
TJ=(T-2451545.0)/36525;

% Berechnung der Parameter Exzentrizität, Neigung und Richtung Erdachse
PSun = 3.828e26;  % Strahlungsleistung der Sonne
RE = 6378000;     
AE = 149597870700;

alpha0= (-12.5); 
alpha = wrapToPi(pi/180*(1.7192*TJ+alpha0)); % Richtung Erdachse
eps0= 23.43929111;
eps   = eps0   + (46.815/3600)*TJ; % Neigung Erdachse
ex0   = 0.016710;
ex    = ex0    - 3804e-8*TJ; % Exzentrizität

% Berechnung über analytische Formel 
DeltaT  = 4*ex.*cos(alpha)/pi - 2*ex.^3.*cos(3*alpha)/3/pi; %(3.187)
PSommer  = PSun*RE^2/8/AE^2./sqrt(1-ex.^2).*(2*sind(eps) + ...
    pi)./(pi+4.*ex.*cos(alpha)); %(3.191)

%------------------------------------------------------------------------------
%%  Ausgabe

header1='Langzeitentwicklung 8000 Jahre';
figure('Name',header1);
plot(TJ/10+2, eps-eps0,'LineWidth',2);
hold on;
plot(TJ/10+2, ex*10,'LineWidth',2);
plot(TJ/10+2, alpha,'LineWidth',2,'LineStyle','-.');
lgd=legend('Ekliptikschiefe {\epsilon}-{\epsilon}_{2000} in °','Exzentrizität e x 10', 'Richtung Erdachse {\alpha} in rad','location','southeast');
lgd.FontSize=16;
legend boxoff;
grid on;
xlim([-2  6]);
lgd=xlabel('Jahrtausend');
set(gca,'FontSize',16);

header1='Langzeitentwicklung 8000 Jahre';
figure('Name',header1);
subplot(2,1,1);
plot(TJ/10+2, DeltaT*365,'Color',Colors(2,:),'LineWidth',2,'LineStyle','-.');
hold on;
grid on;
xlim([-2  6]);
ylim([0  10]);
lgd=xlabel('Jahrtausend');
lgd=ylabel('{\Delta} T in Tagen');
set(gca,'FontSize',16);

subplot(2,1,2);
plot(TJ/10+2, PSommer,'Color',Colors(10,:),'LineWidth',2,'LineStyle','-');
hold on;
grid on;
xlim([-2  6]);
ylim([1.06e17  1.08e17]);
lgd=ylabel('P_{Sommer} in W');
lgd=xlabel('Jahrtausend');
set(gca,'FontSize',16);


%% Sehr lange Zeiten
% Zeiten in Jahrhunderten 
T1=-2500;
T2=+2500;
T=linspace(T1,T2,1000);

% Berechnung der Parameter Exzentrizität, Neigung und Richtung Erdachse für
% lange Zeiten
ex   = 0.0315   + 0.0265*sin(pi2*T/1000-2.549) + 0.00*cos(pi2*T/4130+1050);
eps  = 22.8     + 1.7*sin(pi2*T/410+2.76);
alpha = wrapToPi(pi/180*(1.7192*T+alpha0));

% Berechnung über analytische Formel 
DeltaT  = 4*ex.*cos(alpha)/pi - 2*ex.^3.*cos(3*alpha)/3/pi; %(3.187)
PSommer  = PSun*RE^2/8/AE^2./sqrt(1-ex.^2).*(2*sind(eps) + ...
    pi)./(pi+4.*ex.*cos(alpha)); %(3.191)

%------------------------------------------------------------------------------
%%  Ausgabe

header1='Langzeitentwicklung 0.5 Mio Jahre';
figure('Name',header1);
plot(T/10, eps-22.8,'LineWidth',2);
hold on;
plot(T/10, ex*100,'LineWidth',2);
plot(T/10, alpha,'LineWidth',2,'LineStyle','-.');
lgd=legend('Ekliptikschiefe {\epsilon} in ° (um Mittelwert)','Exzentrizität e x 100', 'Richtung Erdachse {\alpha} in °');
lgd.FontSize=16;
legend boxoff;
lgd=xlabel('Jahrtausend');
xlim([-250  250]);
ylim([-4  8]);
set(gca,'FontSize',18);

figure();
subplot(2,1,1);
plot(T/10, DeltaT*365,'Color',Colors(2,:),'LineWidth',2,'LineStyle','-.');
hold on;
xlim([-250  250]);
lgd=xlabel('Jahrtausend');
lgd=ylabel('{\Delta} T in Tagen');
set(gca,'FontSize',18);


subplot(2,1,2);
plot(T/10, PSommer,'Color',Colors(10,:),'LineWidth',2,'LineStyle','-');
hold on;
xlim([-250  250]);
lgd=ylabel('P_{Sommer} in W');
lgd=xlabel('Jahrtausend');
set(gca,'FontSize',18);

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------
