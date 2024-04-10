% -------------------------------------------------------------------------
% PlanetenbahnenKepler.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet die Bahnen der Planeten Venus bis Jupiter
% auf Basis der Keplergleichung.
% Man berechnet für den 22.11.2065 die Venus-Jupiter-Konjunktion.
% (Äquinoktium und Ekliptik des Datum) 
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

% Initialisierung
AE = 149597870;
dt1 = datetime('2065-01-01 00:00:00');
T1 = juliandate(dt1); % Julianisches Datum
MJuDa1 = juliandate(dt1,'modifiedjuliandate');% Modif. Julianisches Datum
DatumStr1 =  string(dt1,'dd.MM.yyyy');

dt2 = datetime('2065-11-22 00:00:00');
T2 = juliandate(dt2); % Julianisches Datum  %Bedeckungsdatum
MJuDa2 = juliandate(dt2,'modifiedjuliandate');% Modif. Julianisches Datum
DatumStr2 =  string(dt2,'dd.MM.yyyy');

dt3 = datetime('2066-01-01 00:00:00');
T3 = juliandate(dt3); % Julianisches Datum
MJuDa3 = juliandate(dt3,'modifiedjuliandate');% Modif. Julianisches Datum
DatumStr3 =  string(dt3,'dd.MM.yyyy');

Aequi = 'Datum';
Aequi = 'J2000';
[BaPa,BaPadot]=OrbitParameter(T2, Aequi);
%-------------------------------------------------------------------------
% Begin Rechnung
maxP=1000;
T_vector=linspace(T1,T2,maxP);
T_vector2=linspace(T2,T3,maxP);
% Berechnung nach Keplerlösung

for k=2:5
    Planets(k) =PlanetPQR(T_vector, BaPa, BaPadot, k);
    Planets2(k)=PlanetPQR(T_vector2, BaPa, BaPadot, k);
end
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graphische Ausgabe

header1='Jupiterbedeckung durch Venus 2065';
figure('Name',header1);
for iPlot = 2:size(Planets,2)
    plot(Planets(iPlot).xyz(1,:),Planets(iPlot).xyz(2,:),'Color', Colors(iPlot,:),'LineWidth',2);
    hold on
    axis equal
end
%Sonne FP
SonneFP;
%Endpunkte
p(12)=plot([Planets(5).xyz(1,maxP) Planets(3).xyz(1,maxP)], [Planets(5).xyz(2,maxP) Planets(3).xyz(2,maxP)],':','Color', Colors(3,:));
p(12).LineWidth=2;
plot(50,50,'Color','w');  %Leerzeile für Legende
%Positionen Anfangs, Endpunkte, Midpoint
for iPlot = 2:5
    p=plot(Planets(iPlot).xyz(1,1),Planets(iPlot).xyz(2,1),'Color', Colors(iPlot,:));
    p.Marker='+';
    p.LineWidth=2;
    txt= Planets(iPlot).Name;
    text(Planets(iPlot).xyz(1,1)+0.05,Planets(iPlot).xyz(2,1)+0.05,txt,'FontSize',16,'Color',Colors(iPlot,:));
    p=plot(Planets(iPlot).xyz(1,maxP),Planets(iPlot).xyz(2,maxP),'Color', Colors(iPlot,:));
    p.Marker='o';
    p.LineWidth=2;
    p=plot(Planets2(iPlot).xyz(1,maxP),Planets2(iPlot).xyz(2,maxP),'Color', Colors(iPlot,:));
    p.Marker='*';
    p.LineWidth=2;
   
 end
for iPlot = 2:size(Planets,2)
    plot(Planets2(iPlot).xyz(1,:),Planets2(iPlot).xyz(2,:),'Color', Colors(iPlot,:),'LineStyle','--','LineWidth',2);
end
%Linien Anfangs, Endpunkte, Midpoint
for iPlot = 2:size(Planets,2)
    plot([0 Planets(iPlot).xyz(1,1)], [0 Planets(iPlot).xyz(2,1)],':','Color', Colors(iPlot,:));
    plot([0 Planets(iPlot).xyz(1,maxP)], [0 Planets(iPlot).xyz(2,maxP)],'--','Color', Colors(iPlot,:));
    plot([0 Planets2(iPlot).xyz(1,end)], [0 Planets2(iPlot).xyz(2,end)],'-.','Color', Colors(iPlot,:));
 end
ylim([-5 2])
xlim([-5 2]);
grid on;
% grid minor,
header2 = strjoin([DatumStr1,' - ',DatumStr3]);
text(-1.5,1.75,header2,'FontSize',18);
xlabel('x in AE')
ylabel('y in AE');
legend(Planets.Name,'Sonne','Richtung Frühlingspunkt','Venus-Jupiter-Bedeckung',DatumStr2,'location','northwest');
legend boxoff;
set(gca,'FontSize',16);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  3-Dimensionale Darstellung
figure('Name',header1);
for iPlot = 2:size(Planets,2)
    p=plot3(Planets(iPlot).xyz(1,:),Planets(iPlot).xyz(2,:),Planets(iPlot).xyz(3,:)*10,'Color', Colors(iPlot,:));
    p.LineWidth=2;
    hold on
    axis equal
end
SonneFP;
p(12)=plot3([Planets(5).xyz(1,maxP) Planets(3).xyz(1,maxP)], [Planets(5).xyz(2,maxP) Planets(3).xyz(2,maxP)], [10*Planets(5).xyz(3,maxP) 10*Planets(3).xyz(3,maxP)],':','Color', Colors(3,:));
p(12).LineWidth=2;
plot(50,50,'Color','w');  %Leerzeile für Legende
for iPlot = 2:size(Planets,2)
    p=plot3(Planets2(iPlot).xyz(1,:),Planets2(iPlot).xyz(2,:),Planets2(iPlot).xyz(3,:)*10,'Color', Colors(iPlot,:),'LineStyle','--');
    p.LineWidth=2;
    hold on
    axis equal
end
for iPlot = 2:5
    p=plot3(Planets(iPlot).xyz(1,1),Planets(iPlot).xyz(2,1),10*Planets(iPlot).xyz(3,1),'Color', Colors(iPlot,:));
    p.Marker='+';
    p.LineWidth=1;
    p=plot3(Planets(iPlot).xyz(1,maxP),Planets(iPlot).xyz(2,maxP),10*Planets(iPlot).xyz(3,maxP),'Color', Colors(iPlot,:));
    p.Marker='o';
    p.LineWidth=2;
    p=plot3(Planets2(iPlot).xyz(1,maxP),Planets2(iPlot).xyz(2,maxP),10*Planets2(iPlot).xyz(3,maxP),'Color', Colors(iPlot,:));
    p.Marker='*';
    p.LineWidth=1;
   
 end
for iPlot = 2:size(Planets,2)
    plot3(Planets2(iPlot).xyz(1,:),Planets2(iPlot).xyz(2,:),10*Planets2(iPlot).xyz(3,:),'Color', Colors(iPlot,:),'LineStyle','--');
end
for iPlot = 2:size(Planets,2)
    plot3([0 Planets(iPlot).xyz(1,1)], [0 Planets(iPlot).xyz(2,1)],[0 10*Planets(iPlot).xyz(3,1)],':','Color', Colors(iPlot,:));
    plot3([0 Planets(iPlot).xyz(1,maxP)], [0 Planets(iPlot).xyz(2,maxP)],[0 10*Planets(iPlot).xyz(3,maxP)],'--','Color', Colors(iPlot,:));
    plot3([0 Planets2(iPlot).xyz(1,end)], [0 Planets2(iPlot).xyz(2,end)],[0 10*Planets2(iPlot).xyz(3,end)],'-.','Color', Colors(iPlot,:));
 end
xlim([-5 2]);
ylim([-5 2]);
xl = xlim;
yl = ylim;
[X,Y] = meshgrid(xl,yl);
surf(X,Y,zeros(size(X)))
shading flat;
alpha 0.05;
grid on;
text(-4.5,-3,0,header2,'fontsize',16);
xlabel('x in AE')
ylabel('y in AE');
zlabel('z in 0.1 AE');
legend(Planets.Name,'Sonne','Richtung Frühlingspunkt','Venus-Jupiter-Bedeckung',DatumStr2,'location','northeast');
legend boxoff;
set(gca,'FontSize',16);

% Ende Programmteil I
%_________________________________________________________________________
% Genauere Berechnung Bedeckung

dt1 = datetime('2065-11-22 12:00:00');
T1 = juliandate(dt1); % Julianisches Datum
MJuDa1 = juliandate(dt1,'modifiedjuliandate');% Modif. Julianisches Datum
DatumStr1 =  string(dt1,'dd.MM.yyyy');

dt2 = datetime('2065-11-22 14:00:00');
T2 = juliandate(dt2); % Julianisches Datum  %Bedeckungsdatum
MJuDa2 = juliandate(dt2,'modifiedjuliandate');% Modif. Julianisches Datum
DatumStr2 =  string(dt2,'dd.MM.yyyy');

%Bestimmung der Schrittweite zu Zeitberechnung
St_p_H=60;
if day(dt1) == day(dt2)
    Nhour=hour(dt2)-hour(dt1);
else
    Nhour = 24-hour(dt1) + hour(dt2);
end
maxP=Nhour*St_p_H+1; 
Lab_p_H=4;  % jede Stunde/Lab_p_H 4 Ticks setzen auf der Bahn
NrLabSteps = St_p_H/Lab_p_H;
NrLabels=Nhour*Lab_p_H+1;


T_vector=linspace(T1,T2,maxP);
MJD_vector=linspace(MJuDa1,MJuDa2,maxP);
t_vector = Jd2JJht(T_vector);

% Lichtlaufzeit
tau = 8.32/1440/36525; %Lichtlaufzeit Sonne-Erde in JJht;
tauJup = Planets(5).ekl(1,1)*tau; %Lichtlaufzeit Jupiter-Erde;
tauVen = Planets(2).ekl(1,1)*tau; %Lichtlaufzeit Venus-Erde;
epsE = deg2rad(EpsErde(T1));

% Berechnung nach Störungstheoriekoordinaten
VenData = PerturbImport('Data/VenusPos.csv');
SunData = PerturbImport('Data/SonnePos.csv');
JupData = PerturbImport('Data/JupiterPos.csv');
VenPos  = VenusExakt(t_vector-tauVen,VenData);
JupPos  = JupiterExakt(t_vector-tauJup,JupData);
SunPos  = SonneExakt(t_vector-tau,SunData,epsE,'AE');

VenPos.geo  = SunPos.xyz + VenPos.xyz;
JupPos.geo  = SunPos.xyz + JupPos.xyz;
VenPos.equ  = CalcAnglesfromXYZ((mtimes(R_x(-epsE),VenPos.geo)));
JupPos.equ  = CalcAnglesfromXYZ((mtimes(R_x(-epsE),JupPos.geo)));

    
    
% Ausdrucken
 for k=1:maxP
    DateObs = datetime(MJD_vector(k)+0.000001,'ConvertFrom','modifiedjuliandate');
    Datestr = string(DateObs, 'dd-MM-yyyy HH:mm:ss ');
    aV = StrHMS(24+rad2deg(VenPos.equ(2,k))/15);
    dV = StrDMS(rad2deg(VenPos.equ(3,k)));
    rJ = num2str(JupPos.equ(1,k),'%5.3f');
    aJ = StrHMS(24+rad2deg(JupPos.equ(2,k))/15);
    dJ = StrDMS(rad2deg(JupPos.equ(3,k)));
    fprintf('\n %s  RA Venus: %s Jupiter: %s   DEC Venus: %s Jupiter: %s   Abstand Jupiter: %s ',...
       Datestr, aV, aJ, dV, dJ, rJ);
end

%__________________________________________________________________________
% Berechnung der relativen Koordinaten vor der Sonnescheibe
% in RA und Deklination (in °)
for k=1:maxP
    DateObs = datetime(MJD_vector(k),'ConvertFrom','modifiedjuliandate');
    Datestr = string(DateObs, 'dd-MM-yyyy HH:mm ');
    %Umrechnung in relative äquatoriale Koordinaten 
    RAV(k)  = rad2deg(VenPos.equ(2,k)-JupPos.equ(2,k));
    DekV(k) = rad2deg(VenPos.equ(3,k)-JupPos.equ(3,k));
end
figure('Name','Jupiter Bedeckung durch Venus');
[RAmin, RAIndex] = min(abs(RAV));
abstJ = (Planets2(5).xyz(1,1)-Planets2(3).xyz(1,1))^2 + ...
        (Planets2(5).xyz(1,2)-Planets2(3).xyz(1,2))^2 + ...
        (Planets2(5).xyz(1,1)-Planets2(3).xyz(1,1))^2;
abstJ = sqrt(abstJ);
abstV = (Planets2(2).xyz(1,1)-Planets2(3).xyz(1,1))^2 + ...
        (Planets2(2).xyz(1,2)-Planets2(3).xyz(1,2))^2 + ...
        (Planets2(2).xyz(1,1)-Planets2(3).xyz(1,1))^2;
abstV = sqrt(abstV);
% Berechnen der scheinbaren Durchmesser Jupiter
rJ = rad2deg(71398/abstJ/AE);
rV = rad2deg(6051/abstV/AE);
plot(RAV, DekV,'LineWidth',2,'Color',Colors(2,:)); %Bahn der Venus
hold on;
PlotCircle(0,0,rJ,Colors(5,:),3); %Jupiterscheibe
PlotCircle(RAV(RAIndex),DekV(RAIndex),rV,Colors(2,:),3); %Venusscheibe
% Uhrzeitlabels
plot(RAV, DekV,'-+','MarkerIndices',1:NrLabSteps:length(DekV),'MarkerSize',8,'LineWidth',1,'Color',Colors(2,:));
for k=1:NrLabels-2
    mylabels(k,:)=string(datetime(T_vector((k-1)*NrLabSteps +1),'convertfrom','juliandate'),'HH:mm');
    RAVL(k)=RAV((k-1)*NrLabSteps +1);
    DekVL(k)=DekV((k-1)*NrLabSteps +1); 
end
h=LabelPoints(RAVL, DekVL,mylabels,'S',1,1,'FontSize',14,'Color',Colors(2,:));
header2 = strjoin(['Jupiter Bedeckung durch Venus am ',DatumStr1," (UT)"]);
title(header2,'FontSize',12);
legend('Venus Bahn','location','northeast');
legend box off;
grid on;
% grid minor,
xlabel('{\alpha} in °');
ylabel('{\delta} in °');
axis equal;
ylim([-0.04 0.04]);
xlim([-0.04 0.04]);
set(gca,'FontSize',16);

% Berechnung von Alt-Az-Koordinaten und Graphische Darstellung für
% bestimmte Location
lambda = -80;
phi    = 50;
ttlstr = "Jupiterbedeckung durch Venus 22.11.2065 ";
ttlstr = sprintf('%s\n Breite: %+4.1f° Laenge: %+5.1f°',ttlstr,phi,lambda);
phi     = deg2rad(phi);
lambada = deg2rad(lambda);

SunPos.Tau = pi*2*LMST(MJD_vector,lambda)/24-SunPos.equ(2);
SunPos.Alt = asind(sin(phi).*sin(SunPos.equ(3))+cos(phi).*cos(SunPos.equ(3)).*cos(SunPos.Tau));
VenPos.Tau = pi*2*LMST(MJD_vector,lambda)/24-VenPos.equ(2);
JupPos.Tau = pi*2*LMST(MJD_vector,lambda)/24-JupPos.equ(2);
VenPos.Alt= asind(sin(phi).*sin(VenPos.equ(3))+cos(phi).*cos(VenPos.equ(3)).*cos(VenPos.Tau));
JupPos.Alt= asind(sin(phi).*sin(JupPos.equ(3))+cos(phi).*cos(JupPos.equ(3)).*cos(VenPos.Tau));

SunPos.Az  = atan2d(cos(SunPos.equ(3)).*sin(SunPos.Tau),(cos(SunPos.equ(3)).*cos(SunPos.Tau).*sin(phi)-sin(SunPos.equ(3))*cos(phi)))-180;
VenPos.Az = atan2d(cos(VenPos.equ(3)).*sin(VenPos.Tau),(cos(VenPos.equ(3)).*cos(VenPos.Tau).*sin(phi)-sin(VenPos.equ(3))*cos(phi)))-180;
JupPos.Az = atan2d(cos(JupPos.equ(3)).*sin(JupPos.Tau),(cos(JupPos.equ(3)).*cos(JupPos.Tau).*sin(phi)-sin(JupPos.equ(3))*cos(phi)))-180;
figure('Name',header1);
plot(datetime(T_vector,'convertfrom','juliandate'), SunPos.Alt,'-o','MarkerIndices',1:NrLabSteps:length(SunPos.Az),'MarkerSize',8,'LineWidth',1,'Color',Colors(10,:));
hold on;
plot(datetime(T_vector,'convertfrom','juliandate'), VenPos.Alt,'-+','MarkerIndices',1:NrLabSteps:length(SunPos.Az),'MarkerSize',8,'LineWidth',1,'Color',Colors(2,:));
plot(datetime(T_vector,'convertfrom','juliandate'), JupPos.Alt,'-d','MarkerIndices',1:NrLabSteps:length(SunPos.Az),'MarkerSize',8,'LineWidth',1,'Color',Colors(5,:));
plot(datetime(T_vector,'convertfrom','juliandate'), 0*T_vector, 'LineStyle','-.','Color','b','LineWidth',2);
Centerliney = linspace(-30,30,maxP);
Centerlinex = ones(maxP);
plot(datetime(T_vector,'convertfrom','juliandate'), 0*T_vector, 'LineStyle','-.','Color','b','LineWidth',2);
plot(datetime(T_vector(RAIndex)*Centerlinex,'convertfrom','juliandate'), Centerliney, 'LineStyle',':','Color','b','LineWidth',1);
xlim([datetime(T_vector(1),'convertfrom','juliandate') datetime(T_vector(maxP),'convertfrom','juliandate')]); 
if VenPos.Alt(1) < VenPos.Alt(maxP)
    ylim([VenPos.Alt(1) VenPos.Alt(maxP)]);
else
    ylim([VenPos.Alt(maxP) VenPos.Alt(1)]);
end
xlabel('Zeit in UT');
ylabel('Hoehe in °');

title(ttlstr); 
legend('Sonne', 'Venus', 'Jupiter', 'Location', 'south');
legend box off;
grid on;
set(gca,'FontSize',16);

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------