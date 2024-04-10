% -------------------------------------------------------------------------
% MondErdeSonneApollo8.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet die Positionen von Mond, Erde, Sonne und Apollo 8 am
% 24.12.1968
% (Äquinoktium und Ekliptik J2000) 
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

% Zeitraum und Darstellung 
dt1 = datetime('1968-12-19 16:38:00');
dt2 = datetime('1968-12-27 16:38:00');
dt3 = datetime('1968-12-24 16:38:00');

T1 = juliandate(dt1);
T2 = juliandate(dt2);
T3 = juliandate(dt3);

StepsperDay=24; % Zahl der Stützstellen pro Tag
Nges=(T2-T1)*StepsperDay; % Zahl der Stunden im Intervall T1 - T2
% Input
StepsMarker= 1; % Zahl der markierten Tage 
NrStepsMarker=StepsMarker*StepsperDay; % Zahl der markierten Tage
% StepsLabel=2; % Zahl der angezeigten Tage
% NrStepsLabel=StepsLabel*NrStepsMarker; % Zahl der angezeigten Tage
% MaxStepsLabel=floor(Nges/NrStepsLabel)+1;

T=linspace(T1,T2,Nges); % Berechnungsintervall

DatumStr1 =  string(dt1,'dd.MM.yy');
DatumStr2 =  string(dt2,'dd.MM.yy');

Aequi = 'J2000';
[BaPa,BaPadot]=OrbitParameter(T1, Aequi);

%%
%-------------------------------------------------------------------------
% Begin Rechnung

epsErad = deg2rad(EpsErde(T1));
% Berechnung der Erdeposition zum Perihel des Asteroiden

% Berechnung der geoozentrischen Koordinaten von Sonne und Mond
% von T1 bis T2 
Erde = PlanetPQR(T, BaPa, BaPadot,3);
Moon = KeplerMond(T, epsErad);

%Normierung auf mittleren Mondabstand
D_Moon = 385000;
AE = 149550000;

Moon.xyz = +Moon.equ_xyz/D_Moon;
Erde.xyz = +Erde.xyz*AE/D_Moon;
MPEX = Erde.xyz(1,1);
MPEY = Erde.xyz(2,1);
Sun.xyz  = -Erde.xyz;


%------------------------------------------------------------------------------
%% Graphische Ausgabe
% x-y-Ebene

header1='Erde-Sonne-Mond-Apollo 8';
figure('Name',header1);
Ax_L = 1;
subplot(1,2,1);  % Achtung um 90° gedreht, damit zu Internet-Bild passend !
plot(-Sun.xyz(2,:),Sun.xyz(1,:),':+','MarkerIndices',...
     1:NrStepsMarker:length(Moon.xyz(1,:)),...
     'LineWidth',1,'Color',Colors(10,:));
hold on;
plot(-Moon.xyz(2,:),Moon.xyz(1,:),':+','MarkerIndices',...
     1:NrStepsMarker:length(Moon.xyz(1,:)),...
     'LineWidth',1,'Color',Colors(2,:));
plot(0,0,'o','Color',Colors(3,:));
plot(0,0,'o','Color',Colors(4,:));
y=linspace(0,100,100);
x=0*y;
p2=plot(x,y,'Color','black','LineStyle', '--','LineWidth',1); %FP
PlotCircle(0,0,6370/D_Moon,Colors(3,:),1);
PlotCircle (-Moon.xyz(2,(day(dt3)-day(dt1))*StepsperDay+1),...
    Moon.xyz(1,(day(dt3)-day(dt1))*StepsperDay+1),...
    (200+1750)/D_Moon,Colors(4,:),2); 
x = [-Moon.xyz(2,(day(dt3)-day(dt1))*StepsperDay+1) 0];
y = [Moon.xyz(1,(day(dt3)-day(dt1))*StepsperDay+1) 0];
line(x,y,'Color',Colors(4,:),'LineStyle',':', 'LineWidth',2);
for k=1:StepsperDay:Nges
    x = [0 -Sun.xyz(2,k)];
    y = [0 Sun.xyz(1,k)];
    line(x,y,'Color',Colors(10,:),'LineStyle','-');
    text(-Moon.xyz(2,k)+Ax_L/100,Moon.xyz(1,k)+Ax_L/30,...
    string(datetime(T(k),'ConvertFrom','juliandate'),'dd.MM.'),...
    'FontSize',16);
end
 
% Housekeeping
title(header1);
xlim([-0.5*Ax_L 1.5*Ax_L]);
ylim([-0.5*Ax_L 1.5*Ax_L]);
axis square;
grid on;
% grid minor;
xlabel('x in D_{Mond}')
ylabel('y in D_{Mond}');

lgdtxt(1) = "Sonne ";
lgdtxt(2) = "Mond";
lgdtxt(3) = "Erde";
lgdtxt(4) = "FP";
legend(lgdtxt(1),lgdtxt(2),lgdtxt(3),'Apollo 8',lgdtxt(4));
legend box off;
set(gca,'FontSize',18);

%________________________________________________________________________
%
% x-z-Ebene
subplot(1,2,2);
plot(Sun.xyz(1,:),Sun.xyz(3,:),':+','MarkerIndices',...
     1:NrStepsMarker:length(Moon.xyz(1,:)),'LineWidth',1,...
     'Color',Colors(10,:));
hold on;
plot(Moon.xyz(1,:),Moon.xyz(3,:),':+','MarkerIndices',...
     1:NrStepsMarker:length(Moon.xyz(1,:)),'LineWidth',1,...
     'Color',Colors(2,:));
PlotCircle(0,0,6370/D_Moon,Colors(3,:),1);
SonneFP;
PlotCircle(0,0,6370/D_Moon,Colors(3,:),1);
for k=1:StepsperDay:Nges
     text(Moon.xyz(1,k)+Ax_L/20,Moon.xyz(3,k),string(datetime(T(k),...
         'ConvertFrom','juliandate'),'dd.MM.'));
end
PlotCircle (Moon.xyz(1,(day(dt3)-day(dt1))*StepsperDay),...
            Moon.xyz(3,(day(dt3)-day(dt1))*StepsperDay),...
            (200+1750)/D_Moon,Colors(4,:),2); 

% Housekeeping
title(header1);
xlim([-0.5*Ax_L 1.5*Ax_L]);
ylim([-Ax_L Ax_L]);
axis square;
grid on;
grid minor;
xlabel('x in D_{Mond}')
ylabel('y in D_{Mond}');
legend(lgdtxt(1),lgdtxt(2),lgdtxt(3));
legend box off;
set(gca,'FontSize',14);

%_________________________________________________________________
%%
% 3-Dimensionale Darstellung
figure('Name',header1);
plot3(Sun.xyz(1,:),Sun.xyz(2,:),Sun.xyz(3,:),':+','MarkerIndices',...
      1:NrStepsMarker:length(Moon.xyz(1,:)),'LineWidth',1,...
      'Color',Colors(10,:));
xlim([-0.5*Ax_L 1.5*Ax_L]);
ylim([-1.5*Ax_L 0.5*Ax_L]);
zlim([-Ax_L Ax_L]);
axis square;
hold on;
plot3(Moon.xyz(1,:),Moon.xyz(2,:),Moon.xyz(3,:),':+','MarkerIndices',...
      1:NrStepsMarker:length(Moon.xyz(1,:)),'LineWidth',1,...
      'Color',Colors(2,:));
PlotCircle(0,0,6370/D_Moon,Colors(3,:),1);
SonneFP;
PlotCircle(0,0,6370/D_Moon,Colors(3,:),1);
plot3(0,0,0,'o','Color',Colors(3,:),'LineWidth',3);
for k=1:StepsperDay:Nges
    x = [0 Sun.xyz(1,k)];
    y = [0 Sun.xyz(2,k)];
    z = [0,Sun.xyz(3,k)];
    line(x,y,z, 'Color',Colors(10,:),'LineStyle','-');
    text(Moon.xyz(1,k)+Ax_L/20,Moon.xyz(2,k),Moon.xyz(3,k),...
         string(datetime(T(k),'ConvertFrom','juliandate'),'dd.MM.'),...
         'FontSize',12);
end
% Housekeeping
title(header1);
xl = xlim;
yl = ylim;
[X,Y] = meshgrid(xl,yl);
surf(X,Y,zeros(size(X)));
shading flat
alpha 0.1
grid on;
grid minor;
xlabel('x in D_{Mond}')
ylabel('y in D_{Mond}');
zlabel('z in D_Mond');
legend(lgdtxt(1),lgdtxt(2),lgdtxt(3),'location','southeast');
legend box off;
set(gca,'FontSize',14);

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------