% -------------------------------------------------------------------------
% SL9numerisch.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Komet Shoemaker-Levy 9 numerische Integration
% Berechnet die Bahn des Kometen Shoemaker-Levy 9 bis zu seinem 
% Zusammenprall mit dem Jupiter über numerische Integration.
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;


% Initialisierung
BaPaK=KometParameter('Kometen.csv');

% Zeiten
dt1 = datetime('1993-01-01 00:00:00');
dt2 = datetime('1994-07-17 20:35:00');
dt3 = datetime('1994-12-31 00:00:00');
T1 = juliandate(dt1);
T2 = juliandate(dt2);
T3 = juliandate(dt3);

Aequi = 'J2000';
% Einlesen Bahnparameter und deren Ableitung 
[BaPa,BaPadot]=OrbitParameter(T1, Aequi);

% Parameter Jupiter
secperday = 86400;          % Sekunden pro Tag
muG      = 2.95479E-04;     % µG [AE^3/d^2]
muJ      = muG/1047.348644; % Jupitermasse*G
AE       = 149597870 ;      % AE in km
RJ       = 69911;           % Radius in km
rJ       = RJ/AE;           % Radius in AE

%-------------------------------------------------------------------------

%Beginn Rechnung Auswahl der Kometen
NrK=5;  %SL9
NrP=5;  %Jupiter

% Position und Geschwindigkeit des Kometen zur Zeit T1 (JPL)
% Ekliptikale Koordinaten T1  1993-Jan-01 00:00
l1=181.1041;
b1=-0.5639;
r1=5.478189665071;
% Ekliptikale Koordinaten T1+1  1993-Jan-02 00:00
l2 = 181.1747;
b2 =  -0.5688;
r2 =   5.478141861326;
% Kartesische Koordinaten
xs1 = CalcXYZfromAngles([r1;deg2rad(l1);deg2rad(b1)]);
xs2 = CalcXYZfromAngles([r2;deg2rad(l2);deg2rad(b2)]);

% PosJT = JupiterKoord(T);
T=linspace(T1,T2,1000);
PosJT = @(T)JupiterKoord_ex(T, 5, BaPa, BaPadot);
JupPos = JupiterKoord_ex(T, 5, BaPa, BaPadot);

% Startwerte für Differentialgleichung
fac=1;
start1 = [xs1(1);xs1(2);xs1(3);fac*(xs2(1)-xs1(1));fac*(xs2(2)-xs1(2));...
    fac*(xs2(3)-xs1(3))];
disp(start1);
% Verschiedene Startwerte bei Geschwindigkeiten
fac=0.8;
start2 = [xs1(1);xs1(2);xs1(3);fac*(xs2(1)-xs1(1));fac*(xs2(2)-xs1(2));...
    fac*(xs2(3)-xs1(3))];
fac=1.2;
start3 = [xs1(1);xs1(2);xs1(3);fac*(xs2(1)-xs1(1));fac*(xs2(2)-xs1(2));...
    fac*(xs2(3)-xs1(3))];

kurve(1) = " Messwerte 1.1.93";
kurve(2) = " 80%  v_0";
kurve(3) = " 120% v_0";
kurve(4) = " 133% Jupitermasse";
kurve(5) = " 33%  Jupitermasse";

%Berechnung der Kraftunterschiede Sonne vs Jupiter bei T1
xJ0 = JupiterKoord_ex(T1,NrP,BaPa,BaPadot);
F1  = - muG/norm(xs1)^3;
F2  = - muJ/norm(xs1-xJ0)^3;
diffF = F2/F1;

% Berechnung der Planetenbahnen bis Kollision und Ende 1994 (Keplerloesung)

maxP=1000;
T_vector1=linspace(T1,T2,maxP);
T_vector2=linspace(T2,T3,maxP);
for k=3:5
    Planets(k) =PlanetPQR(T_vector1, BaPa, BaPadot, k);
    Planets2(k)=PlanetPQR(T_vector2, BaPa, BaPadot, k);
end

T=linspace(T1,T2,1000);
%  Setze  Optionen  fuer  ode45 options  =  odeset('RelTol',1e-9);
options  =  odeset('RelTol',1e-9);
mJuF =1;  %Variation der Jupitermasse
%  Integriere  Differentialgleichungssystem
[t1 , YB1] =  ode45(@elliptic_orbit,  [T1  T2],  start1,  options, mJuF,...
    PosJT);
[~ , YB2]  =  ode45(@elliptic_orbit,  [T1  T2],  start2,  options, mJuF,...
    PosJT);
[~ , YB3]  =  ode45(@elliptic_orbit,  [T1  T2],  start3,  options, mJuF,...
    PosJT);
mJuF = 1.33333;
[t4 , YB4N]  =  ode45(@elliptic_orbit,  [T1  T2],  start1,  options, ...
    mJuF, PosJT);
for k=1:41 
  YB4(k,:) = YB4N(k,:);  %Nur bis zur Kollision
end
mJuF = 0.33333;
[~ , YB5]  =  ode45(@elliptic_orbit,  [T1  T2],  start1,  options, mJuF, PosJT);

% Endgeschwindigkeiten 17.07.1994
% SL 9
kollP =length(YB1)
v_final_SL9 = [YB1(kollP,4), YB1(kollP,5), YB1(kollP,6)]
% Jupiter
v_final_Jup = [Planets(NrP).xyz(1,maxP)-Planets(NrP).xyz(1,maxP-1),...
    Planets(NrP).xyz(2,maxP)-Planets(NrP).xyz(2,maxP-1),...
    Planets(NrP).xyz(3,maxP)-Planets(NrP).xyz(3,maxP-1)];
v_final_Jup = v_final_Jup/(Planets(NrP).Time(maxP)-...
              Planets(NrP).Time(maxP-1));
v_Jup = norm(v_final_Jup)*AE/secperday

% 2. Kosmische Geschwindigkeit für Jupiter zum Vergleich
v_Jup = sqrt(2*muJ/rJ)
v_Jup = v_Jup*AE/secperday;

%Relativgeschwindigkeit
vrel = v_final_SL9-v_final_Jup;
vrel = norm(vrel);
vrel = vrel*AE/secperday  %in km/s

% Energie
rho = 0.5 *1000   % Dichte in kg pro m³
D   = 1400        % Durchmesser Parent Body in m

E_kin = 0.5*(4/3)*rho*pi*(D/2)^3*vrel^2*1e6 % in kgm2/s² = Ws
TNT = 4.185*1e9;  % in Ws
E_equ = E_kin/TNT 

%_______________________________________________________________________
% Graphische Ausgabe 
% Stelle Ergebnis YB in 2D dar (plot);
header ='Bahnen SL 9 und Jupiter';
figure('Name',header);

% x-y-Koordinaten
subplot(1,2,1);
hold on;
% Kometenbahnen für verschiedene Parameter
p1=plot(YB1(:,1),YB1(:,2),'Color',Colors(2,:),'LineStyle','-', ...
    'LineWidth',2); 
p2=plot(YB2(:,1),YB2(:,2),'Color',Colors(8,:),'LineStyle',':', ...
    'LineWidth',1); 
p3=plot(YB3(:,1),YB3(:,2),'Color',Colors(8,:),'LineStyle','-.', ...
    'LineWidth',1); 
p4=plot(YB4(:,1),YB4(:,2),'Color',Colors(2,:),'LineStyle',':', ...
    'LineWidth',2); 
p5=plot(YB5(:,1),YB5(:,2),'Color',Colors(2,:),'LineStyle','-.', ...
    'LineWidth',2); 

p=plot(start1(1),start1(2),'o','Color',Colors(2,:)); %Position Komet @ T1
p.LineWidth=2;
LabelPoints(start1(1),start1(2),string(dt1,'dd.MM.yy')','E',0.25,0,...
    'FontSize',14,'Color',Colors(2,:));
p=plot(YB1(length(YB1),1),YB1(length(YB1),2),'d','Color',Colors(2,:)); 
p.LineWidth=3;
LabelPoints(YB1(length(YB1),1),YB1(length(YB1),2),string(dt2,...
    'dd.MM.yy')','E',0.25,0,'FontSize',14,'Color',Colors(2,:));
p=plot(YB2(length(YB2),1),YB2(length(YB2),2),'o','Color',Colors(8,:)); 
p.LineWidth=1;
p=plot(YB3(length(YB3),1),YB3(length(YB3),2),'o','Color',Colors(8,:)); 
p.LineWidth=1;
p=plot(YB4(length(YB4),1),YB4(length(YB4),2),'d','Color',Colors(2,:)); 
p.LineWidth=3;
dt4 = datetime(t4(length(YB4)),'ConvertFrom','juliandate');
LabelPoints(YB4(length(YB4),1),YB4(length(YB4),2),string(dt4,...
    'dd.MM.yy')','SW',0.35,0,'FontSize',14,'Color',Colors(2,:));
p=plot(YB5(length(YB5),1),YB5(length(YB5),2),'o','Color',Colors(2,:));
p.LineWidth=1;

p=plot(JupPos(1,:),JupPos(2,:)'); %Bahn Jupiter
p.Color=Colors(NrP,:);
p.LineWidth=2;
%Position Jupiter @ T1
p=plot(JupPos(1,1),JupPos(2,1),'o','Color',Colors(NrP,:)); 
LabelPoints(JupPos(1,1),JupPos(2,1),string(dt1,'dd.MM.yy')','E',0.25,0,...
    'FontSize',14,'Color',Colors(NrP,:));
p.LineWidth=1;
axis equal
grid on
xlabel('x in AE')
ylabel('y in AE');
xlim([-6.0 -3.0]);
ylim([-5 0]);
lg=legend([p1 p2 p3 p4 p5],kurve);
legend box off;
set(gca,'FontSize',18);

% x-z-Koordinaten
subplot(1,2,2);
hold on;
% Kometenbahnen
plot(YB1(:,1),YB1(:,3),'Color',Colors(2,:),'LineStyle','-','LineWidth',2); 
plot(YB2(:,1),YB2(:,3),'Color',Colors(8,:),'LineStyle',':','LineWidth',1); 
plot(YB3(:,1),YB3(:,3),'Color',Colors(8,:),'LineStyle','-.','LineWidth',1); 
plot(YB4(:,1),YB4(:,3),'Color',Colors(2,:),'LineStyle',':','LineWidth',2); 
plot(YB5(:,1),YB5(:,3),'Color',Colors(2,:),'LineStyle','-.','LineWidth',2); 

% Position Komet @ T1
p=plot(start1(1),start1(3),'o','Color',Colors(2,:));
p.LineWidth=2;
p=plot(YB1(length(YB1),1),YB1(length(YB1),3),'d','Color',Colors(2,:)); 
p.LineWidth=3;
p=plot(YB2(length(YB2),1),YB2(length(YB2),3),'o','Color',Colors(8,:));
p.LineWidth=1;
p=plot(YB3(length(YB3),1),YB3(length(YB3),3),'o','Color',Colors(8,:)); 
p.LineWidth=1;
p=plot(YB4(length(YB4),1),YB4(length(YB4),3),'d','Color',Colors(2,:)); 
p.LineWidth=3;
p=plot(YB5(length(YB5),1),YB5(length(YB5),3),'o','Color',Colors(2,:)); 
p.LineWidth=1;

% Position Komet @ T1
p=plot(start1(1),start1(3),'o','Color',Colors(2,:)); 
p.LineWidth=2;
LabelPoints(start1(1),start1(3),string(dt1,'dd.MM.yy')','E',0.25,0, ...
    'FontSize',14,'Color',Colors(2,:));
p=plot(YB1(length(YB1),1),YB1(length(YB1),3),'d','Color',Colors(2,:)); 
p.LineWidth=3;
LabelPoints(YB1(length(YB1),1),YB1(length(YB1),3),string(dt2, ...
    'dd.MM.yy')','E',0.25,0,'FontSize',14,'Color',Colors(2,:));
p=plot(YB2(length(YB2),1),YB2(length(YB2),3),'d','Color',Colors(8,:)); 
p.LineWidth=1;
p=plot(YB3(length(YB3),1),YB3(length(YB3),3),'d','Color',Colors(8,:)); 
p.LineWidth=1;
p=plot(YB4(length(YB4),1),YB4(length(YB4),3),'+','Color',Colors(2,:)); 
p.LineWidth=3;
dt4 = datetime(t4(length(YB4)),'ConvertFrom','juliandate');
LabelPoints(YB4(length(YB4),1),YB4(length(YB4),3),string(dt4, ...
    'dd.MM.yy')','NE',0.35,0,'FontSize',14,'Color',Colors(2,:));
p=plot(YB5(length(YB5),1),YB5(length(YB5),3),'d','Color',Colors(2,:)); 
p.LineWidth=1;

% Sonne
p=plot(0,0,'o');  
p.Marker='*';
p.Color=Colors(10,:);
p.LineWidth=3;
% Jupiter
p=plot(JupPos(1,:),JupPos(3,:)); 
p.Color=Colors(NrP,:);
p.LineWidth=2;
% Position Jupiter @ T1
p=plot(JupPos(1,1),JupPos(3,1),'o','Color',Colors(NrP,:)); 
LabelPoints(JupPos(1,1),JupPos(3,1),string(dt1,'dd.MM.yy')','N',0.1,0, ...
    'FontSize',14,'Color',Colors(NrP,:));
grid on
xlabel('x in AE')
ylabel('z in AE');
xlim([-6 -3]);
ylim([-0.25 0.25]);
set(gca,'FontSize',18);

%___________________________________________________________________

% 3-Dimensionale Darstellung

figure('Name',header);
p=plot3(YB1(:,1),YB1(:,2),YB1(:,3),'Color',Colors(2,:));
p.LineWidth=2;
hold on
p=plot3(YB1(1,1),YB1(1,2),YB1(1,3),'Color', Colors(2,:));
p.Marker='o';
p.LineWidth=2;
text(YB1(1,1),YB1(1,2),YB1(1,3), BaPaK.Name(NrK) + " " + string(dt1, ...
    'dd-MM-yy'),'HorizontalAlignment','center','FontSize',14,'Color', ...
    Colors(2,:));
SonneFP;
for iPlot = 3:5
    p=plot3(Planets(iPlot).xyz(1,1),Planets(iPlot).xyz(2,1), ...
        Planets(iPlot).xyz(3,1),'Color', Colors(iPlot,:));
    p.Marker='o';
    p.LineWidth=1;
    p=plot3(Planets(iPlot).xyz(1,maxP),Planets(iPlot).xyz(2,maxP), ...
        Planets(iPlot).xyz(3,maxP),'Color', Colors(iPlot,:));
    p.Marker='d';
    p.LineWidth=2;
    p=plot3(Planets2(iPlot).xyz(1,maxP),Planets2(iPlot).xyz(2,maxP), ...
        Planets2(iPlot).xyz(3,maxP),'Color', Colors(iPlot,:));
    p.Marker='+';
    p.LineWidth=1;
    p=plot3(Planets(iPlot).xyz(1,:),Planets(iPlot).xyz(2,:), ...
        Planets(iPlot).xyz(3,:),'Color', Colors(iPlot,:));
    p.LineWidth=2;
    p=plot3(Planets2(iPlot).xyz(1,:),Planets2(iPlot).xyz(2,:), ...
        Planets2(iPlot).xyz(3,:),'Color',Colors(iPlot,:),'lineStyle','-.');
    p.LineWidth=2;
    text(Planets(iPlot).xyz(1,1),Planets(iPlot).xyz(2,1), ...
        Planets(iPlot).xyz(3,1)," " + Planets(iPlot).Name + " " ...
        + string(dt1,'dd-MM-yy'),'FontSize',14,'Color',Colors(iPlot,:));
    text(Planets2(iPlot).xyz(1,1),Planets2(iPlot).xyz(2,1), ...
        Planets2(iPlot).xyz(3,1), " " + string(dt2,'dd-MM-yy'), ...
        'FontSize',14,'Color',Colors(iPlot,:));
    text(Planets2(iPlot).xyz(1,maxP),Planets2(iPlot).xyz(2,maxP), ...
        Planets2(iPlot).xyz(3,maxP), " " + string(dt3,'dd-MM-yy'), ...
        'HorizontalAlignment','right','FontSize',14,'Color', ...
        Colors(iPlot,:));
end

xlabel('x in AE')
ylabel('y in AE');
zlabel('z in AE');
xlim([-6 2]);
ylim([-6 2]);
xl = xlim;
yl = ylim;
[X,YB] = meshgrid(xl,yl);
surf(X,YB,zeros(size(X)))
shading flat
alpha 0.1

grid on
grid minor
set(gca,'FontSize',18);

%_____________________________________________________________________

% Funktionen

% Differentialgleichungssystem
function  dy  =  elliptic_orbit(t, y, mJu, PosJ)
dy  =  zeros(6,1);  %  Es  muss  ein  Spaltenvektor  zurückgegeben  werden 
x   =  zeros(3,1);
posJ = feval(PosJ,t);
muG      = 2.95479E-04;     % µG [AE^3/d^2]
muJ      = mJu*muG/1047.348644; % Jupitermasse G
x(1)   =  y(1);
x(2)   =  y(2);
x(3)   =  y(3);
dy(1)  =  y(4);
dy(2)  =  y(5);
dy(3)  =  y(6);
abstand1=norm(x);
abstand2=norm(x-posJ);
dy(4)  =  -muG*y(1)/ abstand1^3-muJ*(y(1)-posJ(1))/ abstand2^3;
dy(5)  =  -muG*y(2)/ abstand1^3-muJ*(y(2)-posJ(2))/ abstand2^3;
dy(6)  =  -muG*y(3)/ abstand1^3-muJ*(y(3)-posJ(3))/ abstand2^3; 
end


% Exakte Berechnung der Jupiterkoordinaten
function PosJ = JupiterKoord_ex(T, k, BaPa, BaPadot)
    % Bahnparameter Jupiter
    PlanetPos(k)=PlanetPQR(T, BaPa, BaPadot, k);
    PosJ=[PlanetPos(k).xyz(1,:);PlanetPos(k).xyz(2,:);...
        PlanetPos(k).xyz(3,:)];
end

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------
