% -------------------------------------------------------------------------
% FlugVenusParabelbahn.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Astrodynamik" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet die Parabel-Flugbahn einer Raumsonde zur Venus
% zu verschiedenen Startzeiten.
% Vorgegeben ist das Startdatum, gesucht ist das Ankunftsdatum mit dem
% schnellsten Transfer unabängig vom notwendigen Antriebsbedarf.
%
% Achtung alle Winkelberechnungen in ° !!!
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

% Initialisierung

% Bitte ein Datum auswählen
dt1 = datetime('2022-06-10 00:00:00');  % energetisch problematisch
dt1 = datetime('2022-09-10 00:00:00');  % energetisch problematisch
dt1 = datetime('2023-09-10 00:00:00');  % energetisch günstig
T1 = juliandate(dt1); % Julianisches Datum % Startdatum
MJuDa1 = juliandate(dt1,'modifiedjuliandate');% Modif. Julianisches Datum
DatumStr1 =  string(dt1,'dd.MM.yy');
T2 = T1+365; 

Aequi = 'Date';
% Einlesen der Bahnparameter
[BaPa,BaPadot]=OrbitParameter(T2,Aequi);
GS  = 2.95479E-04;     % [AE^3/Tage^2] für Sonne


%% Beginn Rechnung
%-------------------------------------------------------------------------%

maxT  = 365; % Berechnungsbereich
scFac = 1;   % Schrittweite 24 h 
T_vector     = linspace(T1,T1+maxT,maxT*scFac+1); %Start bis Scanbereich
calendardays = linspace(0,maxT,maxT*scFac+1);
ianf = 2;
iend = 3;

% Berechnung der Planetenbahnen nach Keplerloesung
for k=1:4
    Planets(k)=PlanetPQR(T_vector, BaPa, BaPadot, k);
end

% Bestimmung minimale TOF zur Venus (Parabelbahn)
% Wir vereinfachen uns das Problem, in  dem wir die Bewegung in der 
% Ekliptik-Ebene annehmen

rV1   = Planets(2).ekl(1,1);
phiV1 = rad2deg(Planets(2).ekl(2,1));
rE1   = Planets(3).ekl(1,1);
phiE1 = rad2deg(Planets(3).ekl(2,1));
% Bereich in dem nach TOF gesucht wird (Schrittweite 6h)
scanwidth = maxT*scFac;
for k=1:maxT*scFac+1
    rV2(k)   = Planets(2).ekl(1,k);
    phiV2(k) = wrapTo360(rad2deg(Planets(2).ekl(2,k)));
    phiE2(k) = wrapTo360(rad2deg(Planets(3).ekl(2,k)));
    if phiV2(k)-phiE1 < 180
        fac = 1;
    else
        fac = -1;
    end
    c(k)     = sqrt(rE1^2+rV2(k)^2-2*rE1*rV2(k)*cosd(phiV2(k)-phiE1));
    s(k)     = 0.5*(rE1+rV2(k)+c(k));
    tau(k)   = (sqrt(2)/3)*sqrt(s(k)^3/GS)*(1-fac*((s(k)-c(k))/s(k))^(3/2));
end
[DeltaTOF, kArr]   = min(abs(tau-calendardays));
kArr = kArr-1;
DayofArrival = kArr/scFac;
TArr = T1 + DayofArrival;
dtArr = datetime(TArr,'ConvertFrom','JulianDate');
DatumStr2 =  string(dtArr,'dd.MM.yy');  %Ankunftsdatum auf Venusbahn


%% Berechnung Transferbahn (Halbparameter p und  Richtung Scheitel u0)
%
% Numerische Bestimmung p, geht auch analytisch, aber warum, wenn man
% MATLAB hat
u = linspace(0,180,10000);
p1 = rE1*(1+cosd(phiE1-u));
p2 = rV2(kArr)*(1+cosd(phiV2(kArr)-u));
gamma = phiV2(kArr)-phiE1;
delp =p1-p2;
kz = 1;
for k=2:length(u)
    if delp(k)*delp(k-1)  <0
        kend(kz) = k;
        ptemp(kz)= p1(k);
        kz =  kz+1;
    end
end
p_par = ptemp(1);
% Analytische Berechnung (siehe Walter, Astronautics, 2018)
p_par = (2*rE1*rV2(kArr)/c(kArr)^2)*(sind(gamma/2))^2;
p_par = p_par*(rE1+rV2(kArr)+2*sqrt(rE1*+rV2(kArr))*cosd(gamma/2));
% Bestimmung u0
u0=phiE1+acosd(p_par/rE1-1);

% Berechnung Transferbahnen 
u = linspace(phiE1,180,1000);
x_par1 = p_par*cosd(u)./(1+cosd(u-u0));
y_par1 = p_par*sind(u)./(1+cosd(u-u0));
u = linspace(phiE1,phiV2(kArr),1000);
x_par = p_par*cosd(u)./(1+cosd(u-u0));
y_par = p_par*sind(u)./(1+cosd(u-u0));


%% Graphische Ausgabe, Print-Out
%18-------------------------------------------------------------------------%
ttlprint = "Lambert-Problem Parabeltransfer zur Venus";
fprintf('\n %s', ttlprint);
fprintf('\n\n')
fprintf('|  rE1 (AE)  |  rV2 (AE)  |  Abflug    | Ankunft    |   p (AE)  |\n');
fprintf('|  %7.4f   |  %7.4f   |  %s  | %s   |  %7.4f  |\n',...
                    rE1, rV2(kArr), DatumStr1, DatumStr2, p_par);
fprintf('\n')
fprintf('|  phiE1(°)  |  phiE2(°) |  phiV1 (°)  |  phiV2 (°) | gamma (°) |\n');
fprintf('|  %7.3f   |  %7.3f  |   %7.3f   |  %7.3f   |  %+7.3f  |\n',...
                    phiE1, phiE2(kArr), phiV1, phiV2(kArr), gamma);
fprintf('\n')

% Trajektorien
header1='Innere Planeten: Projektion Bahnen auf die Ekliptik';
figure('Name',header1);
%Planetenbahnen
for iPlot = ianf:iend
    plot(Planets(iPlot).xyz(1,:),Planets(iPlot).xyz(2,:),'Color', ...
         Colors(iPlot,:));
    hold on
    axis equal
end
% Sonne 
PlotCircle(0,0,0.07,Colors(10,:),2);  
%Positionen Start und Ankunft
k = 1;
for iPlot =  ianf:iend
    p(k)=plot(Planets(iPlot).xyz(1,1),Planets(iPlot).xyz(2,1),'+-',...
              'Color', Colors(iPlot,:),'LineWidth',2,'MarkerSize',8,...
              'MarkerFaceColor',Colors(iPlot,:));
    k=k+1;
    plot([0 Planets(iPlot).xyz(1,1)],[0 Planets(iPlot).xyz(2,1)],...
         ':','Color', Colors(iPlot,:),'LineWidth',1);
    p(k)=plot(Planets(iPlot).xyz(1,kArr),Planets(iPlot).xyz(2,kArr),'d-',...
              'Color', Colors(iPlot,:),'LineWidth',2,'MarkerSize',8,...
              'MarkerFaceColor',Colors(iPlot,:));
    k=k+1;
    plot([0 Planets(iPlot).xyz(1,kArr)],...
         [0 Planets(iPlot).xyz(2,kArr)],...
          ':','Color', Colors(iPlot,:),'LineWidth',1);
end
k=k+1;
% Plot der Transferbahn !
plot(x_par1,y_par1, 'Color', Colors(1,:),'LineWidth',1, ...
    'LineStyle',Style(4));%Transferparabel
p(5)=plot(x_par,y_par, 'Color', Colors(1,:),'LineWidth',2);%Transferparabel
ylim([-1.1 1.1])
xlim([-1.1 1.1]);
grid on;
grid minor,
header2 = strjoin([DatumStr1,' - ', DatumStr2]);
ttl=title(header2);
set(ttl,'FontSize',14, 'FontWeight','normal');
xlabel('x in AE')
ylabel('y in AE');
legend(p, strjoin([Planets(ianf).Name, "(bei Start)"]),...
          strjoin([Planets(ianf).Name, "(bei Ankunft)"]),...
          strjoin([Planets(iend).Name, "(bei Start)"]),...
          strjoin([Planets(iend).Name, "(bei Ankunft)"]),...
          "Transferbahn", 'location','bestoutside','numcolumns',1);
legend boxoff;
set(gca,'FontSize',16);

 
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------
%  
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------
