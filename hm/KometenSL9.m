% -------------------------------------------------------------------------
% KometenSL9.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Komet Shoemaker-Levy 9
% Berechnet die Kepler-Bahn des Kometen D/1993 F2-D (Shoemaker-Levy 9) als 
% Funktion der Zeit und vergleicht diese Rechnungen mit den Werten der 
% Berechnung des JPL. 
% (Aequinoktium und Ekliptik des Datum J2000)
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

% Initialisierung
T0 = juliandate(datetime(2000,1,1,12,0,0));
Aequi = 'J2000';
% Einlesen der Bahnparameter der Planeten und deren Ableitung 
[BaPa,BaPadot]=OrbitParameter(T0, Aequi);
% Einlesen der Bahnparameter der Kometen 
BaPaK=KometParameter('Kometen.csv');

GM   = 2.95479E-04;     % µG [AE^3/d^2]
pi2  = 2*pi;

% Berechnung der Planetenbahnen nach Keplerloesung
u=linspace(0,360,361);
% Iteration durch die betrachteten Planetemn 
for iPlot=3:6
    ePn=BaPa.eP(iPlot);
    aPn=BaPa.aP(iPlot);
    Planets(iPlot).Name=string(BaPa.Name(iPlot));
    fac=aPn*(1-ePn*ePn)./(1+ePn*cosd(u));
    x = cosd(u).*fac;
    y = sind(u).*fac;
    z = cosd(u).*0;
    yin=[x;y;z];
    iPn=deg2rad(BaPa.iP(iPlot));
    OmegaPn=deg2rad(BaPa.OmegaP(iPlot));
    PiqPn=deg2rad(BaPa.PiqP(iPlot));
    wPn=PiqPn-OmegaPn;
    yout=mtimes(PQR(-OmegaPn,-iPn,-wPn),yin);
    Planets(iPlot).xyz(1,:)=yout(1,:);
    Planets(iPlot).xyz(2,:)=yout(2,:);
    Planets(iPlot).xyz(3,:)=yout(3,:);
end

dt1 = datetime('1994-01-01 00:00:00');
T1 = juliandate(dt1); % Julianisches Datum
MJuDa1 = juliandate(dt1,'modifiedjuliandate');% Modif. Julianisches Datum
DatumStr1 =  string(dt1,'dd.MM.yyyy');

dt2 = datetime('1994-07-17 00:00:00');
T2 = juliandate(dt2); % Julianisches Datum  %Bedeckungsdatum
MJuDa2 = juliandate(dt2,'modifiedjuliandate');% Modif. Julianisches Datum
DatumStr2 =  string(dt2,'dd.MM.yyyy');

Year1=year(dt1);
Year2=year(dt2);
Nges=round(T2-T1)+1;    % Anzahl Tage 
T=linspace(T1,T2,Nges)+1;

% Berechnung der Schrittweiten 
StepsperDay=1;
NrStepsMarker=10*StepsperDay;
NrStepsLabel=NrStepsMarker*2;
MaxStepsLabel=floor(Nges/NrStepsLabel);

%Berechnung der Planetenposition
for k=3:6
    PlanetPos(k)=PlanetPQR(T, BaPa, BaPadot, k);
end

%%
%--------------------------------------------------------------------------
% Beginn Rechnung

% Auswahl des Kometen
HK=5;
name(iPlot,:)=BaPaK.Name(HK);

% Berechnung der heliozentrischen Koordinaten des Jupiter und des Kometen
% nach Keplerloesung von T1 bis T2

Komets   = KometPQR(GM, T, BaPaK, HK);
Jup.xyz=PlanetPos(5).xyz;

% Definieren der Optionen zum Einlesen 
opts = delimitedTextImportOptions("NumVariables", 11);
% Bereich und Trennzeichen 
opts.DataLines = [2, Inf];
opts.Delimiter = ";";
% Spaltennamen und - typen 
opts.VariableNames = ["Datum", "JD", "l_SL", "b_SL", "d_SL", "l_Jup",...
    "b_Jup", "d_Jup", "x", "y", "z"];
opts.VariableTypes = ["datetime", "double", "double", "double", ...
    "double", "double", "double", "double", "double", "double", "double"];
opts = setvaropts(opts, 1, "InputFormat", "yyyy-MM-dd");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Einlesen der Datei 
SL1 = readtable("SL9.csv", opts);
clear opts

% Koordinatentransformation 
Jup.ekl=CalcAnglesfromXYZ(Jup.xyz);
SL9.ekl=CalcAnglesfromXYZ(Komets.xyz);
SL9.xyz=Komets.xyz;


%________________________________________________________________________
%
% Ausgabe im Command Window 
for n=1:Nges
   DateStr=string(datetime(T(n),'ConvertFrom','juliandate'),'dd.MM. HH:mm');
   fprintf('\n Datum : %s | %10.1f ', DateStr, floor(T(n)));
   fprintf('|SL9:  {alpha} %5.2f°  {delta} %4.2f° {r} %6.4f AE   ',...
     wrapTo360(rad2deg(SL9.ekl(2,n))), rad2deg(SL9.ekl(3,n)), SL9.ekl(1,n));
   fprintf('|Jup:  {alpha} %5.2f°  {delta} %4.2f° {r} %6.4f AE',...
     wrapTo360(rad2deg(Jup.ekl(2,n))), rad2deg(Jup.ekl(3,n)), Jup.ekl(1,n));
end
fprintf('\n');

%__________________________________________________________________________
%%
% Graphische Ausgabe

%__________________________________________________________________________
% Figure 1

header1='Komet SL9';
figure('Name',header1);

for iPlot=1:1
    p(iPlot)=plot(Komets(iPlot).xyz(1,:),Komets(iPlot).xyz(2,:), ...
        'Color', Colors(iPlot,:));
    hold on;
    plot(SL1.x,SL1.y,'Color', Colors(iPlot+1,:));
    p(iPlot).LineWidth=2;
    hold on;
    plot(Jup.xyz(1,:),Jup.xyz(2,:),':+','MarkerIndices',1:NrStepsMarker:...
        length(Jup.xyz(1,:)),'LineWidth',1,'Color',Colors(5,:));
    plot(Komets(iPlot).xyz(1,:),Komets(iPlot).xyz(2,:),':+', ...
        'MarkerIndices',1:NrStepsMarker:length ...
        (Komets(iPlot).xyz(1,:)),'LineWidth',1,'Color',Colors(iPlot,:));
    plot(SL1.x,SL1.y,':+','MarkerIndices',1:NrStepsMarker:length ...
        (Komets(iPlot).xyz(1,:)),'LineWidth',1,'Color',Colors(iPlot+1,:));
    h=LabelPoints(Jup.xyz(1,1), Jup.xyz(2,1), ...
        'Jupiter -1.1.94-','W',0.04,'FontSize',12,'Color',Colors(5,:));
    h=LabelPoints(Jup.xyz(1,Nges), Jup.xyz(2,Nges), ...
        'Jupiter -17.7.94- ','W',0.04,'FontSize',12,'Color',Colors(5,:));
    h=LabelPoints(SL1.x(1), SL1.y(1),'- SL9 1.1.94','E', ...
        0.04,'FontSize',12,'Color',Colors(iPlot+1,:));
    h=LabelPoints(SL1.x(Nges), SL1.y(Nges),'- SL9 17.7.94 ','E', ...
        0.04,'FontSize',12,'Color',Colors(iPlot+1,:));
    h=LabelPoints(SL9.xyz(1,1), SL9.xyz(2,1),'- SL9 (ber.)','E', ...
        0.04,'FontSize',12,'Color',Colors(iPlot,:));
    h=LabelPoints(SL9.xyz(1,Nges), SL9.xyz(2,Nges),'- SL9 (ber.)','E', ...
        0.04,'FontSize',12,'Color',Colors(iPlot,:));
    plot(SL9.xyz(1,Nges), SL9.xyz(2,Nges),'+','Color', Colors(iPlot,:), ...
        'MarkerSize',5,'MarkerFaceColor',Colors(iPlot,:));
    plot(Jup.xyz(1,1), Jup.xyz(2,1),'o','Color', Colors(5,:), ...
        'MarkerSize',5,'MarkerFaceColor',Colors(5,:));
    plot(Jup.xyz(1,Nges), Jup.xyz(2,Nges),'o','Color', Colors(5,:), ...
        'MarkerSize',5,'MarkerFaceColor',Colors(5,:));

    % Sonne, Ekliptik und Planeten
    SonneFP
    
    a=[0,Jup.xyz(1,1)];
    b=[0,Jup.xyz(2,1)];
    plot(a,b,'--','Color',Colors(5,:));
    a=[0,Jup.xyz(1,Nges)];
    b=[0,Jup.xyz(2,Nges)];
    plot(a,b,'--','Color',Colors(5,:));
end

% Einstellungen der Achsen 
xlim([-5 -3.5]);
ylim([-4 -2.5]);
axis square;

% Iteration durch die Planeten 
for iPlot = 3:5
    plot(Planets(iPlot).xyz(1,:),Planets(iPlot).xyz(2,:),'Color', ...
        Colors(iPlot,:));
 end
for iPlot = 3:4
    plot(PlanetPos(iPlot).xyz(1,1),PlanetPos(iPlot).xyz(2,1),'o', ...
        'Color',Colors(iPlot,:),'MarkerSize',5,'MarkerFaceColor', ...
        Colors(iPlot,:));
    plot(PlanetPos(iPlot).xyz(1,Nges),PlanetPos(iPlot).xyz(2,Nges),'o', ...
        'Color', Colors(iPlot,:),'MarkerSize',5,'MarkerFaceColor', ...
        Colors(iPlot,:));
    h=LabelPoints(PlanetPos(iPlot).xyz(1,1), PlanetPos(iPlot).xyz(2,1), ...
        Planets(iPlot).Name,'S',0.03,'FontSize',12,'Color',Colors(iPlot,:));
end

grid on;
grid minor;

% Beschriftungen 
header2 = strcat(BaPaK.Name(HK),' am Jupiter ');
ttl=title(header2, 'FontSize',12);
ttl=xlabel('x in AE');
ttl=ylabel('z in AE');
set(gca,'FontSize',14);

%__________________________________________________________________________
% Figure 2 

header1='Komet SL9';
figure('Name',header1);

for iPlot=1:1
    p(iPlot)=plot(Komets(iPlot).xyz(1,:),Komets(iPlot).xyz(3,:), ...
        'Color',Colors(iPlot,:));
    hold on;
    plot(SL1.x,SL1.z,'Color', Colors(iPlot+1,:));
    p(iPlot).LineWidth=2;
    plot(Jup.xyz(1,:),Jup.xyz(3,:),':+','MarkerIndices', ...
        1:NrStepsMarker:length(Jup.xyz(1,:)),'LineWidth',1, ...
        'Color',Colors(5,:));
    plot(Komets(iPlot).xyz(1,:),Komets(iPlot).xyz(3,:),':+', ...
        'MarkerIndices',1:NrStepsMarker:length(Komets(iPlot).xyz(1,:)), ...
        'LineWidth',1,'Color',Colors(iPlot,:));
    plot(SL1.x,SL1.z,':+', ...
        'MarkerIndices',1:NrStepsMarker:length(SL1.x), ...
        'LineWidth',1,'Color',Colors(iPlot+1,:));
    h=LabelPoints(Jup.xyz(1,1), Jup.xyz(3,1), ...
        ' - Jupiter 1.1.94','NE',0.04,'FontSize',10,'Color',Colors(5,:));
    h=LabelPoints(Jup.xyz(1,Nges), Jup.xyz(3,Nges), ...
        'Jupiter 17.7.94 - ','W',0.04,'FontSize',10,'Color',Colors(5,:));
    h=LabelPoints(SL1.x(1), SL1.z(1), ...
        ' SL9 1.1.94 -','W',0.04,'FontSize',10,'Color',Colors(iPlot+1,:));
    h=LabelPoints(SL1.x(Nges), SL1.z(Nges), ...
        '- SL9 17.7.94 ','E',0.04,'FontSize',10,'Color',Colors(iPlot+1,:));
    h=LabelPoints(SL9.xyz(1,1), SL9.xyz(3,1), ...
        '- SL9 (ber.) 1.1.94','E',0.04,'FontSize',10, ...
        'Color',Colors(iPlot,:));
    h=LabelPoints(SL9.xyz(1,Nges), SL9.xyz(3,Nges),...
        '- SL9 (ber.) 17.7.94 ','E',0.04,'FontSize',10, ...
        'Color',Colors(iPlot,:));
    plot(SL9.xyz(1,Nges), SL9.xyz(3,Nges),'+','Color',Colors(iPlot,:), ...
        'MarkerSize',5,'MarkerFaceColor',Colors(iPlot,:));
    plot(Jup.xyz(1,1), Jup.xyz(3,1),'o','Color', Colors(5,:), ...
        'MarkerSize',5,'MarkerFaceColor',Colors(5,:));
    plot(Jup.xyz(1,Nges), Jup.xyz(3,Nges),'o','Color', Colors(5,:), ...
        'MarkerSize',5,'MarkerFaceColor',Colors(5,:));

    % Sonne, Ekliptik und Planeten
    SonneFP;
    
    a=[0,Jup.xyz(1,1)];
    b=[0,Jup.xyz(3,1)];
    plot(a,b,'--','Color',Colors(5,:));
    a=[0,Jup.xyz(1,Nges)];
    b=[0,Jup.xyz(3,Nges)];
    plot(a,b,'--','Color',Colors(5,:));
end

% Einstellungen der Achsen 
xlim([-6 2]);
ylim([-0.15 0.15]);
axis square;

% Iteration durch die Planeten 
for iPlot = 3:5
    plot(Planets(iPlot).xyz(1,:),Planets(iPlot).xyz(3,:), ...
        'Color', Colors(iPlot,:));
 end
for iPlot = 3:4
    plot(PlanetPos(iPlot).xyz(1,1),PlanetPos(iPlot).xyz(3,1),'o', ...
        'Color',Colors(iPlot,:),'MarkerSize',5,'MarkerFaceColor', ...
        Colors(iPlot,:));
    h=LabelPoints(PlanetPos(iPlot).xyz(1,1), PlanetPos(iPlot).xyz(3,1), ...
        Planets(iPlot).Name,'S',0.02,'FontSize',8,'Color',Colors(iPlot,:));
end

grid on;
grid minor;

% Beschriftungen 
header2 = strcat(BaPaK.Name(HK),' am Jupiter ');
ttl=title(header2, 'FontSize',12);
ttl=xlabel('x in AE');
ttl=ylabel('z in AE');
set(gca,'FontSize',14);

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------
