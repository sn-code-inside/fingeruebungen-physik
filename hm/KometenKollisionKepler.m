%%-------------------------------------------------------------------------
% KometenKollison02.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Astrodynamik" aus
% "FingerÃ¼bungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% 
% Simulation der Bewegung des Kometen Swift-Tuttle 
% durch Berechnung der Bahn des Kometen als Funktion 
% der Zeit über die Lösung der Kepler-Gleichung (Rechnung in AE)
% unter fallweiser Benutzung einer parabolischen Naeherung.
% (Aequinoktium und Ekliptik Datum bzw. J2000)
% Die Bahn für die Erde wird entweder durch eine einfache Lösung der KG, 
% oder durch eine hochgenaue Reihenentwicklung bestimmt.
%
% Berechnung für Perihelumgebung 1992 und 2126.
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

% Initialisierung
% Einlesen der Bahnparameter der Kometen und Asteroiden 
BaPaK=KometParameter('Data/Kometen.csv');
GM      = 2.95479E-04;     % µG [AE^3/d^2]
AEk   = 1.495978707e08;         % AE in km
PNamen = ["Merkur" "Venus" "Erde" "Mars" "Jupiter" "Saturn" "Uranus" ...
          "Neptun" ... %"9P/Tempel 1" 
          "Swift-Tuttle" "Sonne"];

% Initialisierung & Input

% Input
% Auswahl des Komet (Nr 14 für Swift-Tuttle eingeben)

% Zeitraum und Darstellung Perihelnaehe Swift-Tuttle
% fuer 1992
HK=14;
dt1 = datetime('1992-11-01 00:00:00');
dt2 = datetime('1993-02-28 00:00:00');
Aequi = 'date';

% % Zeitraum und Darstellung Perihelnaehe Swift-Tuttle
% % fuer 2126
HK=14;
dt1 = datetime('2126-04-01 00:00:00');
dt2 = datetime('2126-07-31 00:00:00');
Aequi = 'Date';

T1 = juliandate(dt1);
T2 = juliandate(dt2);

% Einlesen der Bahnparameter der Planeten und deren Ableitungen 
[BaPa,BaPadot]=OrbitParameter(T1, Aequi); 

% Berechnung fuer kurze Zeitraeume 
long=0;
Nges=10000; %Zahl der Stuetzstellen im Intervall T1 - T2
T_vector  = linspace(T1,T2,Nges);      %Berechnungsintervall
T_vector2 = linspace(T1,T2+1000,Nges); %Berechnungsintervall länger
StepsPerDay= Nges/(T2-T1); %Zahl der Stuetzstellen pro Tag
DatumStr1 =  string(dt1,'dd.MM.yy');
DatumStr2 =  string(dt2,'dd.MM.yy');
YearStr1 =  string(dt1,'MMM-yyyy');
YearStr2 =  string(dt2,'MMM-yyyy');

%Input
StepsMarker= 5; %Zahl der markierten Tage 
%(Richtwert 5 bei T2-T1=1 Jahr usw.)
NrStepsMarker=round(StepsMarker*StepsPerDay); %Zahl der markierten Tage
%Input
StepsLabel=5; %Zahl der angezeigten Tage
%((Richtwert 2 T2-T1=1 Jahr usw.)
NrStepsLabel=StepsLabel*NrStepsMarker; %Zahl der angezeigten Tage
MaxStepsLabel=floor(Nges/NrStepsLabel);
XaxAE=1.2;

%%-------------------------------------------------------------------------
%
%Beginn Rechnung

% Berechnung der heliozentrischen Koordinaten der Erde und 
% der anderen Planeten von T1 bis T2 
for k=1:4
    PlanetKep(k)  = PlanetPQR(T_vector, BaPa, BaPadot,k);
    PlanetKep2(k) = PlanetPQR(T_vector2, BaPa, BaPadot,k);
end

% Überschreiben der Bahndaten Erde mit genaueren Koordinaten
SunData  = PerturbImport('SonnePos.csv');
epsE = EpsErde(T_vector(1));               % Ekliptikschiefe @ T1 
t1Jh = Jd2JJht(T_vector);   %Umrechnung Julianische Jahrhunderte
SunPos   = SonneExakt(t1Jh,SunData,epsE,'AE');
PlanetKep(3).xyz  = -SunPos.xyz;
PlanetKep(3).ekl  = CalcAnglesfromXYZ(PlanetKep(3).xyz);

% Berechnung der heliozentrischen Koordinaten 
% des Kometen nach Loesung mit Stumppfschen Funktionen von T1 bis T2 
Komet  = KometPQR(GM, T_vector, BaPaK,HK);

% Berechnung der Planetenbahnen Ellipsengleichung von T1 bis T2 
Planets=PlanetOrbits(BaPa);

%Berechnung Abstand Erde-Komet
rek=Komet.xyz-PlanetKep(3).xyz;
rekbetr=sqrt(rek(1,:).*rek(1,:)+rek(2,:).*rek(2,:)+rek(3,:).*rek(3,:));
[minrekbetr,ErdNaehe] = min(rekbetr);

%Berechnung Abstand Sonne-Komet und Perihel
rsk=Komet.xyz ;
rskbetr=sqrt(rsk(1,:).*rsk(1,:)+rsk(2,:).*rsk(2,:)+rsk(3,:).*rsk(3,:));
[minrskbetr,Perihel] = min(rskbetr);

%%-Graphische Ausgabe ---------------------------------------------------
%-----------------------------------------------------------------------% 

header2 = strjoin([strcat(YearStr1',' bis '),YearStr2]);
header2 = strjoin([BaPaK.Name(HK),': ',header2]);
header1 = 'Komet Swift-Tuttle';
NSteps  = length(PlanetKep(3).xyz(1,:));


% -------------------------------------------------------------------------
% Figure 1: 2-D-Ausgabe Kometenbahnen  
% -------------------------------------------------------------------------

figure('Name',header1);
for ifig = 1:2   
% x-y-Ebene
subplot(1,2,ifig);         % linkes Schaubild in Figure 1, ifig=1
% Plot Komet
plot(Komet.xyz(1,:),Komet.xyz(1+ifig,:),'+','MarkerIndices', ...
    1:NrStepsMarker:NSteps,'LineWidth',2,'Color', ...
    Colors(9,:));
hold on;
plot(PlanetKep(3).xyz(1,:),PlanetKep(3).xyz(1+ifig,:),'+','MarkerIndices', ...
    1:NrStepsMarker:NSteps,'LineWidth',2,'Color', ...
    Colors(3,:));
% Einstellungen der Achsen 
xlim([-XaxAE XaxAE]);
ylim([-XaxAE XaxAE]);
grid on;
grid minor;
% Beschriftungen der Zeitpunkte 
for k=1:MaxStepsLabel
    dt4 = datetime(Komet.Time((k-1)*NrStepsLabel+1),'ConvertFrom', ...
        'juliandate');
    mylabels(k,:)= datestr(dt4,'dd.mm.');
end
for k=1:MaxStepsLabel
    Xlab(k)=Komet.xyz(1,(k-1)*NrStepsLabel+1);
    YZlab(k)=Komet.xyz(ifig+1,(k-1)*NrStepsLabel+1);
    XlabE(k)=PlanetKep(3).xyz(1,(k-1)*NrStepsLabel+1);
    YZlabE(k)=PlanetKep(3).xyz(ifig+1,(k-1)*NrStepsLabel+1);
end
LabelPoints(Xlab,YZlab,mylabels,'W',0.01,0,'FontSize',8,'Color', ...
    Colors(9,:));
LabelPoints(XlabE,YZlabE,mylabels,'N',0.01,0,'FontSize',8,'Color', ...
    Colors(3,:));
 
% Sonne , Richtung Fruehlingspunkt
SonneFP;
% Darstellung Merkur, Venus- und Erdbahn
for iPlot = 1:3
    plot(Planets(iPlot).xyz(1,:),Planets(iPlot).xyz(ifig+1,:),'Color', ...
        Colors(iPlot,:),'LineWidth',2);
end
% Darstellung Position ErdNaehe
    line ([PlanetKep(3).xyz(1,ErdNaehe),Komet.xyz(1,ErdNaehe)], ...
          [PlanetKep(3).xyz(ifig+1,ErdNaehe),Komet.xyz(ifig+1,ErdNaehe)],...
           'Color',Colors(4,:),'lineWidth',2);
    plot(Komet.xyz(1,ErdNaehe),Komet.xyz(ifig+1,ErdNaehe),'o', ...
        'Color',Colors(4,:),'MarkerSize',5,'MarkerFaceColor',Colors(4,:));
    h=LabelPoints(Komet.xyz(1,ErdNaehe), ...
        Komet.xyz(ifig+1,ErdNaehe),'Erdnähe','W',0.01,'FontSize',8, ...
        'Color',Colors(3,:));
% Darstellung Position Perihel
    line ([0,Komet.xyz(1,Perihel)],[0,Komet.xyz(ifig+1,Perihel)],...
           'Color',Colors(10,:),'lineWidth',2);
    plot(Komet.xyz(1,Perihel),Komet.xyz(ifig+1,Perihel),'o', ...
        'Color',Colors(10,:),'MarkerSize',5,'MarkerFaceColor',Colors(10,:));
    h=LabelPoints(Komet.xyz(1,Perihel), ...
        Komet.xyz(ifig+1,Perihel),'Perihel','W',0.01,'FontSize',8, ...
        'Color',Colors(10,:));

% Beschriftungen 
ht=title(header2);
set(ht,'FontSize',12, 'FontWeight', 'normal');
axis square;
xlabel('x in AE');
if ifig ==1 
    ylabel('y in AE');
else
    ylabel('z in AE');
end
set(gca,'FontSize',14);
end


% -------------------------------------------------------------------------
% Figure 2: 3D-Darstellung im Zeitbereich um Perihel 1992 o. 2126 
% -------------------------------------------------------------------------

dtErdNaehe = datestr(datetime(Komet.Time(ErdNaehe),'convertfrom','juliandate'),...
               ' dd.mm.yy ');
dtPerihel  = datestr(datetime(Komet.Time(Perihel),'convertfrom','juliandate'),...
               ' dd.mm.yy');
PNamen2    = PNamen(1:4);
PNamen2(5) = PNamen(9);
PNamen2(6) = PNamen(10);
PNamen2(7) = 'Erdnähe';
PNamen2(8) = 'Perihel';


figure('Name',header1);
xAx = 2;

hp2(5)=plot3(Komet.xyz(1,:),Komet.xyz(2,:), ...
        Komet.xyz(3,:),'Color', ...
        Colors(9,:),'LineWidth',2);
hold on
axis equal
% Achse Sonne-Fruehlingspunkt 
SonneFP;
TimeStep = 1800;
msize    = 10;

for iPlot= 1:4
  hp2(iPlot) = plot3(PlanetKep(iPlot).xyz(1,:), ...
      PlanetKep(iPlot).xyz(2,:), ...
      PlanetKep(iPlot).xyz(3,:),'Color', Colors(iPlot,:), ...
      'linewidth',2);
  plot3(PlanetKep2(iPlot).xyz(1,:),PlanetKep2(iPlot).xyz(2,:), ...
      PlanetKep2(iPlot).xyz(3,:),'Color', Colors(iPlot,:), ...
      'linewidth',1);
end
hp2(6)= plot(0, 0,'o','color',Colors(10,:), ...
'markersize',msize+2,'markerfacecolor',Colors(10,:),'linewidth',2);

labelrange = 1:TimeStep:NSteps;
XLabel   = PlanetKep(3).xyz(1,labelrange);
YLabel   = PlanetKep(3).xyz(2,labelrange);
ZLabel   = PlanetKep(3).xyz(3,labelrange);
XLabelK  = Komet.xyz(1,labelrange);
YLabelK  = Komet.xyz(2,labelrange);
ZLabelK  = Komet.xyz(3,labelrange);

dtKstr = datestr(datetime(PlanetKep(3).Time(labelrange),'ConvertFrom',...
                'juliandate'),' dd.mm.');

%Zeitdaten
plot3(XLabel, YLabel, ZLabel,'+', ...
     'color',Colors(3,:),'markerfacecolor',Colors(9,:),...
     'markersize',msize,'linewidth',2);
text(XLabel, YLabel, ZLabel, dtKstr, 'FontSize',9,...
            'Color', Colors(3,:));

plot3(XLabelK, YLabelK, ZLabelK, '+', ...
     'color',Colors(9,:),'markerfacecolor',Colors(9,:),...
     'markersize',msize,'linewidth',2);
text(XLabelK,YLabelK, ZLabelK, dtKstr,'FontSize',9,...
            'Color', Colors(9,:));

hp2(7)=line([PlanetKep(3).xyz(1,ErdNaehe) Komet.xyz(1,ErdNaehe)],...
     [PlanetKep(3).xyz(2,ErdNaehe) Komet.xyz(2,ErdNaehe)],...
     [PlanetKep(3).xyz(3,ErdNaehe) Komet.xyz(3,ErdNaehe)],...
     'Color', Colors(4,:),'Linewidth',3);
plot3(PlanetKep(3).xyz(1,ErdNaehe), ...
    PlanetKep(3).xyz(2,ErdNaehe),PlanetKep(3).xyz(3,ErdNaehe), ...
    'o','Color', Colors(iPlot,:),'MarkerSize',5, ...
    'MarkerFaceColor',Colors(4,:));
plot3(Komet.xyz(1,ErdNaehe), ...
    Komet.xyz(2,ErdNaehe), Komet.xyz(3,ErdNaehe), ...
    'o','Color', Colors(iPlot,:),'MarkerSize',5, ...
    'MarkerFaceColor',Colors(4,:));
text(Komet.xyz(1,ErdNaehe), ...
     Komet.xyz(2,ErdNaehe), Komet.xyz(3,ErdNaehe),...
     dtErdNaehe, 'FontSize',9,'Color', Colors(4,:));

hp2(8)=line([0 Komet.xyz(1,Perihel)],...
     [0 Komet.xyz(2,Perihel)],...
     [0 Komet.xyz(3,Perihel)],...
     'Color', Colors(10,:),'Linewidth',3);
plot3(Komet.xyz(1,Perihel), ...
    Komet.xyz(2,Perihel), Komet.xyz(3,Perihel), ...
    'o','Color', Colors(10,:),'MarkerSize',5, ...
    'MarkerFaceColor',Colors(10,:));
text(Komet.xyz(1,Perihel), ...
     Komet.xyz(2,Perihel), Komet.xyz(3,Perihel),...
     dtPerihel, 'FontSize',9,'Color', Colors(10,:));

% Einstellungen der Achsen 
xlim([-xAx xAx]);
ylim([-xAx xAx]);
xl = xlim;
yl = ylim;
[X,Y] = meshgrid(xl,yl);
surf(X,Y,zeros(size(X)));
shading flat
alpha 0.1
grid on;
axis equal;
legend(hp2, PNamen2, 'location', 'westoutside')
legend box off
axL = 1.6;
xlim([-axL axL]);
ylim([-axL axL]);
zlim([-axL axL]);
xlabel('x in AE')
ylabel('y in AE');
zlabel('z in AE');
ht = title(strjoin(["Zeitraum : ",DatumStr1,"-",DatumStr2]));
set(ht,'FontSize',14,'FontWeight','normal');
set(gca,'FontSize',14);


%% ----------------------------------------------------------------------
%  Ausgabe Werte für Erdnähe 

dtErdNaehe = datestr(datetime(Komet.Time(ErdNaehe),'convertfrom','juliandate'),...
               ' dd.mm.yyyy');
dtPerihel  = datestr(datetime(Komet.Time(Perihel),'convertfrom','juliandate'),...
               ' dd.mm.yyyy');

fprintf('\n')
fprintf('\n')
fprintf('\n Berechnung für 2126 über Lösung der Kepler-Gleichung')
fprintf('\n')
fprintf('\n Abstand Erdnähe : %12.8f AE ',minrekbetr)
fprintf('\n Abstand Erdnähe : %12.4e km ',minrekbetr*AEk)
fprintf('\n Datum Erdnähe   : %s', dtErdNaehe)
fprintf('\n')
fprintf('\n Abstand Perihel : %12.8f AE ',minrskbetr)
fprintf('\n Abstand Perihel : %12.4e km ',minrskbetr*AEk)
fprintf('\n Datum Perihel   : %s', dtPerihel)
fprintf('\n')


%_________________________________________________________________________
% 

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------
