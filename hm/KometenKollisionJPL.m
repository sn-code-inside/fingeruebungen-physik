%%-------------------------------------------------------------------------
% KometenKollisonJPL.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Astrodynamik" aus
% "FingerÃ¼bungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% 
% Darstellung der Bewegung des Kometen Swift-Tuttle 
% durch JPL Daten 
%
%
% -------------------------------------------------------------------------
clc
clear 
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Style = ["-", "-.", ":", "--", ":"];
Colors=GetColorLines;

AEk   = 1.495978707e08;         % AE in km

% Einlesen der JPL Files
E  = importfile("Erdbahn2126.csv", [2, Inf]);
ST = importfile("SwiftTuttle2126.csv", [2, Inf]);
ST.x = ST.x/AEk;
ST.y = ST.y/AEk;
ST.z = ST.z/AEk;


% % Zeitraum und Darstellung Perihelnaehe Swift-Tuttle
% % fuer 2126
T1 = ST.T(1);
T2 = ST.T(end);
dt1 = datetime(ST.dT(1));
dt2 = datetime(ST.dT(end));
dt1Pstr= datestr(dt1);
dt2Pstr= datestr(dt2);

E.x = E.x/AEk;
E.y = E.y/AEk;
E.z = E.z/AEk;
xAx =1.5;


% Einlesen der Bahnparameter der Planeten und deren Ableitungen 
Aequi = 'Date';
[BaPa,BaPadot]=OrbitParameter(T1, Aequi); 

% Berechnung fuer kurze Zeitraeume 
long=0;
Nges=10000; %Zahl der Stuetzstellen im Intervall T1 - T2
T=linspace(T1,T2,Nges); %Berechnungsintervall
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

TimeStep = 30;
XLabelK = ST.x(1:TimeStep:end);
YLabelK = ST.y(1:TimeStep:end);
ZLabelK = ST.z(1:TimeStep:end);
XLabel = E.x(1:TimeStep:end);
YLabel = E.y(1:TimeStep:end);
ZLabel = E.z(1:TimeStep:end);

dtKstr = datestr(datetime(E.T(1:TimeStep:end),'ConvertFrom',...
                'juliandate'),' dd.mm.');

% Berechnung der heliozentrischen Koordinaten der Erde und 
% der inneren Planeten von T1 bis T2 
for k=1:2
    PlanetKep(k)=PlanetPQR(T, BaPa, BaPadot,k);
end

%Berechnung Abstand Erde-Komet
rek=[E.x-ST.x, E.y-ST.y, E.z-ST.z];
for k= 1:length(rek)
    rekbetr(k) = sqrt(rek(k,1)^2+rek(k,2)^2+rek(k,3)^2);
end
[minrekbetr,ErdNaehe] = min(rekbetr);

%Berechnung Abstand Sonne-Komet und Perihel
rsk=[ST.x, ST.y, ST.z];
for k= 1:length(rsk)
    rskbetr(k) = sqrt(ST.x(k)^2+ST.y(k)^2+ST.z(k)^2);
end
[minrskbetr,Perihel] = min(rskbetr);

dtErdNaehe = datestr(datetime(ST.T(ErdNaehe),'convertfrom','juliandate'),...
               ' dd.mm.yy ');
dtPerihel  = datestr(datetime(ST.T(Perihel),'convertfrom','juliandate'),...
               ' dd.mm.yy ');


%% Graphik 3D

figure()
plot3(E.x,E.y,E.z,'Color', Colors(3,:),'Linewidth',2)
hold on
p(1)=plot3(ST.x,ST.y,ST.z,'Color', Colors(9,:),'Linewidth',2);
SonneFP;
for iPlot= 1:2
  hp(iPlot) = plot3(PlanetKep(iPlot).xyz(1,:), ...
      PlanetKep(iPlot).xyz(2,:), ...
      PlanetKep(iPlot).xyz(3,:),'Color', Colors(iPlot,:), ...
      'linewidth',2);
end
grid on
xlim([-xAx xAx]);
ylim([-xAx xAx]);
xl = xlim;
yl = ylim;
[X,Y] = meshgrid(xl,yl);
surf(X,Y,zeros(size(X)));
shading flat
alpha 0.1
hold on
axis equal
% Achse Sonne-Fruehlingspunkt 
SonneFP;
xlim([-xAx xAx]);
ylim([-xAx xAx]);
zlim([-xAx xAx]);

%Zeitdaten
plot3(XLabelK, YLabelK, ZLabelK,'+', ...
     'color',Colors(9,:),'markerfacecolor',Colors(9,:),...
     'markersize',8,'linewidth',2);
plot3(XLabel, YLabel, ZLabel,'+', ...
     'color',Colors(3,:),'markerfacecolor',Colors(8,:),...
     'markersize',8,'linewidth',2);
text(XLabel, YLabel, ZLabel, dtKstr, 'FontSize',9, 'Color', Colors(3,:));
text(XLabelK, YLabelK, ZLabelK, dtKstr, 'FontSize',9, 'Color', Colors(9,:));

p(2)= line([E.x(ErdNaehe) ST.x(ErdNaehe)],...
     [E.y(ErdNaehe) ST.y(ErdNaehe)],...
     [E.z(ErdNaehe) ST.z(ErdNaehe)],...
     'color', Colors(4,:), 'Linewidth',3);
text(ST.x(ErdNaehe), ...
     ST.y(ErdNaehe), ST.z(ErdNaehe),...
     dtErdNaehe, 'FontSize',9,'Color', Colors(4,:));

plot3(ST.x(Perihel), ST.y(Perihel),ST.z(Perihel), ...
    'o','Color', Colors(10,:),'MarkerSize',5, ...
    'MarkerFaceColor',Colors(10,:));
p(3)=line([0 ST.x(Perihel)], [0 ST.y(Perihel)], [0 ST.z(Perihel)],...
     'Color', Colors(10,:),'Linewidth',2);
plot3(ST.x(Perihel),ST.y(Perihel),ST.z(Perihel), ...
    'o','Color', Colors(10,:),'MarkerSize',5, ...
    'MarkerFaceColor',Colors(10,:));
text(ST.x(Perihel), ...
     ST.y(Perihel), ST.z(Perihel),...
     dtPerihel, 'FontSize',9,'Color', Colors(10,:));

% Einstellungen der Achsen 
xlim([-40 40]);
ylim([-40 40]);
xl = xlim;
yl = ylim;
[X,Y] = meshgrid(xl,yl);
surf(X,Y,zeros(size(X)));
legend(p,'Swift-Tuttle','Erdnähe','Perihel','location','best')
legend box off
zlim([-40 40]);
shading flat
alpha 0.1
grid on;
axis equal;
AEax = 2;
header2 = strjoin([DatumStr1,' bis ',DatumStr2]);
header2 = strjoin(['Swift-Tuttle: ',header2]);

ht=title(header2);
set(ht,'FontSize',12, 'FontWeight', 'normal');
xlim([-AEax AEax]);
ylim([-AEax AEax]);
zlim([-AEax AEax]);
xlabel('x in AE')
ylabel('y in AE');
zlabel('z in AE');
set(gca,'FontSize',14);

%% Ausgabe


fprintf('\n')
fprintf('\n')
fprintf('\n Berechnung für 2126 über JPL-Daten')
fprintf('\n')
fprintf('\n Abstand Erdnähe : %12.8f AE ',minrekbetr)
fprintf('\n Abstand Erdnähe : %12.4e km ',minrekbetr*AEk)
fprintf('\n Datum Erdnähe   : %s', dtErdNaehe)
fprintf('\n')
fprintf('\n AbstandPerihel  : %12.8f AE ',minrskbetr)
fprintf('\n AbstandPerihel  : %12.4e km ',minrskbetr*AEk)
fprintf('\n Datum Perihel   : %s', dtPerihel)
fprintf('\n')


%%

function Bahn = importfile(filename, dataLines)
%  Example:
%  Erdbahn2126 = importfile("\Erdbahn2126.csv", [2, Inf]);
% Input handling
% If dataLines is not specified, define defaults
    if nargin < 2
        dataLines = [2, Inf];
    end
    % Set up the Import Options and import the data
    opts = delimitedTextImportOptions("NumVariables", 5, "Encoding", "UTF-8");
    % Specify range and delimiter
    opts.DataLines = dataLines;
    opts.Delimiter = ";";
    % Specify column names and types
    opts.VariableNames = ["T", "dT", "x", "y", "z"];
    opts.VariableTypes = ["double", "datetime", "double", "double", "double"];
    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    % Specify variable properties
    opts = setvaropts(opts, "dT", "InputFormat", "yyyy-MMM-dd HH:mm:ss");
    % Import the data
    Bahn = readtable(filename, opts);
end
