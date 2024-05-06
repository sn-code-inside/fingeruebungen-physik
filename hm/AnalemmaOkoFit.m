% -------------------------------------------------------------------------
% AnalemmaOkoFit.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Analemma Oberkochen Fit im Horizontal-KOS
% Umrechnung inkl. Bildkorrektur des aufgenommenen Analemma von 
% Oberkochen 2019/20 in das horizontale KOS, 
% sowie Darstellung verschiedener theoretisch berechneter Analemma. 
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

% Initialisierung
lambda = 10.0990193;       % Geographische Länge Oberkochen +-BB.bb  
phi    = 48.7854368;       % Geographische Breite Oberkochen +-BB.bb  

% Parameter für Bildanpassung
Az0  = DD([117 30 19]);  % Zentralwert des Azimuth (aus Eichung)
Alt0 = DD([25 20 23]);   % Zentralwert der Hoehe (aus Eichung)
scalex =1   ;            % Skalenfaktor x-Achse
%Parameterme für Ausgleichsrechnung (separat bestimmt)
M1 = [1.5389,-0.0445, -61.3737];
M2 = [0.0502, 1.6505, -29.0163];

exz  = 2;    %doppelte Exzentrizität
exz0 = 0.01617;

%%
% Vorwärtsrechnung
dt1 = datetime('2019-03-31 07:40:00');
T = juliandate(dt1);
MJuDa = juliandate(dt1,'modifiedjuliandate');
DatumStr1 =  string(dt1,'dd.MM.yyyy');
T0=(T-2451545)/36525;
eps=23.43929111-(46.8150+0.00059*T0-0.001813*T0*T0)*T0/3600;
pi2 = 2*pi;
T_v=T:(T+365);
% DateBer = string(datetime(T_v,'ConvertFrom','juliandate'),'dd.MM.yyyy');
MJuDa_v = MJuDa:(MJuDa+365);
Timey = 1:366;
Null = zeros(1,366);
% Berechnung nach Keplerlösung
[RA_v,Dec_v]  = KeplerSonneEx(T_v,deg2rad(eps),1.0); 
[RA_v1,Dec_v1]= KeplerSonneEx(T_v,deg2rad(eps),exz); %doppelte Exzentrizität
% Berechnung des Stundenwinkels
tau_v  = rad2deg((GMST(MJuDa_v)-RA_v))+lambda;
tau_v1 = rad2deg((GMST(MJuDa_v)-RA_v1))+lambda;
RA_v   = wrapTo360(rad2deg(RA_v));
RA_v1  = wrapTo360(rad2deg(RA_v1));
Dec_v  = rad2deg(Dec_v);
Dec_v1 = rad2deg(Dec_v1);
tau_v  = 360*wrapToPi(pi2*tau_v/360)/pi2; 
tau_v1 = 360*wrapToPi(pi2*tau_v1/360)/pi2; 
%Berechnung Azimut und Höhe
[Az_v, Alt_v]  = Equ2Hor(RA_v,Dec_v,phi,tau_v);
[Az_v1, Alt_v1]= Equ2Hor(RA_v1,Dec_v1,phi,tau_v1);
%Korrektur Refraktion
Alt_vk=Alt_v+1.02./tand(Alt_v+10.3./(Alt_v+5.11))/60;
Alt_v1=Alt_v1+1.02./tand(Alt_v1+10.3./(Alt_v1+5.11))/60;
Az_v = Az_v +180;
Az_v1 = Az_v1 +180;

%% 
% Image des Analemma
% Einlesen von File
opts = delimitedTextImportOptions("NumVariables", 8);
opts.DataLines = [1, Inf];
opts.Delimiter = ";";
% Specify column names and types
opts.VariableNames = ["Datum", "X", "Y", "Nr"];
opts.SelectedVariableNames = ["Datum", "X", "Y"];
opts.VariableTypes = ["string", "double", "double","double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
% Import the data
AnaOko = readtable("AnalemmaOko2019.csv", opts);
clear opts

%%
xPix = 623.147;
yPix = 415.883;
XDotP =  AnaOko.X;
YDotP = -AnaOko.Y;

% Berechne die Bildkoordinaten im Standardkoordinaten
xZero = xPix/2;
yZero = yPix/2;
XDot = (XDotP - xZero)*scalex/xPix;
YDot = (YDotP + yZero)/xPix;  %durch xPix oder yPix ?

% Jetzt direkte Umrechnung in Horizontalkoordinatensystem
[Az1, Alt1]= Bild2Equ(Az0,Alt0, -XDot, YDot);

%Ausgleichsrechnung
% DateMess =string(AnaOko.Datum(:)');
% [M1,M2]=Ausgleich(DateBer,DateMess, Az_v,Az1, Alt_vk, Alt1,n1,n2,n3)
Az = M1(1)*Az1+M1(2)*Alt1+M1(3);
Alt = M2(1)*Az1+M2(2)*Alt1+M2(3);

%------------------------------------------------------------------------------
% Graphische Ausgabe

figure();
header1='Analemma Simulation'; 
%Messwerte
plot(wrapTo360(Az),Alt,'bo',...
    'LineWidth',1,...
    'MarkerSize',5,...
    'MarkerFaceColor',Colors(3,:))              
hold on;
%ber. mit Refraktion
plot(Az_v,Alt_vk,'Color',Colors(4,:),'LineWidth',2, 'LineStyle','-'); 
%ber. andere Exzentrizität
plot(Az_v1,Alt_v1,'Color',Colors(2,:),'LineWidth',2,'LineStyle',':');
legend('Messwerte korrigiert', 'theoretisch \it{e} \rm = 0.01617', 'theoretisch \it{e} \rm = '+string(num2str(2*exz0)));
legend boxoff;
xlim([80 160]);
grid on;
xlabel('Azimut °');
ylabel('Hoehe °');
set(gca,'XGrid','on', 'YGrid', 'on', 'Fontsize', 18, 'linewidth', 1);

%%

% Berechne Deklination, RA aus Bildkoordinaten
% Analog auch für Alt Az !!!!
function [RA, Dec]= Bild2Equ(RA0,Dec0, X, Y) 
% Eingabe:
%   RA0       RA der optischen Achse bze. Az
%   Dec0      Deklination der optischen Achse bzw. Hoehe
%   X         Standardkoordinate X [dimensionslos]
%   Y         Standardkoordinate Y [dimensionslos]
% Ausgabe:
%   RA        RA
%   Dec       Dec
% Beachte: alle Winkel in [Grad]
% ------------------------------------------------------------------------------
  RA = RA0+rad2deg(atan(-X./(cosd(Dec0)-Y.*sind(Dec0))));
  Dec =(asind((sind(Dec0)+ Y.*cosd(Dec0))./sqrt(1+X.*X+Y.*Y)));
end

% Berechne Alt-Az aus Deklination, RA 
function [Az, Alt]= Equ2Hor(RA,Dec,phi,tau) 
% Eingabe:
%   RA       RA 
%   Dec      Deklination 
%   phi      Breite
%   tau      Stundenwinkel
% Ausgabe:
%   Az        Azimuth
%   Alt       Altitude
%% Beachte: alle Winkel in [Grad]
% ------------------------------------------------------------------------------
  Alt=asind(sind(phi).*sind(Dec)+cosd(phi).*cosd(Dec).*cosd(tau));
  Az= asind(cosd(Dec).*sind(tau)./cosd(Alt));
end
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------