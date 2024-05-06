% -------------------------------------------------------------------------
% AnalemmaSonne.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Das vorliegende Programmaterial wurde von Kenneth von Buenau und Michael
% Kaschke erstellt.
% -------------------------------------------------------------------------
% Analemma der Sonne vom Mond
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

% Initialisierung Epoche und Daten
dt1 = datetime('2014-01-01 12:00:00');
T1  = juliandate(dt1);
t1  = Jd2JJht(T1);
epsE    = EpsErde(T1);  
epsErad = deg2rad(epsE); 
Aequi = 'Datum';
% Bahnparameter Erde/Sonne
[BaPa,BaPadot]= OrbitParameter(T1,Aequi);
exz    = BaPa.eP(3);         % Exzentrizität
PiqP   = BaPa.PiqP(3);       % heliozentr. Länge Perihel
M0P    = BaPa.M0P(3);        % mittlere Anomalie
n_trop = BaPadot.L0dot(3);   % Umlaufzeit
a = BaPa.aP(3)*149597871;    % Halbachse in km
LSP    = 282.9373;           % Fester Wert der Erde als Perihel 
                             % der effektiven Mondbahn genommen

% Datumsvektor fuer Zeitbereich
NrDays= 364;   % Erdsonnentage (Sols) für Berechnung
TF = T1+NrDays+1;
Tv          = linspace(T1,TF,NrDays+1);   %Jul. Tage
tv          = Jd2JJht(Tv);                %Julian. Jahrhunderte 
dtv(1:length(Tv)) = datetime(Tv,'convertfrom','juliandate');

% Monatsbeginn
index = 1;
for tag = 1:length(dtv)  %in Julianischen Tagen
    MonatsTag = day(dtv(tag),'dayofmonth');
    if MonatsTag == 1
      BeginnMonat(index) = tag;
      Monat(index) = month(dtv(tag));
      Jahr(index)  = year(dtv(tag))-2000;
      index = index+1;
    end
end

% Vollmondtage 2014,2020
VollmondStr = ...
['16-Jan-2014 04:51:55'; ...
 '14-Feb-2014 23:52:48' ; ...     
 '16-Mar-2014 17:08:11'; ...    
 '15-Apr-2014 07:42:12'; ...
 '14-May-2014 19:15:52'; ...    
 '13-Jun-2014 04:11:28'; ...    
 '12-Jul-2014 11:24:55'; ...    
 '10-Aug-2014 18:09:24'; ... 
 '09-Sep-2014 01:38:09'; ...    
 '08-Oct-2014 10:50:29'; ...  
 '06-Nov-2014 22:22:41'; ...    
 '06-Dec-2014 12:26:40'];    
VollMond     = datetime(VollmondStr,'InputFormat','dd-MMM-yyyy HHH:mm:ss');
VollMondTag  = day(VollMond,'dayofmonth');
VollMondMonat= month(VollMond,'monthofyear');
for index = 1:length(VollMond)  
    for tag = 1: length(dtv(1,:))
        MondTag  = day(dtv(1,tag),'dayofmonth');
        MondMonat  = month(dtv(1,tag),'monthofyear');
        if MondTag == VollMondTag(index) && ...
             MondMonat == VollMondMonat(index)
           IstVollMond(index) = tag;
        end
    end
end


Moon           = KeplerMond(Tv,epsErad);   % Mondkoordinaten in Näherung
MoonPos        = Moon.equ_xyz;             % Geozentr-aequatorial
[RAsun,DecSun] = KeplerSonneEx(Tv,epsErad,1); % Sonnenkoordinaten geo-aequator.
upsilon = RAsun - (PiqP-180)/180*pi;       % wahre Anomalie der Sonne
r = a*(1-exz^2)./(1+exz*cos(upsilon));     % Abstand Sonne Erde
Sun_sph = [r' RAsun' DecSun'];              
SunPos    = CalcXYZfromAngles(Sun_sph'); % Sonnenkoordinaten von Erde kart.
SunfromMoon = SunPos - MoonPos;          % Sonnenkorodinaten von Mond kart. aeq

for tag = 1:NrDays+1
    tz = tv(tag);
    % Berechnung der Rotationsachse des Mondes zum Zeitpunkt tz
    [alpha0,delta0,Sense] = RotAxis(tz);
    % Umrechnung in geozentrische-ekliptikale Koordinaten
    axis_b = asind(-cosd(delta0)*sind(alpha0)*sind(epsE)+ ...
             sind(delta0)*cosd(epsE));
    axis_l = wrapTo360(atan2d(cosd(delta0)*sind(alpha0)*cosd(epsE)+ ... 
             sind(delta0)*sind(epsE),cosd(delta0)*cosd(alpha0))); 
    phi    = (axis_l+90)/180*pi;
    lambda = (90-axis_b)/180*pi;
    
    rot1 = R_z(-phi);
    rot2 = R_x(+lambda);
    rot3 = R_z(+phi);
    
    rot = R_x(epsErad);  % Rotation Referenzebene Aequator -> Ekliptik
    SunfromMoonDay = rot*SunfromMoon(:,tag);
    
    SunAequ_xyz    = rot3*rot2*rot1*(SunfromMoonDay);
    SunAequ(tag,:) = CalcAnglesfromXYZ(SunAequ_xyz);
    
    axis = CalcXYZfromAngles([1,axis_l/180*pi,axis_b/180*pi]');
    eps(tag) = acosd(sum(axis .* [0 0 1]' )); % Achsneigung des Monds
end
% alpha Mittlere Sonne
alpha_FMS = atand(cosd(eps)*tand(LSP+M0P)); % Projektion auf Mondäquator
alphaMS = (alpha_FMS + n_trop*tv);
% alpha wahre Sonne
alpha_S = SunAequ(:,2)/pi*180;
% ZGL
ZGL = 4*(wrapTo180(alphaMS - alpha_S'));
Dec = SunAequ(:,3)/pi*180;


%%
% Graphische Ausgabe

% ZGL
figure(); 
plot(dtv,ZGL,'Color', Colors(4,:),'LineStyle','-','LineWidth',2);
hold on;
ylim([-10 +10]);
xlabel('Datum')
ylabel('ZGL in min');
title('Zeitgleichung Sonne vom Mond');
grid on;
set(gca,'Fontsize',16);

% Analemma
figure();
plot(ZGL, Dec,'-+','MarkerIndices',BeginnMonat,'LineWidth',2,...
      'Color',Colors(4,:),'LineStyle','-');
hold on;
plot(ZGL, Dec,':o','MarkerIndices',IstVollMond,'LineWidth',2,...
      'Color',Colors(10,:),'LineStyle',':');
plot(ZGL, Dec, 'Color', Colors(4,:),'LineStyle','-','LineWidth',2);

for k=1:length(BeginnMonat)
        mylabels(k,:)=sprintf('1.%02u.',Monat(k));
        xL(k)=ZGL(BeginnMonat(k));
        yL(k)=Dec(BeginnMonat(k));
        LabelPoints(xL(k), yL(k) ,mylabels(k,:),'e',0.3,0,'FontSize',12,...
                   'Color',Colors(4,:));
end
for k=1:length(VollMond)
    xL(k)=ZGL(IstVollMond(k));
    yL(k)=Dec(IstVollMond(k));
end
plot(xL,yL,'Color',Colors(10,:),'Linestyle',':','Linewidth',2);
grid on
ylim([-3 3]);
xlim([-10 10]);
title('Analemma Sonne vom Mond');

% for k=1:length(ZGL)
% % clf
% hold on
%  p(2) = plot(ZGL(k),Dec(k));
%  p(2).Color = Colors(10,:);
%  p(2).LineWidth = 4;
%  p(2).Marker = 'o';
%  p(2).Visible = 'on';
%  p(3) = text(15,20,string(dtv(k),'dd-MMM-yyyy'));
%  p(3).FontSize = 12;
%  p(3).Visible = 'on';
%  pause(0.01)
%  p(2).Visible = 'off';
%  p(3).Visible = 'off';
% end
legend(strcat(string(dtv(1),'dd-MMM-yyyy'),' - ',...
        string(dtv(length(dtv)),' dd-MMM-yyyy')),...
        'Vollmondtage','location','southwest');
legend box off;

ylabel('Deklination in °')
xlabel('ZGL in min');
set(gca,'Fontsize',16);


%----------------------------------------------------------------------
% RotAxis: 
%      Berechnet die Orientierung des planetozentrischen Bezugssystems +
%      liefert die Rotationsrichtung
% Eingabe:
%   t         Zeit in Julianischen Jahrhunderten seit J2000
% Ausgabe:
%   alphaA,
%   deltaA    Richtunsgvektroren der Polachse des Monds
%   Sense     Rotationsrichtung (direkt oder retrograd)
%--------------------------------------------------------------------------
function [alphaA, deltaA, Sense] = RotAxis(t)
  % Position der Mondachse in geozentrisch-aequatorialen Koordinaten
  % 

  d = t*36525;
  
  E1 = 125.045-0.0529921*d;
  E2 = 250.089-0.1059842*d;
  E3 = 260.008+13.0120009*d;
  E4 = 176.625 + 13.3407154*d;
  E6 = 311.589 + 26.4057084*d;
  E7 = 134.963 + 13.0649930*d;
  E10 = 15.134 - 0.1589763*d;
  E13 =  25.053 + 12.9590088*d;
  
  RA = 269.9949 + 0.0031*t - 3.8787*sind(E1) - 0.1204*sind(E2)...
      +0.0700*sind(E3) - 0.0172 * sind(E4) + 0.0072* sind(E6)...
      -0.0052*sind(E10) + 0.0043*sind(E13);
  
  Dec =  66.5392 + 0.0130*t + 1.5419*cosd(E1) + 0.0239*cosd(E2)...
      - 0.0278*cosd(E3) + 0.0068*cosd(E4) - 0.0029*cosd(E6)...
      +0.0009*cosd(E7) + 0.0008*cosd(E10) - 0.0009*cosd(E13);
  alphaA  = RA; 
  deltaA  = Dec; 
  % Rotationsrichtung (spielt keine Rolle bei unseren Rechnungen)
  Sense = -1; 
end 
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------