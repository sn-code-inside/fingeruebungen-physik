% -------------------------------------------------------------------------
% PositionSonneVonPlanet.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Das vorliegende Programmaterial basiert in Teilen 
% auf C++ Programmstrukturen und der Beschreibung/Erläuterung in  
%
% "Astronomie mit dem Personal Computer"
% von Oliver Montenbruck und Thomas Pfleger (Springer 1999). 
% Genehmigung des Verlages und der Autoren liegt vor.
% -------------------------------------------------------------------------
% Programm berechnet die planetenzentrischen Positionen der
% Sonne von den verschiedenen Planeten aus 
% auf Basis der Koordinaten aus der Keplergleichung und der 
% Störungstheorie (Venus bis Jupiter)
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

% Initialisierung
c_light   = 1/0.00578;   % AE/d

% dt1 = datetime('2004-06-08 15:00:00');
dt1 = datetime('1999-08-30 00:00:00');
dt1 = datetime('2000-03-20 06:00:00');
T1 = juliandate(dt1); % Julianisches Datum  ET
t1 = Jd2JJht(T1);
epsErad = deg2rad(EpsErde(T1));
Aequi = 'J2000';
Aequi = 'Datum';
[BaPa,BaPadot]=OrbitParameter(T1,Aequi);

for k =1:9
    Planets(k)=PlanetPQR(T1, BaPa, BaPadot, k);
end
Erde = Planets(3);
fprintf('\n %s Keplerlösung \n', string(dt1) );
fprintf(' Planet \t \t longitude \t \t \t  latitude  \t  \n');
r_Earth       = mtimes(R_x(-epsErad),Erde.xyz);
for PlanetNr = 1:9
    r       = mtimes(R_x(-epsErad),Planets(PlanetNr).xyz);
    Delta   = norm(r-r_Earth);
    [R_equ, f] = Shape (PlanetNr);
    %Lichtlaufzeitkorrektur
    t2(PlanetNr)      =  t1 - Delta/c_light/36525; 
    tkorr             =  t1 - Delta/c_light/36525; 
    Planets(PlanetNr)=PlanetPQR(JJht2Jd(tkorr),BaPa,BaPadot,PlanetNr);
    r       = mtimes(R_x(-epsErad),Planets(PlanetNr).xyz);
    r_geoc = r - r_Earth;
    [E, Sense]  = Orient (PlanetNr,tkorr);
    [lonS, latS, lat2S] = Rotation (-r_geoc,E,Sense,f);
    fprintf(' %s \t \t %s \t  %s  \n', ... 
                    Planets(PlanetNr).Name,...
                    StrDMS(wrapTo360(rad2deg(lonS))),...
                    StrDMS(rad2deg(latS)));
end 

%Exakte Koordinaten aus der Störungstheorie
SunData  = PerturbImport('SonnePos.csv');
VenData =  PerturbImport('VenusPos.csv');
MarsData = PerturbImport('MarsPos.csv');
JupData  = PerturbImport('JupiterPos.csv');
SunPos   = SonneExakt(t1,SunData,epsErad,'AE');

fprintf('\n %s Störungstheorie \n', string(dt1) );
fprintf(' Planet \t \t longitude \t \t \t  latitude  \t  \n');
r_Earth = mtimes(R_x(-epsErad),Erde.xyz);
Delta   = norm(r_Earth);

for PlanetNr = 2:5
    if PlanetNr  == 3
        r_ex       = SunPos.xyz;
    else
        switch PlanetNr
            case 2 
                r   = VenusExakt(t2(PlanetNr),VenData);
            case 4 
                r   = MarsExakt(t2(PlanetNr),MarsData);
            case 5 
                r   = JupiterExakt(t2(PlanetNr),JupData);
        end
        r_ex       = mtimes(R_x(-epsErad),r.xyz);
    end
    % Geometrische geozentrische Position Planet
    Delta   = norm(r_ex);
    [E, Sense]  = Orient (PlanetNr,t2(PlanetNr));
    [lonS, latS, lat2S] = Rotation (-r_ex,E,Sense,f);
    fprintf(' %s \t \t %s \t  %s  \n', ... 
                    Planets(PlanetNr).Name,...
                    StrDMS(wrapTo360(rad2deg(lonS))),...
                    StrDMS(rad2deg(latS)));
end
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------



%--------------------------------------------------------------------------
% Shape: Liefert den Aequatorradius und die Abplattung der Planeten
% Eingabe:
%   PlanetNr  Identifiziert den Planeten 
% Ausgabe:
%   R_equ     Aequatorradius in [km]
%   f         Geometrische Abplattung
%--------------------------------------------------------------------------

function [R_equ, f] = Shape (PlanetNr)
 switch (PlanetNr) 
    case 10    
        R_equ = 696000.0; f = 0.0;       
    case 1     
        R_equ =  2439.0; f = 0.0;       
    case 2     
        R_equ =  6051.0; f = 0.0;        
    case 3
        R_equ =  6378.14; f = 0.00335281; 
    case 4
        R_equ =  3393.4;  f = 0.0051865;  
    case 5
        R_equ = 71398.0;  f = 0.0648088;  
    case 6
        R_equ = 60000.0;  f = 0.1076209;  
    case 7
        R_equ = 25400.0;  f = 0.030;      
    case 8
        R_equ = 24300.0;  f = 0.0259;     
    case 9
        R_equ =  1500.0;  f = 0.0;
 end
end

%--------------------------------------------------------------------------
% Orient: Berechnet die Orientierung des planetozentrischen Bezugssystems +
%         liefert die Rotationsrichtung
% Eingabe:
%   PlanetNr  Identifiziert den Planeten (s. APC_Planets.h)
%   t         Zeit in Julianischen Jahrhunderten seit J2000
% Ausgabe:
%   ET        Transformationsmatrix vom geozentrischen, auf den mittleren 
%             Aequator und das Aequinoktium J2000 bezogenen KOS
%             zum koerperfesten aequatorialen Bezugssystem
%   Sense     Rotationsrichtung (direkt oder retrograd)
%--------------------------------------------------------------------------
function [ET, Sense] = Orient (PlanetNr,t)
  % Rektaszension und Deklination der Rotationsachse bezogen auf
  % Aequator und Ekliptik zur Epoche J2000; Orientierung des
  % Nullmeridians
  d = t*36525;  % Tage seit J2000
  switch (PlanetNr)
    case 10
        RA  = 286.13; Dec =  63.87; W   =  84.182 +  14.1844000*d;  
    case 1
        RA  = 281.01  -   0.033*t;  Dec =  61.45  -   0.005*t;    
        W   = 329.68  +   6.1385025*d;  
    case 2
        RA  = 272.76; Dec =  67.16; W   = 160.20  -   1.4813688*d;  
    case 3
        RA  =   0.00  -   0.641*t;  Dec =  90.00  -   0.557*t;
        W   = 190.16  + 360.9856235*d;  
    case 4
        RA  = 317.681 -   0.108*t; Dec =  52.886 -   0.061*t;
        W   = 176.901 + 350.8919830*d;  
    case 5
        RA  = 268.05  -   0.009*t; Dec =  64.49  +   0.003*t;
        W =  67.10  + 877.900*d;  
    case 6
        RA  =  40.589 -   0.036*t; Dec =  83.537 -   0.004*t;  
        W   = 227.2037 + 844.3*d;
    case 7
        RA  = 257.311; Dec = -15.175;
        W   = 203.81  - 501.1600928*d;   
    case 8
        N   = deg2rad(357.85+52.316*t);    
        RA  = 299.36  + 0.70*sin(N); 
        Dec =  43.46  - 0.51*cos(N); 
        W   = 253.18  + 536.3128492*d - 0.48*sin(N);  
    case 9
        RA  = 313.02; Dec =   9.09;
        W   = 236.77  -  56.3623195*d;
  end
  W   = wrapTo360(W);
  RA  = deg2rad(RA); 
  Dec = deg2rad(Dec); 
  W   = deg2rad(W);  
  % Transformation vom mittleren Erdaequator und Fruehlingspunkt J2000
  % zum koerperfesten Bezugssystem von Aequator und Nullmeridian
  ET = PQR (W, pi/2-Dec, pi/2+RA);
  % Rotationsrichtung
  if ( PlanetNr==2 || PlanetNr==7 || PlanetNr==10 )
    Sense = +1;
  else
    Sense = -1;
  end
end 
 

%--------------------------------------------------------------------------
% Rotation: Berechnet die planetographischen Koordinaten
% Eingabe:
%   r        Planetozentrische Koordinaten der Erde 
%            bezogen auf Erdaequator 
%            und Aequinoktium J2000
%   E        Transformationsmatrix
%   Sense    Rotationsrichtung (direkt oder retrograd)
%   f        Geometrische Abplattung
% Ausgabe:
%   lon      Planetographische Laenge der Erde in [rad]
%   lat_g    Planetographische Breite der Erde in [rad]
%   lat_c    Planetozentrische Breite der Erde in [rad]
%--------------------------------------------------------------------------

function [lon, lat, lat2] =  Rotation (r, E, Sense, f)
  % Planetozentrische Breite und Laenge
  s     = mtimes(E,r);
  plako = CalcAnglesfromXYZ(s);
  lon   = Sense*plako(2); 
  lon   = Modulo(lon,2*pi);
  lat   = plako(3);
  lat2  = atan(tan(lat)/((1.0-f)*(1.0-f)));
end

