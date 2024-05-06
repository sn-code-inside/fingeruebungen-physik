% -------------------------------------------------------------------------
% PositionPlanet.m
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
% Programm berechnet die planetenzentrischen Positionen von
% Erde (oder Sonne) vom Mars auf Basis der Keplergleichung und 
% der der Störungstheorie
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

% Initialisierung
AE = 149597870;
Arcs    = 3600*180/pi;
dt1 = datetime('2004-06-08 10:00:00');
dt1 = datetime('1999-08-30 00:00:00');
T1 = juliandate(dt1); % Julianisches Datum  ET
t1 = Jd2JJht(T1);
epsE = EpsErde(T1);

% Heliozentrische Position der Erde, mittlerer Aequator und Aequinoktium
% des Datums


%Einlesen und Berechnung genauer Koordinaten
SunData  = PerturbImport('Data/SonnePos.csv');
MarsData = PerturbImport('Data/MarsPos.csv');
SunPos   = SonneExakt(t1,SunData,epsE,'AE');
MarsPos  = MarsExakt(t1,MarsData);
tau =  8.32/1440/36525; %Lichtlaufzeit Sonne-Erde in JJht;

ErdePos.Name = 'Erde';
ErdePos.xyz = -SunPos.xyz;
ErdePos.ekl = CalcAnglesfromXYZ(ErdePos.xyz);

Mars = Convert2Equ(MarsPos, ErdePos, deg2rad(epsE));
Erde = Convert2Equ(ErdePos, ErdePos, deg2rad(epsE));

r_Earth = Erde.geo;

% Wenn ich die planetozentrischen Sonnenkorodinaten haben will muss ich meE
% hier die Sonnenkoordinaten eingeben oder ????
% r_Earth = r_Sun = SunPos.xyz;
% Stimmt aber nicht, da am 8.6.2004 (siehe Abb. 3.28) die lon von Sonne und
% Erde ziemlich unterschiedlich sein muss. 
% Finden Sie den Fehler ?


PlanetNr = 4; % Mars

% Geometrische geozentrische Position und Lichtlaufzeit
r = Mars.geo;
Delta = norm(r-r_Earth);
tauP  = Delta*tau; %Lichtlaufzeit Mars-Erde;
t1 = t1 - tauP;
% Planetenposition zum Zeitpunkt der Lichtaussendung
MarsPos = MarsExakt(t1,MarsData);
Mars = Convert2Equ(MarsPos, ErdePos, deg2rad(epsE));
r = Mars.geo;
% Geozentrische Position
r_geoc = r - r_Earth;
% Scheinbarer Halbmesser (in ["])
[R_equ,f] =Shape (PlanetNr);
D_equ = Arcs*2.0*asin(R_equ/(Delta*AE));
% Transformation vom mittleren Aequator und Aequinoktium des
% Datums zum koerperfesten Bezugssystem
[E, Sense]       = Orient (PlanetNr,t1);
[lon, lat, lat2] = Rotation (r,E,Sense,f);
fprintf('Mars: \t lon:  %s  lat:  %s   lat  %s  \n', ... 
            StrDMS(wrapTo360( rad2deg(lon) ) ),...
            StrDMS(rad2deg(lat)), StrDMS(rad2deg(lat2)));

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Shape: Liefert den Aequatorradius und die Abplattung der Planeten
% Eingabe:
%   PlanetNr  Identifiziert den Planeten (s. APC_Planets.h)
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
%   E         Transformationsmatrix vom geozentrischen, auf den mittleren 
%             Aequator und das Aequinoktium J2000 bezogenen KOS
%             zum koerperfesten aequatorialen Bezugssystem
%   Sense     Rotationsrichtung (direkt oder retrograd)
%--------------------------------------------------------------------------
function [E, Sense] = Orient (PlanetNr,t)
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
    case 7
        RA  = 257.311; Dec = -15.175;
        W   = 203.81  - 501.1600928*d;   
    case 8
        N   = Rad*(357.85+52.316*t);    
        RA  = 299.36  + 0.70*sin(N); 
        Dec =  43.46  - 0.51*cos(N); 
        W   = 253.18  + 536.3128492*d - 0.48*sin(N);  
    case 9
        RA  = 313.02; Dec =   9.09;
        W   = 236.77  -  56.3623195*d;
  end

  RA  = deg2rad(RA); 
  Dec = deg2rad(Dec); 
  W   = deg2rad(Modulo(W,360.0));  
  
  % Transformation vom mittleren Erdaequator und Fruehlingspunkt J2000
  % zum koerperfesten Bezugssystem von Aequator und Nullmeridian
  E = PQR (W, pi/2-Dec,pi/2+RA);
  % Rotationsrichtung
  if ( PlanetNr==2 || PlanetNr==7 || PlanetNr==10 )
    Sense = 1;
  else
    Sense = -1;
  end
end 
 

%--------------------------------------------------------------------------
% Rotation: Berechnet die planetographischen Koordinaten
% Eingabe:
%   r        Planetozentrische Koordinaten der Erde bezogen auf Erdaequator 
%            und Aequinoktium J2000
%   E        Transformationsmatrix vom geozentrischen, auf den mittleren 
%            Aequator und das Aequinoktium J2000 bezogenen KOS
%            zum koerperfesten aequatorialen Bezugssystem
%   Sense    Rotationsrichtung (direkt oder retrograd)
%   f        Geometrische Abplattung
% Ausgabe:
%   lon      Planetographische Laenge der Erde in [rad]
%   lat_g    Planetographische Breite der Erde in [rad]
%   lat_c    Planetozentrische Breite der Erde in [rad]
%--------------------------------------------------------------------------

function [lon, lat, lat2] =  Rotation ( r, E, Sense, f)
  % Planetozentrische Breite und Laenge
  s =   mtimes( E,r);

  lat = s(3);
  lon   = Sense*s(2);
  lon = Modulo(lon,2*pi);
  lat2 = atan(tan(lat)/((1.0-f)*(1.0-f)));
end



%------------------------------------------------------------------------------
% Illum: Berechnet Beleuchtungsparameter der Planeten
% Eingabe:
%   r         Heliozentrischer Ort des Planeten
%   r_Earth   Heliozentrischer Ort der Erde
% Ausgabe:
%   Elong     Elongation in [rad]
%   phi       Phasenwinkel in [rad]
% % Beachte:    Koordinatensystem und Epoche fuer r und r_Earth muessen 
%             uebereinstimmen
%------------------------------------------------------------------------------
function  [Elong, phi] = Illum (r, r_Earth)
  % Geozentrischer Ort des Planeten
  r_geoc = r - r_Earth;
  % Entfernungen
  R  = norm(r);         % Entfernung Sonne-Planet
  RE = norm(r_Earth);   % Entfernung Sonne-Erde
  D  = norm(r_geoc);    % Entfernung Erde-Planet
  % Elongation, Phasenwinkel und Phase
  Elong = acos ((D*D + RE*RE - R*R) / (2*D*RE) );
  c_phi = (D*D + R*R - RE*RE) / (2*D*R);
  phi   = acos (c_phi);
end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------
