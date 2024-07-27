% -------------------------------------------------------------------------
% PlanetenPositionen.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet die Positionen der Planeten Venus bis Jupiter
% auf Basis der Keplergleichung und der Störungstheorie
% -------------------------------------------------------------------------

% Initialisierung
clc
clear 
close all 
addpath('./IncludeFolder/')
addpath('./Data/')

dt1 = datetime('1989-01-01 00:00:00');
dt1 = datetime('2004-06-08 10:00:00');
dt1 = datetime('1999-08-30 00:00:00');

T1 = juliandate(dt1); % Julianisches Datum  ET
t1 = Jd2JJht(T1);
Aequi = 'J2000';
Aequi = 'Datum';
[BaPa,BaPadot]=OrbitParameter(T1,Aequi);
fprintf('\n');
epsE = EpsErde(T1);
for k =1:9
    temp(k)=PlanetPQR(T1, BaPa, BaPadot, k);
end
Erde = temp(3);
fprintf('\n %s \n', string(dt1) );
for k =1:9
   Planets(k) = Convert2Equ(temp(k), Erde, deg2rad(epsE)); 
   fprintf('Planet: %s\t lambda:  %s  beta:  %s   Abstand  %12.4f \t RA: %s \t DEC: %s \t Abstand  %12.4f \n', ... 
          Planets(k).Name,StrDMS(wrapTo360(rad2deg(Planets(k).ekl(2)))),...
          StrDMS(rad2deg(Planets(k).ekl(3))), Planets(k).ekl(1), ...
          StrHMS(wrapTo360(rad2deg(Planets(k).equ(2)))/15),...
          StrDMS(rad2deg(Planets(k).equ(3))), Planets(k).equ(1));
end

tau =  8.32/1440/36525; % Lichtlaufzeit Sonne-Erde in JJht;
tauJup = Planets(5).ekl(1,1)*tau; % Lichtlaufzeit Jupiter-Erde;
tauVen = Planets(2).ekl(1,1)*tau; % Lichtlaufzeit Venus-Erde;
tauMars= Planets(4).ekl(1,1)*tau; % Lichtlaufzeit Mars-Erde;

VenData  = PerturbImport('VenusPos.csv');
SunData  = PerturbImport('SonnePos.csv');
MarsData = PerturbImport('MarsPos.csv');
JupData  = PerturbImport('JupiterPos.csv');
VenPos   = VenusExakt(t1-tauVen,VenData);
SunPos   = SonneExakt(t1,SunData,epsE,'AE');
MarsPos  = MarsExakt(t1-tauMars,MarsData);
JupPos   = JupiterExakt(t1-tauJup,JupData);

% PlanetPos(3) = SunPos;

PlanetPos(2) = VenPos;
PlanetPos(4) = MarsPos;
PlanetPos(5) = JupPos;

PlanetPos(3).xyz = -SunPos.xyz;
PlanetPos(3).ekl = CalcAnglesfromXYZ(PlanetPos(3).xyz);
PlanetPos(3).Name = 'Erde';
ErdePos  = PlanetPos(3);


for k=2:5 
    Planets2(k) = Convert2Equ(PlanetPos(k), ErdePos, deg2rad(epsE));
end


fprintf('\n %s \n', string(dt1) );
for k = 2:5
    fprintf('Planet: %s\t lambda:  %s  beta:  %s   Abstand  %12.4f \t RA: %s \t DEC: %s \t Abstand  %12.4f \n', ... 
            Planets(k).Name,StrDMS(wrapTo360(rad2deg(PlanetPos(k).ekl(2)))),...
            StrDMS(rad2deg(PlanetPos(k).ekl(3))), PlanetPos(k).ekl(1), ...
            StrHMS(wrapTo360(rad2deg(Planets2(k).equ(2)))/15),...
            StrDMS(rad2deg(Planets2(k).equ(3))), Planets2(k).equ(1));
end

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------