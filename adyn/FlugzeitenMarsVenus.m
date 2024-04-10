% -------------------------------------------------------------------------
% FlugzeitenMarsVenus.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Astrodynamik" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet die Flugzeiten fÃ¼r Hohmann-Ellipsen zu einem der 
% Nachbar-Planeten (Venus, Mars).
% 
% Näherung: Kreisförmige Bahnen der Planeten.
%
% -------------------------------------------------------------------------
%
clc
clear all
close all 

% Gl. 3.29
%% Mars
aE = 149598022e3;% m
aM = 227990000e3; % m
aH = (aE+aM)/2;
G = 6.67430e-11; % m3 kg/s2
mS = 1.989e30; % kg
muS = G*mS;
tH = pi*sqrt(aH*aH*aH/muS);
tHd = tH/60/60/24;
nM = 19140.563; % deg/Jh
nE = 35999.371; % deg/Jh
nMSI = nM * pi/180 /100/365.25/60/60/24;
nESI = nE * pi/180 /100/365.25/60/60/24;
psi_0 = pi - nMSI*tH;
psi_0deg = psi_0* 180/pi;
psi_F = pi - nESI*tH;
psi_Fdeg = psi_F* 180/pi;
tMW = (2*psi_F - 2*(-1)*pi)/(nESI-nMSI);
tMWd = tMW / 60/60/24;
tMgesd = tMWd + 2*tHd;

ttlprint = "Flugzeiten zum Mars (Tabellenwerte)";
fprintf('\n %s', ttlprint);
fprintf('\n\n')
fprintf('|  aE (m)    | aM (m)     | nE (rad/d) | nM (rad/d) |   tH (d)   |\n');
fprintf('| %7.4e | %7.4e | %7.4e | %7.4e |  %7.4f  |\n',...
                    aE, aM, nESI*86400, nMSI*86400, tHd);
fprintf('\n')
fprintf('|   psi0(°)  |  psi_F(°) |   tMW (d)   |  2*tH (d)  |  tMges (d) |\n');
fprintf('|  %7.3f   |  %7.3f  |   %7.3f   |  %7.3f   |  %+7.3f  |\n',...
                    psi_0deg, psi_Fdeg, tMWd, 2*tHd, tMgesd);
fprintf('\n')

% Mittlere Bewegungen berechnet

muS  = 2.95479E-04;  % [AE^3/Tage^2] für Sonne 
AE   = 149597870e3;  % AE in m 
aE = 1.0;
aM = 1.5237;

nE = sqrt(muS/aE^3); % in rad/Tag
nM = sqrt(muS/aM^3);

TsynodM = +2*pi/(nE-nM);

aHM = (aE+aM)/2;
tHM = pi*sqrt(aHM^3/muS);
psi0M = pi- nM*tHM;
psiFM = pi- nE*tHM;
tWM =  +(2*psiFM+2*pi)/(nE-nM);


ttlprint = "Flugzeiten zum Mars (berechnet)";
fprintf('\n T_synodM = %7.4f d',TsynodM);
fprintf('\n\n');
fprintf('\n %s', ttlprint);
fprintf('\n\n')
fprintf('|  aE (m)    | aM (m)     | nE (rad/d) | nM (rad/d) |   tH (d)   |\n');
fprintf('| %7.4e | %7.4e | %7.4e | %7.4e |  %7.4f  |\n',...
                    aE*AE, aM*AE, nE, nM, tHM);
fprintf('\n')
fprintf('|   psi0(°)  |  psi_F(°) |   tMW (d)   |  2*tH (d)  |  tMges (d) |\n');
fprintf('|  %7.3f   |  %7.3f  |   %7.3f   |  %7.3f   |  %+7.3f  |\n',...
                    rad2deg(psi0M), rad2deg(psiFM), tWM, 2*tHM, 2*tHM+tWM);
fprintf('\n')



%% Venus

G = 6.67430e-11; % m3 kg/s2
mS = 1.989e30; % kg
muS = G*mS;
aE = 149598022e3;% m
aV = 108160000e3;
aH = (aE+aV)/2;
tH = pi*sqrt(aH*aH*aH/muS);
tHd = tH/60/60/24;
nV = 58518.220;
nVSI = nV * pi/180 /100/365.25/60/60/24;
psi_0 = pi - nVSI*tH;
psi_0deg = psi_0* 180/pi;
psi_F = pi - nESI*tH;
psi_Fdeg = psi_F* 180/pi;
tVW = (2*psi_F - 2*pi*(0))/(nVSI-nESI);  % -1 ???
tVWd = tVW / 60/60/24;
tVgesd = tVWd + 2*tHd;

fprintf('\n ------------------------------------------------------------')
fprintf('\n')
ttlprint = "Flugzeiten zur Venus (Tabellenwerte)";
fprintf('\n %s', ttlprint);
fprintf('\n\n')
fprintf('|  aE (m)    | aV (m)     | nE (rad/d) | nV (rad/d) |   tH (d)   |\n');
fprintf('| %7.4e | %7.4e | %7.4e | %7.4e |  %7.4f  |\n',...
                    aE, aV, nESI*86400, nVSI*86400, tHd);
fprintf('\n')
fprintf('|   psi0(°)  |  psi_F(°) |   tVW (d)   |  2*tH (d)  |  tVges (d) |\n');
fprintf('|  %7.3f   |  %7.3f  |   %7.3f   |  %7.3f   |  %+7.3f  |\n',...
                    psi_0deg, psi_Fdeg, tVWd, 2*tHd, tVgesd);
fprintf('\n')

muS  = 2.95479E-04;  % [AE^3/Tage^2] für Sonne 
AE   = 149597870e3;  % AE in m 
aV = 0.7233;         % in AE
aE = 1.0;
nE = sqrt(muS/aE^3); % in rad/Tag
nV = sqrt(muS/aV^3);

TsynodV = +2*pi/(nV-nE);

aHV = (aE+aV)/2;
tHV = pi*sqrt(aHV^3/muS);
psi0V = pi- nV*tHV;
psiFV = pi- nE*tHV;
tWV =  +(2*psiFV+0*pi)/(nV-nE);


fprintf('\n T_synodV = %7.4f d',TsynodV);
fprintf('\n\n');
ttlprint = "Flugzeiten zur Venus (berechnet)";
fprintf('\n %s', ttlprint);
fprintf('\n\n')
fprintf('|  aE (m)    | aV (m)     | nE (rad/d) | nV (rad/d) |   tH (d)   |\n');
fprintf('| %7.4e | %7.4e | %7.4e | %7.4e |  %7.4f  |\n',...
                    aE*AE, aV*AE, nE, nV, tHV);
fprintf('\n')
fprintf('|   psi0(°)  |  psi_F(°) |   tVW (d)   |  2*tH (d)  |  tVges (d) |\n');
fprintf('|  %7.3f   |  %7.3f  |   %7.3f   |  %7.3f   |  %+7.3f  |\n',...
                    rad2deg(psi0V), rad2deg(psiFV), tWV, 2*tHV, 2*tHV+tWV);
fprintf('\n')


