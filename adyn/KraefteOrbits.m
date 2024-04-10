% -------------------------------------------------------------------------
% KraefteOrbits.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Astrodynamik" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% 
% Programm berechnet die Größenordnungen verschiedener externer Kräfte im
% erdnahen Orbits.
%
% -------------------------------------------------------------------------

% Initialisierung
clc
clear 
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style  = ["-", "-.", ":", "--", ":"];

% Parameter, hier alles in m, kg, s
G   = 6.671e-11;      % G in in m^3/s^2/kg
AE  = 149597870700;   % m
aV  = 0.72333199*AE;  % Bahnradius Venus in m
aV  = 0.723*AE;       % Bahnradius Venus in m
MS  = 1.989e30;       % Masse Sonne in kg 
MM  = 4.867e24;       % Masse Mond in kg
rM  = 384000000;

% Planetendaten, ab hier alles in km,s,kg
GS   = MS*G*1e-9;   % Sonne
RV   = 6051;        % Planetenradius Venus in km
R0   = aV/1000;     % Radius SOI in km
v0   = sqrt(GS/R0);         % Umlaufgeschwindigkeit Venus @ SOI

f      = 1/298.257;  % Abplattung
RE     = 6378;       % in km
GE     = 398600.4;   % in km^3 s^-2
omegaE = 0.7292e-4;  % in s^-1 

%% Berechnung zweite zonale Harmonische
r1 = RE+400;
r2 = RE+2000;
r  = [r1 r2];
aG = (RE/r1)^2-1;

J2     = 2*(f - omegaE^2*RE^2/(2*GE/RE))/3;
aJ2    = 2*GE*RE^2*J2./r.^4;

aJ4 = -1.3e-6;

aS = G*MS*RE*1000/AE^3;
aM = G*MM*RE*1000/rM^3;



%% maximal möglicher Streuwinkel und Winkel zw. vinf und vP
%  (alles in km /s)
vF   = sqrt(2*GV./(RV+500)); % Perizentrums-Fluchtgeschwindigkeit in km/s
% Winkel zwischen vinf und vP
thetad = 2*asind(1./(1+2*vinf.^2./vF.^2));
totalthetad = sum(thetad);

%% Berechnung realer Streuwinkel und Winkel zw. vinf und vP
%  (alles in km /s)
RVreal = [1.67 1.49 1.15 1.43 1.70 1.79 1.37 100000]*RV;
vFreal   = sqrt(2*GV./RVreal); % Perizentrums-Fluchtgeschwindigkeit in km/s
% Winkel zwischen vinf und vP
thetad_real = 2*asind(1./(1+2*vinf.^2./vFreal.^2));
totalthetad_real = sum(thetad_real);

%% Ausgabe Werte

nrorbits = ["a" "b" "c" "d" "e" "f" "g" "h"];
fprintf('\n Orbit |  r_P  (AE)  | r_A  (AE)  |  beta_1 (°) |    e     |   T_p (Tage) | v_p (km/s)\n');
for k=1:length(vP)
    fprintf('  %s    |   %5.3f     |   %5.3f    |   %5.2f     |  %5.3f   |   %6.1f     |  %6.2f \n',...
            nrorbits(k), rP(k)*1000/AE, rA(k)*1000/AE, beta1d(k), ecc(k), Tpd(k), vP(k));
end
fprintf('\n');
fprintf('\n GA-Manöver  |  theta (°)  |  theta_total (°)  | R_V,real/R_V ');
for k=1:length(vP)-1
    fprintf('\n %u  %s -> %s   |    %5.3f    |     %6.3f        |  %5.2f  ',...
            k, nrorbits(k), nrorbits(k+1), thetad_real(k), sum(thetad_real(1:k)), RVreal(k)/RV);
end
fprintf('\n');




%% ------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


%% Funktionen
%Ende Funktionen
% -------------------------------------------------------------------------

