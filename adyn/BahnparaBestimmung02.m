% -------------------------------------------------------------------------
% BahnparaBestimmung02.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Astrodynamik" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% 
% Programm berechnet die Bahnparameter
% von  Satelliten um die Erde 
% aus zwei Beobachtungen zu unterschiedlichen Zeiten
%
% -------------------------------------------------------------------------

% Initialisierung
clc
clear 
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors=GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

muE = 398600.4415;   % G*ME in km^3/s^2
RE  = 6378;          % Erdradius in km

%% Berechnung Bahnparameter (HEO) 

% Beobachtungsstation G�rlitz

phiB    = 51;
lambdaB = 15;
vecrB   = RE*[cosd(lambdaB)*cosd(phiB); sind(lambdaB)*cosd(phiB);...
             sind(phiB)];
R1      = 2350;
R2      = 1732;
h1      = 15;
h2      = 29;
Az1     = -60;
Az2     = 33;

eO = [-sind(lambdaB);cosd(lambdaB);0];
eN = [-sind(phiB).*cosd(lambdaB);...
      -sind(phiB).*sind(lambdaB);...
      cosd(phiB)];
eZ = [cosd(phiB).*cosd(lambdaB);....
      cosd(phiB).*sind(lambdaB);...
      sind(phiB)];

U  = [eO,eN,eZ]';
Uinv = inv(U);

dt1    = datetime('2020-01-01 01:14:00');
T1     = juliandate(dt1); % Julianisches Datum  ET
theta1 = GMSTsat(T1);
dt2    = datetime('2020-01-01 01:20:00');
T2     = juliandate(dt2); % Julianisches Datum  ET
theta2 = GMSTsat(T2);

vecR1  = R1*[cosd(h1)*cosd(Az1); cosd(h1)*sind(Az1); sind(h1)]; 
vecR2  = R2*[cosd(h2)*cosd(Az2); cosd(h2)*sind(Az2); sind(h2)]; 

vecr1  = mtimes(R_z(-deg2rad(theta1)),mtimes(Uinv,vecR1)+vecrB);
vecr2  = mtimes(R_z(-deg2rad(theta2)),mtimes(Uinv,vecR2)+vecrB);
vece1  = vecr1/vecnorm(vecr1);
vecr0  = vecr2-dot(vecr2,vece1)*vece1;
r1     = vecnorm(vecr1);
r2     = vecnorm(vecr2);
r0     = vecnorm(vecr0);

W      = cross(vece1,vecr0);

inc    = atan2d(sqrt(W(1)^2+W(2)^2),W(3));
Omega  = atan2d(W(1),(-W(2)));
u1     = atan2d(vecr1(3),-vecr1(1)*W(2)+vecr1(2)*W(1));

tau    = sqrt(muE)*(T2-T1)*86400;
m      = tau^2/(sqrt(2*(r1*r2+dot(vecr1,vecr2))))^3;
k      = (r1+r2)/2/sqrt(2*(r1*r2+dot(vecr1,vecr2)))-1/2;

eta0   = 12/22 +10/22*sqrt(1+(44/9)*(m/(k+5/6)));
eta1   = eta0 + 0.1;
eta2   = eta1 - func(eta1,m,k)*(eta1-eta0)/(func(eta1,m,k)-func(eta0,m,k));
AD     = 0.5*r1*r0;
p      = (2*AD*eta2/tau)^2;

%Ausgabe

fprintf('\n');
fprintf('\n Satellit #1 (HEO)');
fprintf('\n');
fprintf('\n Halbachse         a : %8.2f km',a);
fprintf('\n Exzentrizität     e : %8.6f ',exz);
fprintf('\n Bahnneigung       i : %8.3f °',i);
fprintf('\n Knoten        Omega : %8.3f °',Omega);
fprintf('\n Arg.Perigäum  omega : %8.3f °',omega);
fprintf('\n Mittl. Anomalie   M : %8.3f °',M);
fprintf('\n Mittl. Bewegung   n : %8.3f Umrundungen/d',nBew*86400/2/pi);
fprintf('\n');
fprintf('\n');


%Ausgabe

fprintf('\n');
fprintf('\n Satellit #2 (LEO) (exakte) Berechnung');
fprintf('\n');
fprintf('\n Halbachse         a : %8.2f km',a);
fprintf('\n Bahnhöhe          h : %8.2f km',a-RE);
fprintf('\n Exzentrizität     e : %8.6f ',exz);
fprintf('\n Bahnneigung       i : %8.3f °',i);
fprintf('\n Knoten        Omega : %8.3f °',Omega);
fprintf('\n Arg.Perigäum  omega : %8.3f °',omega);
fprintf('\n Mittl. Anomalie   M : %8.3f °',M);
fprintf('\n Mittl. Bewegung   n : %8.3f Umrundungen/d',nBew*86400/2/pi);
fprintf('\n');
fprintf('\n');



% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------

function yout=R_z(phi)
% Eingabe:  
%   phi     = Rotationswinkel im Bogenmaß [rad]
% Ausgabe:
%   yout    = Rotationsmatrix R_z 
   yout=[cos(phi) sin(phi) 0 ;...
        -sin(phi) cos(phi) 0 ;...
        0 0 1];
end

function yout=func(etain, m, k)
% Eingabe:  
%   eta_in    
% Ausgabe:
%   yout    
   yout= 1 -etain + (m/etain^2)*(4/3 + (24/15)*(m/etain^2 - k) + ...
                             (24*8/105)*(m/etain^2 - k)^2);
end

% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------
