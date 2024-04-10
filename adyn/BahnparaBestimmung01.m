% -------------------------------------------------------------------------
% BahnparaBestimmung01.m
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
% aus einem bekannten Ort- und Geschwindigkeitsvektor
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

GME = 398600.4415;   % G*ME in km^3/s^2
RE  = 6378;          % Erdradius in km

%% Berechnung Bahnparameter (HEO) 

% Satellit #1 (HEO)

 rv  = [12000;32000;-6000];
 vv  = [-1.25;1.3;-0.2];
 
%  rv  = [10000;40000;-5000];
%  vv  = [-1.5;1.0;-0.1];

 r   = sqrt(dot(rv,rv));
 v   = sqrt(dot(vv,vv));

hv      = cross(rv,vv);                         % Drehimpuls
h       = sqrt(dot(hv,hv));
Wv      = hv/h;
i       = atand(sqrt(hv(1)^2+hv(2)^2)/hv(3));   % Inklination
Omega   = atan2d(Wv(1),-Wv(2));                 % Knoten
p       = h^2/GME;                              % semi-latus rectum
a       = 1/(2/r-v^2/GME);                      % große Halbachse
nBew    = sqrt(GME/a^3);                        % mittl. Bewegung, 3. KG
exz     = sqrt(1-p/a);                          % Exzentrizität
E       = atan2(dot(rv,vv)/(a^2*nBew),(1-r/a)); % Exzentr. Anomalie
M       = E-exz*sin(E);                         
M       = rad2deg(M);                           % Mittl. Anomalie
u       = atand(rv(3)/(-rv(1)*Wv(2)+rv(2)*Wv(1)));    % Länge
ups     = atand(sqrt(1-exz^2)*sin(E)/(cos(E)-exz));   % wahre Anomalie
omega   = u-ups;

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



%% Berechnung Bahnparameter (LEO) (fast kreisförmig, geringe Neigung)

% Satellit #2 (fast kreisförmig, geringe Neigung)
%  rv  = [6600;5800;-1000];
%  vv  = [-4.2;4.9;-0.2];

% ISS @ 2023-03-20-21 UTC 12:00.000 
 rv = [-4804.680472139880;-4787.203226753670;418.697383172504];
 vv = [3.084667418034520;-3.63637959795974;-5.99889734723948];

% ISS @ 2023-03-20 UTC 19:00:
%  rv = [4533.406242979590;5058.438123240110;282.9339789484634];
%  vv = [-3.70877528776069;2.97298546351788;6.00647691190360];

% ISS @ 2023-05-21 UTC 04:44:00.000 
%  rv = [2689.182355714180; -4337.832989289070;   4480.047872433510];
%  vv = [-3.84964001446308;    -5.767333716883710;  -3.25983631894374];

 r   = sqrt(dot(rv,rv));
 v   = sqrt(dot(vv,vv));

% Schritt 1
Wv  = cross(rv,vv);
Wv  = Wv/sqrt(dot(Wv,Wv));
Wx = Wv(1); Wy = Wv(2); Wz = Wv(3);

p  =  Wx/(1+Wz);
q  = -Wy/(1+Wz);

% Schritt 2

Av = cross(vv,cross(rv,vv))- GME * rv/r;

fv = [1-p^2+q^2;  2*p*q;     -2*p]/(1+p^2+q^2);
gv = [2*p*q;      1+p^2-q^2;  2*q]/(1+p^2+q^2);

k  = dot(Av,fv)/GME;
h  = dot(Av,gv)/GME;

% Schritt 3

X  = dot(rv,fv);
Y  = dot(rv,gv);

% Schritt 4

a    = 1/(2/r - v^2/GME);
nBew = sqrt(GME/a^3);
NN   = sqrt(1-h*h-k*k);
beta = 1/(1+NN);
cosL = k + ((1-k^2*beta)*X-h*k*beta*Y)/NN/a;
sinL = h + ((1-h^2*beta)*Y-h*k*beta*X)/NN/a;

Lambda = asin(sinL);
Lambda = acos(cosL);
lambda = Lambda - k*sinL + h*cosL;

% Schritt 5

Omega = atan2d(p,q);
omega = atan2d(h,k)- Omega;
exz   = sqrt(h^2+k^2);
i     = 2*atand(p/sind(Omega));
M     = wrapTo360(rad2deg(lambda)-Omega-omega);
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



hv      = cross(rv,vv);                         % Drehimpuls
h       = sqrt(dot(hv,hv));
Wv      = hv/h;
i       = atand(sqrt(hv(1)^2+hv(2)^2)/hv(3));   % Inklination
Omega   = atan2d(Wv(1),-Wv(2));                 % Knoten
p       = h^2/GME;                              % semi-latus rectum
a       = 1/(2/r-v^2/GME);                      % große Halbachse
nBew    = sqrt(GME/a^3);                        % mittl. Bewegung, 3. KG
exz     = sqrt(1-p/a);                          % Exzentrizität
E       = atan2(dot(rv,vv)/(a^2*nBew),(1-r/a)); % Exzentr. Anomalie
M       = E-exz*sin(E);                 
M       = rad2deg(M);                           % mittl. Anomalie
u       = atand(rv(3)/(-rv(1)*Wv(2)+rv(2)*Wv(1)));    % Länge
ups     = atand(sqrt(1-exz^2)*sin(E)/(cos(E)-exz));   % wahre Anomalie
% omega   = u-ups;

%Ausgabe

fprintf('\n');
fprintf('\n Satellit #2 (LEO) Standard-Berechnung');
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
