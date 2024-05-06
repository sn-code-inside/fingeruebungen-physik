% -------------------------------------------------------------------------
% PlanetoCentricCoords.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Derivation of planetocentric longitude of perihelion and right ascensnion
% of fictitious mean sun using geocentric equatorial position of Pole Vector
% -------------------------------------------------------------------------

clc
clear all
close all 
Colors=GetColorLines;

% Initialisierung 
dt1 = datetime('2000-01-02 12:00:00');
% Julianische Startdatum in UT
T1      = juliandate(dt1);
epsE    = EpsErde(T1);  
epsErad = deg2rad(epsE);  

Aequi = 'Datum';
%Aequi = 'J2000';
[BaPa,BaPadot]=OrbitParameter(T1,Aequi);

fprintf('\n');

for PlNr=1:8
    %orbit parameters
    iP     = BaPa.iP(PlNr);      %orbital inclination mars
    OmegaP = BaPa.OmegaP(PlNr);  %longitude of node
    PiqP   = BaPa.PiqP(PlNr);    %longitude of perihelion
    omega  = PiqP-OmegaP;        %angle of perihelion
    Omega90= OmegaP-90;
    M0     = BaPa.M0P(PlNr);     %mean anomalie at epoch

   [alpha0,delta0,Sense] = RotAxis(PlNr,Jd2JJht(T1));
    %convert to geocentric ecliptical coords
    axis_b = asind(-cosd(delta0)*sind(alpha0)*sind(epsE)+ ...
             sind(delta0)*cosd(epsE));
    axis_l = wrapTo360(atan2d(cosd(delta0)*sind(alpha0)*cosd(epsE)+ ... 
             sind(delta0)*sind(epsE),cosd(delta0)*cosd(alpha0))); 
    %rotation matrices
    Rz_O   = [cosd(OmegaP) -sind(OmegaP) 0; sind(OmegaP) cosd(OmegaP) 0; 0 0 1];
    Rz_mO  = [cosd(-OmegaP) -sind(-OmegaP) 0; sind(-OmegaP) cosd(-OmegaP) 0;0 0 1];
    Rz_O90 = [cosd(Omega90) -sind(Omega90) 0;...
             sind(Omega90) cosd(Omega90) 0;0 0 1];
    Rz_mO90= [cosd(-Omega90) -sind(-Omega90) 0;...
             sind(-Omega90) cosd(-Omega90) 0;0 0 1];
    Rz_o   = [cosd(omega) -sind(omega) 0; sind(omega) cosd(omega) 0;0 0 1];
    Rz_mo  = [cosd(-omega) -sind(-omega) 0; sind(-omega) cosd(-omega) 0;0 0 1];
    Rx_i   = [1 0 0;0 cosd(iP) -sind(iP); 0 sind(iP) cosd(iP)];
    Rx_mi  = [1 0 0;0 cosd(-iP) -sind(-iP); 0 sind(-iP) cosd(-iP)];
    
    %convert polar axis to cart. geocentric coords
    p      = [0,0,0];
    [p(1),p(2),p(3)] = sph2cart(axis_l/180*pi,axis_b/180*pi,1);
    
    %convert polar axis cart. geocentric coords to planeto-orbital cart. coords
    %can be achieved by two subsequent rotations of about z-axis with -omega 
    %and about x-axis with inclination -i
    p_bahn = Rx_mi*Rz_mO*p';
    
    %calculate angular momentum vector of orbit in orbital cart. coordinates
    L_bahn = [0,0,1];
    
    %calculate vernal point in orbit coordinates
    %vernal point is perpendicular to angular momentum and pependicular to
    %polar axis as it is the crossing point of equator and orbit
    vernal = cross(p_bahn',L_bahn);
  
    %convert vernal point orbital cart. coords to spherical coordinates
    [a,b,c]=cart2sph(vernal(1),vernal(2),vernal(3));
      
    %determine planeto-centric longitude of perihelion
    %by LSP = omega - a + 180° (180° because of reversal sun-planet coords
    LSP = wrapTo360(omega-rad2deg(a)+180);
    %determine right ascension of fictitious mean sun by
    % alpha_FMS = LSP + M0
    alpha_FMS = LSP+M0;
    fprintf('%s:   \t \t LSP: % 9.3f°\t AlphaFMS: % 9.3f° \n', ...
            string(BaPa.Name(PlNr)), LSP , wrapTo360(alpha_FMS));
end
fprintf('\n');


%%
% Functions
%
%--------------------------------------------------------------------------
% RotAxis: Berechnet die Orientierung des planetozentrischen Bezugssystems+
%         liefert die Rotationsrichtung
% Eingabe:
%   PlNr      Identifiziert den Planeten 
%   t         Zeit in Julianischen Jahrhunderten seit J2000
% Ausgabe:
%   alphaA,
%   deltaA    Richtunsgvekroren der Polachse des Planeten
%   Sense     Rotationsrichtung (direkt oder retrograd)
%--------------------------------------------------------------------------
function [alphaA, deltaA, Sense] = RotAxis (PlNr,t)
  % Rektaszension und Deklination der Rotationsachse bezogen auf
  % Aequator und Ekliptik zur Epoche J2000; Orientierung des
  % Nullmeridians
  % Position der Planetenachse in geozentrisch äquatorialen Koord. Ref. [3]

  switch (PlNr)
    case 1
        RA  = 281.01  -   0.033*t;  Dec =  61.45  -   0.005*t;    
    case 2
        RA  = 272.76; Dec =  67.16; 
    case 3
        RA  =   0.00  -   0.641*t;  Dec =  90.00  -   0.557*t;
    case 4
        RA  = 317.681 -   0.108*t; Dec =  52.886 -   0.061*t;
    case 5
        RA  = 268.05  -   0.009*t; Dec =  64.49  +   0.003*t;
    case 6
        RA  =  40.589 -   0.036*t; Dec =  83.537 -   0.004*t;  
    case 7
        RA  = 257.311; Dec = -15.175;
    case 8
        N   = deg2rad(357.85+52.316*t);    
        RA  = 299.36  + 0.70*sin(N); 
        Dec =  43.46  - 0.51*cos(N); 
    case 9
        RA  = 313.02; Dec =   9.09;
  end
  alphaA  = RA; 
  deltaA  = Dec; 
  % Rotationsrichtung
  if ( PlNr==2 || PlNr==7 || PlNr==10 )
    Sense = +1;
  else
    Sense = -1;
  end
end 

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------
