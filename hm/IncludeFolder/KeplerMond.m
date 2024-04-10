% -------------------------------------------------------------------------
% KeplerMond.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Berechnet Rektaszension und Deklination des Mondes mit relativ
% geringer Genauigkeit 

% Eingabe:
%   T        Zeit in Julianischen Jahrenn Jahrhunderten 
%   eps      Polneigung in Bogenmass
% Ausgabe:
%  RA        Rektaszension der Sonne in [rad]
%  Dec       Deklination der Sonne in [rad]
% -------------------------------------------------------------------------

function Moon = KeplerMond(T,epsErad) 
  pi2=2*pi;
  T=(T-2451545.0)/36525;
  ARC=206264.8062;
  % COSEPS=0.91748; SINEPS=0.39778; % cos/sin(Ekliptikschiefe) *)
  % mittlere Elemente der Mondbahn 
  L0 =  Frac(0.606433+1336.855225*T); % mittl. Laenge des Mondes (in rad)
  L  =  pi2*Frac(0.374897+1325.552410*T); % mittl. Anomalie des Mondes 
  LS =  pi2*Frac(0.993133+  99.997361*T); % mittl. Anomalie Sonne 
  D  =  pi2*Frac(0.827361+1236.853086*T); % Diff. Laenge Mond-Sonne
  F  =  pi2*Frac(0.259086+1342.227825*T); % Knotenabstand 
  DL =  22640*sin(L) - 4586*sin(L-2*D) + 2370*sin(2*D) + 769*sin(2*L)... 
          - 668*sin(LS) - 412*sin(2*F) - 212*sin(2*L-2*D) ...
          - 206*sin(L+LS-2*D) + 192*sin(L+2*D) - 165*sin(LS-2*D) ...
          - 125*sin(D) - 110*sin(L+LS) +148*sin(L-LS) - 55*sin(2*F-2*D);
  S = F + (DL+412*sin(2*F)+541*sin(LS))/ARC; 
  H = F-2*D;
  N = -526*sin(H) + 44*sin(L+H) - 31*sin(-L+H) - 23*sin(LS+H) ...
         + 11*sin(-LS+H) -25*sin(-2*L+F) + 21*sin(-L+F);
  % Geozentrisch-ekliptikale Koordinaten
  L_MOON = pi2 * Frac ( L0 + DL/(1296*1000)); % in rad Länge
  B_MOON = ( 18520.0*sin(S) + N ) / ARC;      % in rad Breite
  D_MOON = 385000-20905*cos(L)-570*cos(2*L)-3699*cos(2*D-L)-2956*cos(2*D);
                                              % in km Abstand
  D_MOON = 1*D_MOON;                                            
  % + eventuell weitere Terme
  % Geozentrisch-aequatoriale Koordinaten 
  Moon.ekl = [D_MOON; L_MOON; B_MOON];
  Moon.ekl_xyz = CalcXYZfromAngles(Moon.ekl);
  % Geozentrisch-aequatoriale Koordinaten 
  % Koordinatentrafo in Equ-System
  Moon.equ_xyz = mtimes(R_x(-epsErad),Moon.ekl_xyz); 
  Moon.equ = CalcAnglesfromXYZ(Moon.equ_xyz);
end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------
