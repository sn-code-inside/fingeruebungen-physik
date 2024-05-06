% -------------------------------------------------------------------------
% KeplerMondM.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
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

function Moon = KeplerMondM(T,epsErad) 
%--------------------------------------------------------------------------
  pi2=2*pi;
  T=(T-2451545.0)/36525;
  ARC=206264.8062;
  % mittlere Elemente der Mondbahn 
  L0 =  Frac(0.606433+1336.855225*T); % mittl. Laenge des Mondes (in rad)
  L  =  pi2*Frac(0.374897+1325.552410*T); % mittl. Anomalie des Mondes 
  % Geozentrisch-ekliptikale Koordinaten
  % Geozentrisch-aequatoriale Koordinaten 
  B_MOON = 0;
  D_MOON = 385000;
  L_MOON =  pi2 * Frac (L0);
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
