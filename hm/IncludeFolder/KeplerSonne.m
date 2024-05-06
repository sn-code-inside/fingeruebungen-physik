% -------------------------------------------------------------------------
% KeplerSonne.m
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
%
% von Oliver Montenbruck und Thomas Pfleger (Springer 1999). 
% Genehmigung des Verlages und der Autoren liegt vor.
% -------------------------------------------------------------------------
% Berechnet Rektaszension und Deklination der Sonne unter Verwendung
% einer analytischen Reihenentwicklung der Keplerloesung mit relativ
% geringer Genauigkeit. 
%
% Eingabe:
%   T         Zeit in Julianischem Datum
%   eps       Ekliptikschiefe
% Ausgabe:
%   RA        Rektaszension der Sonne in [rad]
%   Dec       Deklination der Sonne in [rad]
% -------------------------------------------------------------------------

function [RA,Dec]=KeplerSonne(T,eps) 
  pi2=2*pi;
  t=Jd2JJht(T);
% Mittlere Anomalie und ekliptikale Laenge
  M  = pi2 * Frac(0.993133 + 99.997361*t); 
  L  = pi2 * Frac(0.7859453+M/pi2+((6892.72*sin(M)+...
       71.98*sin(2*M)+1.04*sin(3*M))+6191.2*t)/(1296.0*1000));
% Aequatoriale Position
  ekl_Sun = CalcXYZfromAngles([ones(1,length(t)); L; zeros(1,length(t))]);
  ekl_Sun = mtimes(R_x(-eps),ekl_Sun); % Koordinatentrafo in Equ-System
  ekl_Sun = CalcAnglesfromXYZ(ekl_Sun);
  RA     = ekl_Sun(2,:);
  Dec    = ekl_Sun(3,:); 
end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------