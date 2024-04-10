% -------------------------------------------------------------------------
% KeplerSonneEx.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Berechnet Rektaszension und Deklination der Sonne unter Verwendung
% einer analytischen Reihenentwicklung der Keplerloesung mit relativ
% geringer Genauigkeit. Im Unterschied zu KeplerSonne.m werden hier
% unterschiedliche Exzentrizitaeten beruecksichtigt. 
%
% Eingabe:
%   T         Zeit in Julianischem Datum
%   eps       Ekliptikschiefe
%   exz       Verhaeltnis der Exzentrizaet zu e = 0.016709
% Ausgabe:
%   RA        Rektaszension der Sonne in [rad]
%   Dec       Deklination der Sonne in [rad]
% -------------------------------------------------------------------------

function [RA,Dec] = KeplerSonneEx(T,eps,exz) 
  pi2=2*pi;
  t=Jd2JJht(T);
% Mittlere Anomalie und ekliptikale Laenge
  M  = pi2 * Frac(0.993133 + 99.997361*t); 
  L  = pi2 * Frac(0.7859453+M/pi2+(exz*(6892.72*sin(M)+...
       71.98*exz*sin(2*M)+1.04*exz*exz*sin(3*M))+6191.2*t)/(1296.0*1000));
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