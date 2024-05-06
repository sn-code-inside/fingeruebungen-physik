% -------------------------------------------------------------------------
% GMSTsat.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Finger√ºbungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Bestimmt die GMST in einer guten N‰herung
% -------------------------------------------------------------------------

function GMST=GMSTsat(Jd)
  d      = MJD(Jd)-51544.5;
  GMST   = wrapTo360(280.4606+360.9856473*d);
end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------
