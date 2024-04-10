% -------------------------------------------------------------------------
% JJht2Jd.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Bestimmt das Julianische Datum aus dem Modifiziertem JuDatum
% -------------------------------------------------------------------------

function yout=JJht2Jd(t)
  yout=t*36525+2451545.0;
end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------
