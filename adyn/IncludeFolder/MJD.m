% -------------------------------------------------------------------------
% MJD.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Bestimmt das Modifizierte Julianische Datum aus dem Julianischem Datum
% -------------------------------------------------------------------------

function yout=MJd(Jd)
  yout=Jd-2400000.5;
end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------
