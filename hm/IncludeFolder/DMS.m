% -------------------------------------------------------------------------
% DMS.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Ermittelt Grad, Minuten und Sekunden zu in Dezimalform (°) gegebenem  
% Winkel 
% Eingabe:
%   DG        Grad in dezimaler Darstellung
% Ausgabe:    x[G,M,S]
%   G         Grad
%   M         Minuten
%   S         Sekunden
% -------------------------------------------------------------------------

function xDMS = DMS(DG)
% Ermitteln des Vorzeichens mit der Signum-Funktion 
  VZ = sign(DG);
% Umrechnung und Ausgabe 
  Sec  = 3600*abs(DG);
  xDMS = [VZ*floor(Sec/3600),floor(rem(Sec,3600)/60),...
          floor(rem(rem(Sec,3600),60))];
end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------

