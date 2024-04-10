% -------------------------------------------------------------------------
% PhaseF.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Phases zur Berechnung von Finsternisses
% Vergleicht die ekliptikalen Laengen von Sonne und Mond zum Zeitpunkt t
% -------------------------------------------------------------------------

function Phase = PhaseF(t,SDL,SDR)
   tau = 8.32/1440/36525;
   T = t*36525+2451545.0;   % Zeit in Julian. Datum 
   eps0 = deg2rad(EpsErde(T));
   MoonPos = MondExakt(t,eps0);
   SunPos = SonneExakt(t-tau,SDL,SDR,eps0);
   % Gibt den Unterschied der eklipzikalen Längen zw. Mond und Sonne wieder
   Phase = wrapTo360( MoonPos.ekl1(2)-SunPos.ekl1(2)+180)-180;
end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------