% -------------------------------------------------------------------------
% PhaseFunc.m
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
% PhaseFunc(t,SunData, MPhase)
% Differenz zwischen dem Laengenabstand des Mondes von der Sonne
% und dem dazugehoerigen Wert fuer eine gegebene Phase. Vergleicht dazu 
% die ekliptikalen Laengen von Sonne und Mond zum Zeitpunkt t.
%
% Eingabe: 
%   t :          Zeit in Julianischen Jahrhunderten
%   SundData :   Parameter-Tables für Berechnung der Sonnenkoordinaten
%   MPhase:      Mondphase 1 = Neumond, 2 = Vollmond
% Ausgabe:
%   Phase:       Identifiziert die gesuchte Phase:   
% -------------------------------------------------------------------------

function Phase = PhaseFunc(t,SunData, MPhase)
% Phases zur Berechnung von Finsternissen
   tau = 8.32/1440/36525;  % Lichtlaufzet
   T = JJht2Jd(t);  %Zeit in Julian. Datum 
   epsE = deg2rad(EpsErde(T));
   MoonPos = MondExakt(t,epsE,'RE'); % Mondkoordinaten in 
   SunPos  = SonneExakt(t-tau,SunData,epsE,'RE'); % Sonnenkoordinaten
   % Gibt den Unterschied der ekliptikalen Laengen zw. Mond und Sonne wieder
   Phase = wrapTo360(rad2deg(MoonPos.ekl(2)-SunPos.ekl(2))- ...
                      MPhase*90+180)-180;
end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------