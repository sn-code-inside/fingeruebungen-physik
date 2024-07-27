% -------------------------------------------------------------------------
% GMST.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Berechnet die mittlere Greenwich-Sternzeit aus der Universal Time (UT)
% -------------------------------------------------------------------------

function outGMST = GMST(xMJD)
% Eingabe:
%   MJD Zeit (UT) als Modifiziertes Julianisches Datum
% Ausgabe:
%   Mittlere Greenwich-Sternzeit in [rad]
  Secs = 86400.0;        % Anzahl der Sekunden je Sonnentag
  pi2 = 2*pi;
% Berechnung der bereits vollstaendig verstrichenen Tage durch Abrunden 
  xMJD_0 = floor (xMJD);
% Berechnung der Universal Time UT 
  UT    = Secs*(xMJD-xMJD_0);     % [sec]
% Umrechnung in Greenwich Mean Siderial Time GMST 
  T_0   = (xMJD_0-51544.5)/36525.0; 
  T     = (xMJD  -51544.5)/36525.0; 
  sec_outGMST  = 24110.54841+8640184.812866*T_0+1.0027379093*UT+...
                 (0.093104-0.0000562*T).*T.*T;          % [sec]
  outGMST      = (pi2/Secs)*Modulo(sec_outGMST,Secs);   % [Rad]
end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------
