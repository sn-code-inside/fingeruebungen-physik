% -------------------------------------------------------------------------
% LMST.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% mittlere Ortssternzeit (local mean sidereal time)                   
% -------------------------------------------------------------------------

function outLMST = LMST(xMJD,LAMBDA)
  Secs = 86400.0;        % Anzahl der Sekunden je Sonnentag
  pi2 = 2*pi;
  xMJD_0 = floor (xMJD);
  UT    = 24*(xMJD-xMJD_0);     % [sec]
  T_0   = (xMJD_0-51544.5)/36525.0; 
  GMST  = 6.697374558 + 1.0027379093*UT +(8640184.812866+...
         (0.093104-6.2e-6*T_0).*T_0).*T_0/3600.0;
  outLMST = 24.0*Frac((GMST+LAMBDA/15.0)/24.0);
end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------