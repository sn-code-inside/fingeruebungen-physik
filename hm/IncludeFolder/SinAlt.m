% -------------------------------------------------------------------------
% SinAlt.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------

function sinusalt=SinAlt(MJT, hour, lambda, cosphi, sinphi)
  MJD = 0;
  MJD = MJT+hour/24;
  T = MJD+2400000.5;
  eps = deg2rad(23.43929111);
  [xRA,xDec] = KeplerSonne(T,eps);
  tau = GMST(MJD)+deg2rad(lambda)-xRA;
  tau = wrapToPi(tau);
  sinusalt = sinphi*sin(xDec)+cosphi*cos(xDec)*cos(tau);
end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------
