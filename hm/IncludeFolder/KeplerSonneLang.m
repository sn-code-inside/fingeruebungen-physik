% -------------------------------------------------------------------------
% KeplerSonneLang.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Berechnet Rektaszension und Deklination der Sonne unter Verwendung
% einer analytischen Reihenentwicklung der Keplerloesung mit relativ
% geringer Genauigkeit. Im Unterschied zu KeplerSonne.m werden hier
% verschiedene Exzentrizitaeten beruecksichtigt und es wird ein langer
% Beobachtungszeitraum gewaehlt. 
%
% Eingabe:
%   TV       Zeit in Julianischen Datum
%   eps      Ekliptikschiefe
%   exz      Verhaeltnis der Exzentrizaet zu e = 0.016709
% Ausgabe:
%   L0       Laenge der mittleren Sonne in [rad]
%   RAV      Rektaszension der Sonne in [rad]
%   DecV     Deklination der Sonne in [rad]
% -------------------------------------------------------------------------

% KeplerSonne - Variable Exzentrizität:
function [L0, RAV,DecV]=KeplerSonneLang(T,eps,exz) 
%--------------------------------------------------------------------------
  pi2=2*pi;
  t=Jd2JJht(T);
% Mittlere Anomalie 
  M   = pi2*Frac((357.52772+35999.050340*t+0.0001603*t.*t-...
        t.*t.*t/300000)/360); 
% Laenge der mittleren Sonne
  L0  = pi2*Frac(280.46664567/360 +(36000.76982779*t+...
        0.0003032028*t.*t+t.*t.*t/4993000-t.*t.*t.*t/153000000-...
        t.*t.*t.*t.*t/2e11)/360);
% Mittelpunktsgleichung
  C0  = (2*exz-(exz^3)/4)*sin(M)+(5*exz*exz/4-11/24*(exz^4))*sin(2*M)+...
        13*exz^3*sin(3*M)/12;
% Ekliptikale Laenge der wahren Sonne
  L = L0+C0;
  ekl_Sun = CalcXYZfromAngles([ones(1,length(t));L;zeros(1,length(t))]);
  ekl_Sun = mtimes(R_x(-eps),ekl_Sun); % Koordinatentrafo in Equ-System
  ekl_Sun = CalcAnglesfromXYZ(ekl_Sun);
  RAV  = ekl_Sun(2,:);
  DecV = ekl_Sun(3,:); 
end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------