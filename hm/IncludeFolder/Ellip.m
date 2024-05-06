% -------------------------------------------------------------------------
% Ellip.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Berechnet Ort und Geschwindigkeit fuer elliptische Bahnen
% -------------------------------------------------------------------------

function Planets=Ellip(aP,eP,M,mu,k)
% Eingabe: 
%   ap = Laenge der großen Halbachse 
%   eP = Exzentrizität des Planeten 
%   M  = Mittlere Anomalie 
%   mu = Gravitationskonstante *(Masse_Zentralkörper + Masse_Planet) 
%   k  = Nummer des Planeten 
% Ausgabe: 
% 
   ka=sqrt(mu/aP);          % vgl. Formel 3.38 
   Ecc=EAnom(M,eP);         % exzentrische Anomalie 
   cosE=cos(Ecc);
   sinE=sin(Ecc);
   fac=sqrt((1-eP)*(1+eP));
   rho=1-eP*cosE;
   % x-Abstand zwischen Planet und Zentralkoerper 
   x=aP*(cosE-eP); 
   % y-Abstand zwischen Planet und Zentralkoerper 
   y=aP*fac*sinE; 
   nu=2*atan(sqrt((1+eP)/(1-eP))*tan(Ecc/2));
   r=aP*rho;
   % Position des Planeten
   Planet(k).Pos=[x;y];
   % Orbit des Planeten 
   Planet(k).Orb=[r;nu];
   % Geschwindigkeit des Planeten 
   Planet(k).Vel=[-ka*sinE/rho; ka*fac*cosE/rho];
   Planets=Planet(k);
end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------
