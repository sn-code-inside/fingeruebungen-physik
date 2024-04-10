% -------------------------------------------------------------------------
% Site.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% SITE:  berechnet geozentrische aus geographischen Beobachterkoordinaten     
%        Site.c:  r * cos(phi') (geozentrisch; in Erdradien)                   
%        Site.s:  r * sin(phi') (geozentrisch; in Erdradien)                   
% Input: PHI:    geographische Breite (in Grad)                               
% -------------------------------------------------------------------------

function SITE = Site(PHI)
    ex = 0.006694;        % ex^2=f(2-f) mit Erdabplattung f=1/298.257 
    sinphi = sind(PHI);   
    cosphi = cosd(PHI);
    N = 1.0/sqrt(1.0-ex*sinphi*sinphi);
    SITE.c = N*cosphi;  
    SITE.s = (1.0-ex)*N*sinphi;
end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------
