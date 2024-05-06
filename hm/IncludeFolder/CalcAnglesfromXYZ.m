% -------------------------------------------------------------------------
%  CalcAnglesfromXYZ.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Berechnet die sphaerischen Koordinaten aus kartesischen Koordinaten 
% -------------------------------------------------------------------------

function VecOut=CalcAnglesfromXYZ(VecIn)
% Eingabe: 
%  VecIn :  kartesische Koordinaten als Vektor [X Y Z]
% Ausgabe: 
%  VecOut :  sphaerische Koordinaten als Vektor [r theta phi]
  VecOut=[  sqrt(VecIn(1,:).^2+VecIn(2,:).^2+VecIn(3,:).^2);
            atan2(VecIn(2,:),VecIn(1,:));
            atan2(VecIn(3,:),sqrt(VecIn(1,:).^2+VecIn(2,:).^2))];
end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------
