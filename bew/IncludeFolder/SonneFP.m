% -------------------------------------------------------------------------
% SonneFP.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Graphische Ausgabe der Sonne und der Richtung zum Frühlingpunkt (FP).
% In  : ------
% Out : ------
% -------------------------------------------------------------------------
 
%_________________________________________________________________________
%%

function SonneFP
    Colors = GetColorLines;
    plot(0,0,'o', 'Color', Colors(10,:),'LineWidth',3);  %Sonne
    x=linspace(0,100,100);
    y=0*x;
    p2=plot(x,y,'Color','black','LineStyle', '--','LineWidth',1); %FP
end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------
