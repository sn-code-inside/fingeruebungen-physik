% -------------------------------------------------------------------------
% SonneFP.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Graphische Ausgabe der Sonne und der Richtung zum Frühlingpunkt (FP).
% -------------------------------------------------------------------------


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