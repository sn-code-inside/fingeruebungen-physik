% -------------------------------------------------------------------------
% PlotTrebuchet.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Hilfsfunktion zum Zeichnen Trebuchet
% -------------------------------------------------------------------------

function PlotTrebuchet(h,base,groundl,groundr,Colors,Style)
    line([groundl -base],[-h -h],'Color', Colors(8,:),...
         'LineStyle',Style(2), 'LineWidth',2);
    line([-base groundr],[-h -h],'Color', Colors(8,:),...
         'LineStyle',Style(2), 'LineWidth',2); 
    line([-base 0],[-h 0],'Color', Colors(8,:),...
         'LineStyle',Style(2), 'LineWidth',2); 
    line([+base 0],[-h 0],'Color', Colors(8,:),...
         'LineStyle',Style(2), 'LineWidth',2); 
    line([+base -base],[-h -h],'Color', Colors(8,:),...
         'LineStyle',Style(2), 'LineWidth',2); 
end
% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------
