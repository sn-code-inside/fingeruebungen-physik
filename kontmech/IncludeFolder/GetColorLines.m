% -------------------------------------------------------------------------
% GetColorLines.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Lädt Farbtabelle für Linien
% -------------------------------------------------------------------------

function Colors=GetColorLines
% Ausgabe:
% Farbtabelle
    fid=fopen('PlanetenFarben.csv','r');
    if fid ==1 
        disp('File open not successful');
    else
        TempColors = readtable('PlanetenFarben.csv');
    closeresult =fclose(fid);
    if closeresult ==0
    %    disp('Data successfully loaded');
    else
        disp('Dataload not successful');
    end
    Colors=zeros(15,3);
    for m=1:15 
        Colors(m,1) = TempColors.R(m);
        Colors(m,2) = TempColors.G(m);
        Colors(m,3) = TempColors.B(m);
    end
end
% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------
