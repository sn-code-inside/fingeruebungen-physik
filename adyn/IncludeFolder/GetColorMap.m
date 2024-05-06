% -------------------------------------------------------------------------
% GetColorMap.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Lädt Farbtabelle für zweidimensionale Gitter
% -------------------------------------------------------------------------

function cmap=GetColorMap
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
    tempcmap = [TempColors.R(2) TempColors.G(2) TempColors.B(2)
    TempColors.R(3) TempColors.G(3) TempColors.B(3)
    TempColors.R(4) TempColors.G(4) TempColors.B(4)];

    cmap = imresize(tempcmap, [1000, 3]);
    cmap = min(max(cmap, 0), 1);
end
% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------
