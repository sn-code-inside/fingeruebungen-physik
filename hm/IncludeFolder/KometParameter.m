% -------------------------------------------------------------------------
% KometParameter.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Bahnparameter der Kometen und Asteroiden einlesen
% -------------------------------------------------------------------------

% Einlesen der Daten von Standish et al.
function Komet=KometParameter(fname)
    fid=fopen(fname,'r');
    if fid ==1 
        disp('File open not successful');
    else
        Komet = readtable(fname)
    end
    closeresult =fclose(fid);
    if closeresult ==0
    %    disp('Data successfully loaded');
    else
        disp('Dataload not successful');
    end
end
% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------
