% -------------------------------------------------------------------------
% OrbitParameter.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Einlesen der Bahnparameter der Planeten (Tabelle 3.8. und 3.9).
% In  : T           .. Float        ..  Julianisches Datum
%       Aequi       .. String       ..  Gewähltes Aequinoktium 'J2000' oder
%                                       'Date'
% Out : BaPa        .. Table        ..  Bahnparameter (Tabelle)
%       BaPadot     .. Table        ..  Bahnparameter Ableitung (Tabelle)
% -------------------------------------------------------------------------

function [BaPa,BaPadot]=OrbitParameter(T1,Aequi)
   if strcmp(Aequi,'Date')
      [BaPa,BaPadot]=OrbitParameterEinlesen(Aequi);
   else
      [BaPa,BaPadot]=OrbitParameterEinlesen(Aequi);
   end
end

function [BaPa,BaPadot]=OrbitParameterEinlesen(Aequi)
%-------------------------------------------------------------------------
% Einlesen der Bahnparameter der Planeten (Tabelle 3.8. und 3.9).
% In:   Aequi       .. String       ..  Gewähltes Aequinoktium 
% Out : BaPa        .. Table        ..  Bahnparameter (Tabelle)
%       BaPadot     .. Table        ..  Bahnparameter Ableitung (Tabelle)
%-------------------------------------------------------------------------
    if strcmp(Aequi,'Datum') 
        fname = 'Bahnparameter.csv';
    else
        fname = 'BahnparameterJ2000.csv';
    end
    fid=fopen(fname,'r');
    if fid ==1 
        disp('File open not successful');
    else
        BaPa = readtable(fname);
    end
    closeresult =fclose(fid);
    if closeresult ==0
    %     disp('Data successfully loaded');
    else
        disp('Dataload not successful');
    end
    if strcmp(Aequi,'Datum') 
        fname = 'BahnparameterAbleitung.csv';
    else
        fname = 'BahnparameterJ2000Ableitung.csv';
    end
    fid=fopen(fname,'r');
    if fid ==1 
    %     disp('File open not successful')
    else
        BaPadot = readtable(fname);
    end
    closeresult =fclose(fid);
    if closeresult ==0
    %     disp('Data successfully loaded');
    else
        disp('Dataload not successful');
    end
end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------
