% -------------------------------------------------------------------------
% ImportPlanetParameter.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Einlesen typischer physikalischer Bahnparameter der Planeten 
% und des Mondes.
%
%  Import data from a text file
%  PlanetenParameter = ImportPlanetParameter(FILENAME) reads data from 
%  text file
%  FILENAME for the default selection.  Returns the data as a table.
%
%  PlanetenParameter = ImportPlanetParameter(FILE, DATALINES) 
%  reads data for the
%  specified row interval(s) of text file FILENAME. 
%  Specify DATALINES as
%  a positive scalar integer or a N-by-2 array of positive scalar
%  integers for dis-contiguous row intervals.
% -------------------------------------------------------------------------

function PlanetenParameter = ImportPlanetParameter(filename, dataLines)
%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [2, Inf];
end

%% Setup the Import Options
opts = delimitedTextImportOptions("NumVariables", 5);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ";";

% Specify column names and types
opts.VariableNames = ["Name", "Masse", "aP", "rhoP", "RP"];
opts.VariableTypes = ["string", "double", "double", "double", ...
                      "double"];
opts = setvaropts(opts, 1, "WhitespaceRule", "preserve");
opts = setvaropts(opts, 1, "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
PlanetenParameter = readtable(filename, opts);

end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------