% -------------------------------------------------------------------------
% PerturbImport.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
%  JUPITERPOS1 = IMPORTFILE(FILENAME) reads data from text file FILENAME
%  for the default selection.  Returns the numeric data.
%
%  JUPITERPOS1 = IMPORTFILE(FILE, DATALINES) reads data for the
%  specified row interval(s) of text file FILENAME. Specify DATALINES as
%  a positive scalar integer or a N-by-2 array of positive scalar
%  integers for dis-contiguous row intervals.
%
%  Example:
%  JupiterPos1 = ImportFile02('JupiterPos.csv', [3, Inf]);
%
%  See also READTABLE.
% -------------------------------------------------------------------------

function Jup = PerturbImport(filename, dataLines)

% Input handling
% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [3, Inf];
end

% Setup the Import Options
opts = delimitedTextImportOptions("NumVariables", 10);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ";";

% Specify column names and types
opts.VariableNames = ["JupiterPos", "VarName2", "VarName3", "VarName4", "VarName5", "VarName6", "VarName7", "VarName8", "VarName9", "VarName10"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
Temp = readtable(filename, opts);

% Convert to output type
Jup = table2array(Temp);
    
end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------