% -------------------------------------------------------------------------
% ImportFile02.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Import data from a text file
%  DK01LGL = IMPORTFILE(FILENAME) reads data from text file FILENAME for
%  the default selection.  Returns the data as a table.
%
%  DK01LGL = IMPORTFILE(FILE, DATALINES) reads data for the specified
%  row interval(s) of text file FILENAME. Specify DATALINES as a
%  positive scalar integer or a N-by-2 array of positive scalar integers
%  for dis-contiguous row intervals.
%
%  Example:
%  DK01LGL = ImportFile02('DK01yzLGL.txt', [2, 94]);
%
%  See also READTABLE.
% -------------------------------------------------------------------------

function datafile = ImportFile02(filename, dataLines)

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [2, Inf];
end

%% Setup the Import Options
opts = delimitedTextImportOptions("NumVariables", 6);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = " ";

% Specify column names and types
opts.VariableNames = ["q", "dq", "Var3", "Var4", "Var5", "Var6"];
opts.SelectedVariableNames = ["q", "dq"];
opts.VariableTypes = ["double", "double", "string", "string", "string", "string"];
opts = setvaropts(opts, [3, 4, 5, 6], "WhitespaceRule", "preserve");
opts = setvaropts(opts, 1, "TrimNonNumeric", true);
opts = setvaropts(opts, 1, "ThousandsSeparator", ",");
opts = setvaropts(opts, [3, 4, 5, 6], "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";

% Import the data
datafile = readtable(filename, opts);

end
% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------

