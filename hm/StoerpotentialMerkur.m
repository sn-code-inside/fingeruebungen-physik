% -------------------------------------------------------------------------
% StoerpotentialMerkur.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet Stoerpotentiale auf Merkur 
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

% Import Data
MassenPlaneten = importfile('MassenPlaneten.csv')

radius = 149.6e09*MassenPlaneten.aP;
mass   = 1e24*MassenPlaneten.Masse;
G      = 6.672e-11;
a =1;
MS     = 1.989e30;

% Begin Rechnungen

lambda=mass/2/pi./radius;
for i=2:8
    FP(i)=G*pi*mass(1)*lambda(i)*a/(radius(i)*radius(i)-a*a);
end

F0=-G*MS*mass(1)/a^2;
DeltaF = sum(FP);
a=linspace(0.01*radius(1),1.5*radius(1),10000);
lambda=mass/2/pi./radius;
for i=2:8
    FPr(i,:) = G*pi*mass(1)*lambda(i).*a./(radius(i)*radius(i)-a.*a);
end

FPr_0=G*MS*mass(1)./a.^2;
FPr_all=sum(FPr);


%--------------------------------------------------------------------------
% Graphische Ausgabe
figure()
LgdStr = MassenPlaneten.Name;
LgdStr(1,:) = 'Sonne';
LgdStr(7,:) = 'Summe';

h(1) = semilogy(a/radius(3),FPr_0,'Color', Colors(10,:),'LineWidth',3);
hold on;
for iPlot = 2:6
  h(iPlot) = semilogy(a/radius(3),FPr(iPlot,:),'Color', Colors(iPlot,:), ...
      'LineStyle','-.','LineWidth',2);
  hold on;
end
h(7)=semilogy(a/radius(3),FPr_all,'Color', Colors(9,:),'LineWidth',3);

%   loglog(a/radius(1),FPr_0,'Color', Colors(10,:));
%   loglog(a/radius(1),FPr_all,'Color', Colors(10,:));
grid on;
% grid minor;
header2 = sprintf('Gravitationskraefte auf Merkur');
lgd = title(header2);
lgd.FontSize = 16;
lgd.FontWeight = 'normal';
xlabel('r in AE')
ylabel('Anziehungskraefte in N');
ylim([1e12 1e24]);

lgd=legend(h(1:7),LgdStr(1:7),'Location','best','NumColumns',1);
lgd.FontSize=16;
legend boxoff;
ax=gca;
set(gca,'XGrid','on', 'YGrid', 'on', 'Fontsize', 14, 'linewidth', 1);

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Functionen
% -------------------------------------------------------------------------

function MassenPlaneten = importfile(filename, dataLines)
%IMPORTFILE Import data from a text file
%  MASSENPLANETEN1 = IMPORTFILE(FILENAME) reads data from text file
%  FILENAME for the default selection.  Returns the data as a table.
%
%  MASSENPLANETEN1 = IMPORTFILE(FILE, DATALINES) reads data for the
%  specified row interval(s) of text file FILENAME. Specify DATALINES as
%  a positive scalar integer or a N-by-2 array of positive scalar
%  integers for dis-contiguous row intervals.
%
%  Example:
%  MassenPlaneten1 = importfile("E:\GITHubs\PhysicsUncut\MatlabDateien\hm\Data\MassenPlaneten.csv", [2, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 29-Jan-2023 10:26:10

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [2, Inf];
end

%% Setup the Import Options
opts = delimitedTextImportOptions("NumVariables", 3);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ";";

% Specify column names and types
opts.VariableNames = ["Name", "Masse", "aP"];
opts.VariableTypes = ["string", "double", "double"];
opts = setvaropts(opts, 1, "WhitespaceRule", "preserve");
opts = setvaropts(opts, 1, "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
MassenPlaneten = readtable(filename, opts);

end
