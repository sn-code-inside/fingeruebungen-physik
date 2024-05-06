% -------------------------------------------------------------------------
% Chladni03.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik des Kontinuums" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet die eine Lösung für eine schwingende, rechteckige
% Scheibe durch die Lösung der zugehörigen stationären partiellen
% Diefferentialgleichung.
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
cmap = GetColorMap;
Style = ["-", "-.", ":", "--", ":"];

%% Parameter
a       = 10;                         % erste Kantenlänge in cm
b       = 20;                         % zweite Kantenlänge in cm
rho     = 9.5e-4;                     % Dichte in kg/cm^3
mu      = 10;                         % Torsionsmodul in N/dm^2

%% Erstellen des Lösungsobjekts
% PDGl-Objekt erstellen und Geometrie festlegen
model = createpde();
R1 = [3,4,-a/2,a/2,a/2,-a/2,a/2,a/2,-a/2,-a/2]';
g1 = decsg(R1);
geometryFromEdges(model,g1);

% Bei Bedarf, Rechteck darstellen
%pdegplot(model,EdgeLabels="on")

% Gitter zeigen
hmax = 5e-2;
generateMesh(model,"Hmax",hmax);

% Bei Bedarf zeigen des Gitters
%figure
%pdemesh(model); 
%axis equal

%% Lösen der partiellen Differentialgleichung

% Randedingung definieren: Dirichlet-Randbedungung u = 0
applyBoundaryCondition(model,"dirichlet","Edge", ...
                       1:model.Geometry.NumEdges,"u",0);

% Aufstellen der Differentialgleichung
specifyCoefficients(model,"m",rho,"d",0,"c",mu,"a",0,"f",0);

% Bereich, in dem Eigenwerte gesucht werden
evr = [0,5e4];

% Lösung berechnen
results = solvepdeeig(model,evr);
resultvector = results.Eigenvectors;

%% Darstellen der Lösung
figure()
for i=1:6
    subplot(3,2,i)
    u = resultvector(:,i,1);
    pdeplot(model,"XYData",u,"Mesh","off")
    axis([-0.505*a 0.505*a -0.505*b 0.505*b],'square');
    colormap(cmap)
    xlabel("{\itx} in cm");
    ylabel("{\ity} in cm");
    c = colorbar;
    c.Label.String = '{\itz} in cm';
    set(gca,'FontSize',16,'FontName','Times');
end

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------

