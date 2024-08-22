% -------------------------------------------------------------------------
% Rechteckrohr.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik des Kontinuums" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet die den Fluss einer Flüssigkeit durch ein
% rechteckiges Rohr. Das Programm nutzt die Partial Differential Equation
% Toolbox.
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
cmap = GetColorMap;
Style = ["-", "-.", ":", "--", ":"];

%% Parameter
a       = 0.2;                        % erste Kantenlänge
b       = 0.3;                        % zweite Kantenlänge
eta     = 1.;                         % Viskosität
dpdz    = -20.;                       % Druckgradient

%% Erstellen des Lösungsobjekts

% PDGl-Objekt erstellen und Geometrie festlegen
model = createpde();
R1 = [3,4,-a/2,a/2,a/2,-a/2,b/2,b/2,-b/2,-b/2]';
g1 = decsg(R1);
geometryFromEdges(model,g1);

% Bei Bedarf, Rechteck darstellen
% pdegplot(model,EdgeLabels="on")

% Gitter zeigen
hmax = 5e-3;
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
specifyCoefficients(model,"m",0,"d",0,"c",-eta,"a",0,"f",dpdz);

% Lösung berechnen
results = solvepde(model);

%% Darstellen der Lösung
u = results.NodalSolution;
pdeplot(model,"XYData",u,"ZData",u,"Mesh","off")
axis([-0.15 0.15 -0.15 0.15]);
colormap(cmap)
xlabel("{\itx} in m");
ylabel("{\ity} in m");
zlabel("{\itv_z} in m/s");
c = colorbar;
c.Label.String = '{\itv_z} in m/s';
set(gca,'FontSize',16,'FontName','Times');
view([60 20]);

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------

