% -------------------------------------------------------------------------
% Kreisrohr.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik des Kontinuums" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet die den Fluss einer Flüssigkeit durch ein
% kreisförmiges Rohr und vergleicht die numerische Lösung mit der
% analytischen aus dem Beispiel im Buch. Das Programm nutzt die Partial
% Differential Equation Toolbox.
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
cmap = GetColorMap;
Style = ["-", "-.", ":", "--", ":"];

%% Parameter
eta     = 1.;                         % Viskosität
dpdz    = -20.;                       % Druckgradient

%% Erstellen des Lösungsobjekts
% PDGl-Objekt erstellen und Geometrie festlegen
model = createpde();
geometryFromEdges(model,@circleg);

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
specifyCoefficients(model,"m",0,"d",0,"c",-eta,"a",0,"f",dpdz);

% Lösung berechnen
results = solvepde(model);

%% Berechnen der Differenz zur analytischen Lösung
xy = model.Mesh.Nodes;
analytic = dpdz/(4*eta)*(xy(1,:).^2+xy(2,:).^2-1);

%% Darstellen der Lösung
figure()
subplot(1,2,1)
u = results.NodalSolution;
pdeplot(model,"XYData",u,"ZData",u,"Mesh","on")
axis([-1.2 1.2 -1.2 1.2]);
colormap(cmap)
xlabel("{\itx} in m");
ylabel("{\ity} in m");
zlabel("{\itv_z} in m/s");
c = colorbar;
c.Label.String = '{\itv_z} in m/s';
set(gca,'FontSize',16,'FontName','Times');
view([60 20]);

subplot(1,2,2)
pdeplot(model,"XYData",u(:)-analytic(:),"Mesh","off")
axis([-1.2 1.2 -1.2 1.2 0 6]);
colormap(cmap)
xlabel("{\itx} in m");
ylabel("{\ity} in m");
c = colorbar;
c.Label.String = '{\itv_z}^{numerisch}-{\itv_z}^{analytisch} in {\itm}/{\its})';
set(gca,'FontSize',16,'FontName','Times');

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------

