% -------------------------------------------------------------------------
% PlotCircle.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Hilfsfunktion zum Zeichnen eines Kreises
% -------------------------------------------------------------------------

function PlotCircle (xM,yM,rC,Col,LW)
    u=linspace(0,360,360);
    nx=zeros(360);
    ny=zeros(360);
    Nmax = round(30/rC);
    if Nmax > 100 Nmax = 100; end
    for k=1:Nmax
        nx= rC*k*cosd(u)/Nmax+xM;
        ny= rC*k*sind(u)/Nmax+yM;
        plot(nx,ny,'LineWidth',LW,'Color',Col);
    end
end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------
