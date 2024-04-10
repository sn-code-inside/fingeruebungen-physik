%%-----------------------------------------------------------------------%%
%%-----------------------------------------------------------------------%%
% StrDMS.m
% MATLAB Programm zu Fingeruebungen fuer Physik von Michael Kaschke
% Kapitel Himmelsmechanik
% Alle Rechte beim Autor
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% Das vorliegende Programmaterial ist mit keiner Verpflichtung oder 
% Garantie irgendwelcher Art verbunden. Die Benutzung der Software ist nur 
% zu privaten oder nicht-kommerziellen Zwecken, wie Bildung und Lehre 
% gestattet.
%%-----------------------------------------------------------------------%%
%%-----------------------------------------------------------------------%%

%-------------------------------------------------------------------------
% StrDMS
% Umrechnung eines Winkels von dezimaler Darstellung zu Grad, Bogenminute 
% und Bogensekunde
%-------------------------------------------------------------------------
%% 

function xStrDMS = StrDMS(DG)
% Eingabe:
%   DG        Grad in dezimaler Darstellung
% Ausgabe:  xStrDMS[G,M,S]
%   G         Grad
%   M         Minuten
%   S         Sekunden 
  VZ = sign(DG);
  Sec  = 3600*abs(DG);
  xDMS = [floor(Sec/3600),floor(rem(Sec,3600)/60),rem(rem(Sec,3600),60)];
  xStrDMS = string(num2str(VZ*xDMS(1),'%+04d°'))+ ...
            string(num2str(abs(xDMS(2)),'%02d´'))+ ...
            string(num2str(abs(xDMS(3)),'%05.2f´´')); 
end
%%-----------------------------------------------------------------------%%
% Ende Funktion
%-------------------------------------------------------------------------
