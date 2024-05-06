% -------------------------------------------------------------------------
% EclipseFlag.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
%-------------------------------------------------------------------------
% Berechnet die Zeiten des Neu/Vollmondes in einem bestimmten Zeitraum und
% bestimmt, ob zu der Neumond/Vollmondzeit eine SoFi/MoFi stattfinden kann.
% Für die Koordinaten Sonne und Mond werden genauere Reihen 
% verwendet.
%
% Das vorliegende Programmaterial basiert in Teilen 
% auf C++ Programmstrukturen und der Beschreibung/Erläuterung in  
%
% "Astronomie mit dem Personal Computer"
%
% von Oliver Montenbruck und Thomas Pfleger (Springer 1999). 
% Genehmigung des Verlages und der Autoren liegt vor.
% -------------------------------------------------------------------------
% EclipseFlag(beta,MPhase)
% Flag zur Berechnung von Finsternissen. Gibt einen String zurueck, der die 
% Moeglichkeit bzw. Art einer Mondfinsternis beschreibt
%
% Eingabe: 
%   beta :       Ekliptikale Breite des Mondes in [rad]
%   MPhase:      Mondphase 1 = Neumond, 2 = Vollmond
% Ausgabe:
%   Flag:        String beschreibt die Finsternis
% -------------------------------------------------------------------------

function Flag = EclipseFlag(beta,MPhase)
    b = abs(beta); 
    switch MPhase
        case 1 % SoFi-Test
           if (b<0.015223) 
               Flag = "Total"; % totale Sonnenfinsternis sicher
           else
             if (b<0.018209) 
               Flag = "Total ?"; % totale/ringf. Sonnenfinsternis moeglich
             else
                 if (b<0.024594)
                    Flag = "Partiell"; %partielle SoFi moeglich
                 else                         
                    Flag = "  "; %keine SoFi moeglich
                 end
             end
           end
        case 2 %MoFiTest
           if (b<0.006351) 
               Flag = "Total"; % totale Mondfinsternis sicher
           else
             if (b<0.009376) 
               Flag = "Total ?"; % totale Mondfinsternis moeglich
             else
                 if (b<0.015533)
                    Flag = "Partiell"; %partielle MoFi sicher
                 else                         
                    Flag = "  "; %keine MoFi moeglich
                 end
             end
           end           
    end
end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------
