% -------------------------------------------------------------------------
% DGLPendel.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% gibt die Ableitungen für die DGL des Pendels zurück
% -------------------------------------------------------------------------

function derivs = DGLPendel( t, x, flag,w0)
    % AbleitungPendel: returns the derivatives for the pendulum's full solution
    % The function pen2_der describes the equations of motion for a 
    % pendulum. The parameter w0, is part of the input
    % Entries in the vector of dependent  variables are:  
    % x(1)-position, x(2)-angular velocity
    derivs = [ x(2); -w0^2*sin(x(1))];
end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------
