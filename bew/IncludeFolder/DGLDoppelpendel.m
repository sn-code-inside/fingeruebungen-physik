% -------------------------------------------------------------------------
% DGLDoppelPendel.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% gibt die Ableitungen für die DGL des Doppelpendels zurück
% -------------------------------------------------------------------------

function ders = DGLDoppelpendel(t,w,flag,L1,L2,m1,m2,g)
    % w(1):theta1, w(2):theta1_dot, w(3):theta2, w(4):theta2_dot
    tp = w(1)-w(3); cs=sin(tp); cc=cos(tp);
    ta = m2*L2*cc/(m1+m2)/L1;
    tb = (m2*L2*w(4).^2.*cs+(m1+m2)*g*sin(w(1)))/(m1+m2)/L1;
    tc = (L1*tb.*cc+L1*w(2).^2.*cs-g*sin(w(3)))./(L2-ta.*L1.*cc);
    ders = [w(2);-ta.*tc-tb;w(4);tc];
end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------
