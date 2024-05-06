% -------------------------------------------------------------------------
% Ausgabewerte.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Hilfsfunktion 
% -------------------------------------------------------------------------
function Ausgabewerte(t,a,v,W,eps1,eps2,P)
    fprintf("------------------------------------\nEffizienzparameter ");
    fprintf("\n------------------------------------\n");
    fprintf(1, 't_abw \t %6.2f\t  s \n', t);
    fprintf(1, 'alpha \t  %5.1f\t  ° \n', a);
    fprintf(1, 'v_abw \t  %5.1f\t  m/s \n', v);
    fprintf(1, 'W     \t  %5.1f\t  m \n', W);
    fprintf(1, 'eps_W \t  %5.1f\t  %%\n', eps1*100);
    fprintf(1, 'eps_T \t  %5.1f\t  %%\n', eps2*100);
    fprintf("\n");
    fprintf(1, 'LG    \t  %5.1f\t  m  \n', P.LG);
    fprintf(1, 'lW    \t  %5.1f\t  m  \n', P.lW);
    fprintf("\n");
    fprintf(1, 'h     \t  %5.1f\t  m  \n', P.h);
    fprintf(1, 'B     \t  %5.1f\t  m  \n', P.B);
    fprintf(1, 'b     \t  %5.1f\t  m  \n', P.b);
    fprintf(1, 'LT    \t  %5.1f\t  m  \n', P.b+P.B);
    fprintf("\n");
    fprintf(1, 'MG \t    %7.1f\t  kg \n', P.MG);
    fprintf(1, 'mW    \t  %5.1f\t  kg \n', P.mW);
    fprintf(1, 'mT    \t %5.1f\t  kg \n', P.mT);
    fprintf("\n");
end

