% -------------------------------------------------------------------------
% OberthEffektAnalogie.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Astrodynamik" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% 
% Analogie zum Oberth Effekt als Bewegung durch eien Potentialmulde
%
% -------------------------------------------------------------------------

% Initialisierung
clc
clear 
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors=GetColorLines;
Style = ["-", "-.", ":", "--", ":"];


%Parameter
g = 10;         % in m/s


%% Berechnung des minimalen  Antriebsbedarfs

hA = linspace(0,10,6);
Dh = linspace(0,5,101);

for m= 1:length(hA)
    for k=1:length(Dh)
      hB(m,k) = hA(m)+Dh(k);
      DV(m,k) = -sqrt(2*g*(hB(m,k)-hA(m)))+...
                    sqrt(4*g^2*(hB(m,k)-hA(m))^2+2*g*hA(m));
     end
    tempstr = string(strcat('h_A = ',num2str(hA(m),' %4.1f')));
    lgdstr(m,:) = string(strcat(tempstr,' m'));
end


figure()
for m=1:length(hA)
    plot(hB(m,:) + hA(m),DV(m,:));
    hold on;
end
ylabel('\Delta v')
xlabel('h_B')
grid on
legend(lgdstr,'location','northeast');   
legend box off
title('\Deltav (h_A, \Delta h)')
set(gca,'FontSize',14);

%% Berechnung der Endgeschwindigkeit  Antriebsbedarfs

Dh = linspace(0,20,6);
vB = sqrt(2*g*Dh);
Dv = linspace(10,20,1001);
hA0 = 10;
for m= 1:length(Dh)
    for k=1:length(Dv)
      vfinal(m,k) = sqrt(Dv(k)*Dv(k) - 2*g*hA0 +...
                    2*Dv(k)*sqrt(vB(m)));
       if ~isreal(vfinal(m,k)) 
                  vfinal(m,k) = NaN;  
       end
    end
    tempstr = string(strcat('v_B = ',num2str(vB(m),' %4.1f')));
    lgdstr(m,:) = string(strcat(tempstr,' m/s'));
end


figure()
for m=1:length(vB)
    plot(Dv(:), vfinal(m,:));
    hold on;
end
ylabel('v_{final} in m/s')
xlabel('\Deltav in m/s')
grid on
legend(lgdstr,'location','northwest', 'NumColumns',2);   
legend box off
% title('\Deltav (_A, \Delta h)')
set(gca,'FontSize',14);


% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------