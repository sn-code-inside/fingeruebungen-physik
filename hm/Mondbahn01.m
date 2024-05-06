% -------------------------------------------------------------------------
% Mondbahn01.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet für ein einfaches koplanares ZK-Planet-Mond-System
% die Mondbahn bezogen auf den Zentralkörper für verschiedene Verhältnisse
% der Umlaufzeiten TM (Mond um Planet) und TP (Planet um ZK).
% Das Verhältnis der Radien der Umlaufbahnen ist aM/aP = 0.1.
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

% syms aP aM Omega omega t x y rBahn rho 'real'
% 
% 
% %% Symbolische Lösung Krümmungsradius
% 
% % Koordinaten und Ableitungen
% 
% xM = aP*cos(Omega*t) + aM*cos(omega*t);
% yM = aP*sin(Omega*t) + aM*sin(omega*t);
% 
% dx = diff(xMoon,t);
% dy = diff(yMoon,t);
% 
% ddx = diff(dx,t);
% ddy = diff(dy,t);
% 
% rBahn = sqrt(xM.^2+yM.^2);
% 
% rho = sqrt((dx.^2+dy.^2).^3)/(dx.*ddy-ddx.*dy);

%% Einfache koplanare Modellrechnung
%
AE     = 149.7;        % AE in Mio km

for k=1:3
    switch k   
        case 1 % Fall 1 aP*OmegP^2  > aM*omegM^2 (vglbar Erde-Mond)
            TP     = 365.25;                % in Tagen
            TM     = 29.53;
            aP     = 1;                     % Halbachse Planet in AE
            aM     = 0.384/AE;              % Halbachse Mond in AE
            aP     = 1;                     % Halbachse Planet in AE
            aM     = 0.1;                   % Halbachse Mond in AE
            TP     = 20;
            TM     = 10;
            str1 = num2str(TM/TP,'%5.3f'); 
            lgdstr(1,:) = sprintf(strcat('Fall 1 TM/TP = ',str1));
        case 2 % Fall 2 aP*OmegP  < aM*omegM (vglbar Jupiter-Io)
            TP     = 20;
            TM     = 1;
            aP     = 1;                     % Halbachse Planet in AE
            aM     = 0.1;                   % Halbachse Mond in AE
            ratio  = (aP/TP)/(aM/TM);       % muss < 1 sein
            str1 = num2str(TM/TP,'%5.3f'); 
            lgdstr(2,:) = sprintf(strcat('Fall 2 TM/TP = ',str1));
         case 3 % Fall 3 aP*OmegP^2  < aM*omegM^2  (vglbar Saturn-Titan)
            TP     = 20;
            TM     = 2;
            aP     = 1;                     % Halbachse Planet in AE
            aM     = 0.1;                   % Halbachse Mond in AE
            str1 = num2str(TM/TP,'%5.3f'); 
            lgdstr(3,:) = sprintf(strcat('Fall 3 TM/TP = ',str1));
    end 
%     ratio1  = TM^2/TP^2   
%     ratio2  = aM/aP   
%     ratio3  = TM/TP   
    % Berechnung
    Omega  = 2*pi/TP;    % Omega Planet
    omega  = 2*pi/TM;    % omega Mond
    t      = linspace(0,2*pi/min(Omega,omega),1000);
    xP(k,:)  = aP*cos(Omega*t(:));
    yP(k,:)  = aP*sin(Omega*t(:));
   
    xM(k,:) = aP*cos(Omega*t(:)) + aM*cos(omega*t(:));
    yM(k,:) = aP*sin(Omega*t(:)) + aM*sin(omega*t(:));
    dx  = - Omega*aP*sin(Omega*t) - omega*aM*sin(omega*t);
    dy  = + Omega*aP*cos(Omega*t) + omega*aM*cos(omega*t);
    ddx = - omega^2*aM*cos(omega*t) - Omega^2*aP*cos(Omega*t);
    ddy = - omega^2*aM*sin(omega*t) - Omega^2*aP*sin(Omega*t);
    det = (dx.*ddy-ddx.*dy);
    % Krümmungsradius
    rho(k,:) = abs((dx.*dx + dy.*dy).^(3/2)./det);
    % Planetenbahn, Mondbahn
    r(k,:)   = sqrt((xM(k,:).^2 + yM(k,:).^2));
    rP(k,:)  = sqrt((xP(k,:).^2 + yP(k,:).^2));
    % Krümmungsmittelpunk
    xKM(k,:) = xM(k,:) - dy.*(dx.^2+dy.^2)./det;
    yKM(k,:) = yM(k,:) + dx.*(dx.^2+dy.^2)./det;
    % Konvexität
    y2s(k,:) = det./dx.^3;
end
%% Ausgabe
t =t/max(t);

figure('Name',' Mondbahn zu ZK');
hold on
h(1) = plot(t,rP(1,:),'color',Colors(3,:),'LineWidth',2);
for k=1 :3 
    h(1+k) = plot(t,r(k,:),'color',Colors(9+k,:),'LineWidth',2,...
        'LineStyle',Style(k));
end
legend(h,'Planetenbahn',lgdstr(1,:),lgdstr(2,:),lgdstr(3,:),...
       'location','south')
legend box off
title('Abstand Mond zu ZK (koplanare Kreisbahn)')
grid on;
xlabel('Zeit');
ylim([0.8 1.2]);
ylabel('Abstand Mond-ZK in km');
set(gca,'XGrid','on', 'YGrid', 'on', 'Fontsize', 14, 'linewidth', 1);

figure('Name',' Mondbahn ');
hold on
h(1) = plot(xP(1,:),yP(1,:),'color',Colors(3,:),'LineWidth',2);
for k=1 :3 
    h(1+k) = plot(xM(k,:),yM(k,:),'color',Colors(9+k,:),'LineWidth',2,...
        'LineStyle',Style(k));
end
grid on;
axis equal
axis square
title('Mondbahn')
PlotCircle (0,0,0.05,Colors(10,:),3);  % große Masse
legend(h,'Planetenbahn',lgdstr(1,:),lgdstr(2,:),lgdstr(3,:),...
       'location','eastoutside')
legend box off
xlim([-1.2 1.2]);
ylim([-1.2 1.2]);
xlabel('x in AE');
ylabel('y in AE');
set(gca,'XGrid','on', 'YGrid', 'on', 'Fontsize', 14, 'linewidth', 1);
figure('Name',' Evolute ');
hold on
for k=1:3 
    hk(k) = plot(xKM(k,:),yKM(k,:),'color',Colors(9+k,:),'LineWidth',2,...
        'LineStyle',Style(k));
end
grid on;
axis equal
axis square
title('Krümmungsmittelpunkte (Evolute der Mondbahn)')
legend(hk,lgdstr(1,:),lgdstr(2,:),lgdstr(3,:),...
       'location','eastoutside')
legend box off
xlabel('x in AE');
ylabel('y in AE');
set(gca,'XGrid','on', 'YGrid', 'on', 'Fontsize', 14, 'linewidth', 1);



figure('Name',' Krümmungsradius');
title('lokaler Krümmungsradius')
hold on
for k=1 :3 
    plot(t,rho(k,:),'color',Colors(9+k,:),'LineWidth',2,...
        'LineStyle',Style(k));
end
legend(lgdstr(1:3,:), 'location','eastoutside')
legend box off
ylabel('rho in AE');
set(gca,'XGrid','on', 'YGrid', 'on', 'Fontsize', 14, 'linewidth', 1);


figure('Name',' Konvexität ');
hold on
Width = [2,1,1];
for k=1:3 
    hk(k)= plot(t,y2s(k,:),'color',Colors(9+k,:),'LineWidth',...
        Width(k),'LineStyle',Style(k));
end
grid on;
title('Konvexität')
legend(hk,lgdstr(1,:),lgdstr(2,:),lgdstr(3,:),...
       'location','eastoutside')
legend box off
ylim([-1000 1000]);
xlim([0 0.5]);
xlabel('Zeit');
ylabel('Konvexität');
set(gca,'XGrid','on', 'YGrid', 'on', 'Fontsize', 14, 'linewidth', 1);

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------

