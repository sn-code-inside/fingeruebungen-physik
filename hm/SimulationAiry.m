% -------------------------------------------------------------------------
% AiryModel.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm simuliert das Airy-Kanal-Modell der dynamischen Theorie der
% Gezeiten
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

% Initialisierung
% alles in Einheiten von km und h
RE      = 6370;
TGez    = 12+25/60;
omegGez = 2*pi/TGez;
hT      = 3.5;
g       = (9.81e-3)*3600^2;
cW      = sqrt(hT*g);
lambda  = cW*TGez;

theta = linspace(0,24,24);  %in h
t     = linspace(0,240,240); %in h/10, aller 6 min
q1    = sin(omegGez*t/10);
R     = 1;
delta = pi;

for m=1:length(theta)
    Delh1(m,:) = R*cos(2*(omegGez*t(:)/10-omegGez*theta(m)));
    Delh2(m,:) = R*cos(2*(omegGez*t(:)/10-omegGez*theta(m))-delta/2);
    moon(m)    = omegGez*theta(m);
end

% Simulation Schwingendes System mit und ohen Phasenverschiebung
figure()
xlim([0 24]);
ylim([-2*R 2*R]);
grid on
hold on
xlabel('a.u.');
ylabel('a.u.');
set(gca,'XGrid','on', 'YGrid', 'on', 'Fontsize', 14, 'linewidth', 1);
plot(theta(1),2*R,'o',...
        'MarkerSize',10 , 'Color',Colors(10,:),'LineWidth',2);
h = plot(t(:)/10, Delh1(1,:),'Color',Colors(3,:),'LineWidth',2);
h1= plot(t(:)/10, Delh2(1,:),'Color',Colors(2,:),...
        'LineWidth',2, 'LineStyle',Style(mod(1,2)+1));
h2 = line([theta(1) theta(1)],[-R 2*R],...
         'Color',Colors(4,:),'LineWidth',1,'LineStyle',Style(mod(1,2)+1)); 
for m=2:24
    pause(0.15)
    set(h,'Visible','off');     
    set(h1,'Visible','off');     
    set(h2,'Visible','off');
    h = plot(t(:)/10, Delh1(m,:),'Color',Colors(3,:),...
        'LineWidth',2, 'LineStyle',Style(mod(m,2)+1));
    h1 = plot(t(:)/10, Delh2(m,:),'Color',Colors(2,:),...
        'LineWidth',2, 'LineStyle',Style(mod(m,2)+1));
    plot(theta(m),2*R,'o',...
        'MarkerSize',10 , 'Color',Colors(10,:),'LineWidth',2);
    h2 = line([theta(m) theta(m)],[-R 2*R],...
         'Color',Colors(4,:),'LineWidth',1,'LineStyle',Style(mod(m,2)+1));   
end

% Simulation Schwingendes System mit und ohen Phasenverschiebung
% Ableitung der Phasenverscheibung aus den Parametern Dämpfung und
% Eigenfrequenz
lambda1 = pi*RE;
lambda2 = pi*RE;
omega1 = 2*sqrt(hT*g)/RE;
omega2 = 2*sqrt(hT*g)/RE;
gamma1 = 1.25*omega1;
gamma2 = 1.25*omega2;
A1 =R*omega1^2/sqrt((omega1^2-omegGez^2)^2+4*gamma1^2*omegGez^2);
A2 =R*omega2^2/sqrt((omega2^2-omegGez^2)^2+4*gamma2^2*omegGez^2);
delta1 = atan2(2*gamma1*omegGez,(omega1^2-omegGez^2));
delta2 = atan2(2*gamma2*omegGez,(omega2^2-omegGez^2));
for m=1:length(theta)
  q1(m,:) = cos(2*omegGez*t(:)/10-delta1)*cos(2*theta(m)*omegGez);
  q2(m,:) = sin(2*omegGez*t(:)/10-delta2)*sin(2*theta(m)*omegGez);
  Delh3(m,:) = q1(m,:) + q2(m,:);
end
figure()
xlim([0 24]);
ylim([-2*R 2*R]);
grid on
hold on
xlabel('a.u.');
ylabel('a.u.');
set(gca,'XGrid','on', 'YGrid', 'on', 'Fontsize', 14, 'linewidth', 1);
plot(theta(1),2*R,'o',...
        'MarkerSize',10 , 'Color',Colors(10,:),'LineWidth',2);
h = plot(t(:)/10, Delh3(1,:),'Color',Colors(3,:),'LineWidth',2);
h2 = line([theta(1) theta(1)],[-R 2*R],...
         'Color',Colors(4,:),'LineWidth',1,'LineStyle',Style(mod(1,2)+1)); 
for m=2:24
    pause(0.15)
    set(h,'Visible','off');     
    set(h2,'Visible','off');
    h = plot(t(:)/10, Delh3(m,:),'Color',Colors(3,:),...
        'LineWidth',2, 'LineStyle',Style(mod(m,2)+1));
    plot(theta(m),2*R,'o',...
        'MarkerSize',10 , 'Color',Colors(10,:),'LineWidth',2);
    h2 = line([theta(m) theta(m)],[-R 2*R],...
         'Color',Colors(4,:),'LineWidth',1,'LineStyle',Style(mod(m,2)+1));   
end

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------