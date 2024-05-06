% -------------------------------------------------------------------------
% Schaukel01.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Schaukel als Parametrischer Oszillator
% 
% Lösung der DGL für die Schaukel als Parametrischer Oszillator
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

TA    = 1;              % Schwingungsperiode (1 sec)
NPer  = 45;             % Zahl der Perioden  
tmax  = NPer*TA;        % Maximalzeit  (1 min)

omeg0 = 2*pi/TA;        % Eigenfrequenz
omegA = 2*omeg0;        % Frequenz der Pendellängenänderung
phiA  = pi/4 ;          % Phasenmismatch der Pendellängenänderung
Domeg = 0.25*omeg0;     % Änderung der Pendellänge --> Änderung der Eigenfrequenz
eta   = 1e-2;           % quadratische Daempfung

% Anfansgwerte
x0    = 0.01;           % sehr klein, aber verschieden von Null
x0dot = 0;
AB=[x0;x0dot];       % AB für ode45

% Numerische Lösung
opt=odeset('AbsTol',1.e-7,'RelTol',1.e-5);
% Lineare Lösung, keine Reibung
[t1,x1]=ode45(@dgl_schaukel_L,[0.0,tmax],AB,opt,omeg0,0,omegA,Domeg,0); 
minx1=max(rad2deg(min(x1(:,1))),-180);
maxx1=min(rad2deg(max(x1(:,1))),180);

% Nichtlineare Lösung, keine Reibung
[t2,x2]=ode45(@dgl_schaukel_NL,[0.0,tmax],AB,opt,omeg0,0,omegA,Domeg,0); 
minx2=max(rad2deg(min(x2(:,1))),-180);
maxx2=min(rad2deg(max(x2(:,1))),180);

% Nichtlineare Lösung, Reibung
[t3,x3]=ode45(@dgl_schaukel_NL,[0.0,tmax],AB,opt,omeg0,eta,omegA,Domeg,0); 
minx3=max(rad2deg(min(x3(:,1))),-180);
maxx3=min(rad2deg(max(x3(:,1))),180);

% Nichtlineare Lösung, Reibung, Mismatch (Frequenz und Phase)
[t4,x4]=ode45(@dgl_schaukel_NL,[0.0,tmax],AB,opt,omeg0,eta,1.01*omegA,Domeg,phiA); 
minx4=max(rad2deg(min(x4(:,1))),-180);
maxx4=min(rad2deg(max(x4(:,1))),180);


figure();
subplot(2,2,1)
plot(t1,wrapTo180(rad2deg(x1(:,1))),'Color',Colors(2,:), 'LineWidth',2);
hold on
axis([0 tmax minx1 maxx1])
grid on
ylabel('$\theta$ in Grad','interpreter','latex','FontSize',14)
xlabel('$t$ in s','interpreter','latex','FontSize',14)
h=title('Lineare Lösung,keine Reibung');
set(h,'FontSize',12,'FontWeight','normal'); 
set(gca,'FontSize',16);

subplot(2,2,2)
plot(t2,wrapTo180(rad2deg(x2(:,1))),'Color',Colors(3,:), 'LineWidth',2);
hold on
axis([0 tmax minx1 maxx1])
xlabel('$t$ in s','interpreter','latex','FontSize',14)
grid on
h=title('Nichtlineare Lösung, keine Reibung');
set(h,'FontSize',12,'FontWeight','normal'); 
set(gca,'FontSize',16);

subplot(2,2,3)
plot(t3,wrapTo180(rad2deg(x3(:,1))),'Color',Colors(4,:), 'LineWidth',2);
hold on
axis([0 tmax minx3 maxx3])
xlabel('$t$ in s','interpreter','latex','FontSize',14)
ylabel('$\theta$ in Grad','interpreter','latex','FontSize',14)
grid on
h=title('Nichtlineare Lösung, Reibung');
set(h,'FontSize',12,'FontWeight','normal'); 
set(gca,'FontSize',16);

subplot(2,2,4)
plot(t4,wrapTo180(rad2deg(x4(:,1))),'Color',Colors(7,:), 'LineWidth',2);
hold on
axis([0 tmax minx3 maxx3])
xlabel('$t$ in s','interpreter','latex','FontSize',14)
grid on
h=title('Nichtlineare Lösung, Reibung, Mismatch Frequenz und Phase');
set(h,'FontSize',12,'FontWeight','normal'); 
set(gca,'FontSize',16);
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------



%%
function dY = dgl_schaukel_NL(t, Y, omeg0, eta, omegA, Domeg, phiA)
    % Y(1)-Winkel
    % Y(2)-Winkelgeschwindigkeit
    dY = [Y(2); -(omeg0^2-Domeg^2*cos(omegA*t+phiA))*sin(Y(1))-eta*sign(Y(2))*Y(2)^2];
end

function dY = dgl_schaukel_L(t, Y, omeg0, eta, omegA, Domeg, phiA)
    % Y(1)-Winkel
    % Y(2)-Winkelgeschwindigkeit
    dY = [Y(2); -(omeg0^2-Domeg^2*cos(omegA*t+phiA))*Y(1)-eta*sign(Y(2))*Y(2)^2];
end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------
