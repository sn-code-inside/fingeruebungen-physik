% -------------------------------------------------------------------------
% KfSymmKreiselSimul.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm zur Simulation 
% des kräftefreien symmetrischen Kreisels im körperfesten System
% 
% Gestrichene Größen sind ohne Strich bezeichnet
% z.B.  omega' = omega
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];
Marker = ['o','d','o','s','+'];

%% Parameter
%
J1=1; J2=J1; J3=1.5; gam=-(J1-J3)/J1;           % Trägheitsmomente
omega3=1; 
omega20 = 0.2;
omega10 = 0.2;
omegaperp = sqrt(omega10^2+omega20^2);
omegaNu=gam*omega3;                             % OmegaNutation = Omega
phi0=0; 
tmax=2*pi/omegaNu; NPts=100; 
tv=linspace(0,tmax,NPts);                       % Zeit
omega1=omegaperp*cos(omegaNu*tv+phi0);
omega2=omegaperp*sin(omegaNu*tv+phi0);          % omega1, omega2
% Drehimpulse
L1=J1*omega1; 
L2=J2*omega2; 
L3=J3*omega3;                
v1=max(L1); v2=max(L2); v3=max(L3*(1+0.1));     % Darstellungsbereich
% Winkel zw Omega und Omega3 (Polkegel)
alpha_PK=atan(omegaperp/omega3);                         
% Winkel zw Omega und L (Rastpolkegel/Spurkegel)
alpha_RPK =acos((J1*omegaperp^2+J3*omega3^2)/...               
         (sqrt((omegaperp^2+omega3^2)*(J1^2*omegaperp^2+J3^2*omega3^2))));
% Winkel zw Omega3 und L
theta_N=atan(J1*tan(alpha_PK)/J3);                  
% Präzession  von Omega um L (Larmorfrequenz)
OmegaS = J3*omega3/J1/cos(theta_N);             % Omega'
omega=sqrt(omegaperp^2+omega3^2);               % Betrag von omega
L=sqrt(J1^2*omegaperp^2+J3^2*omega3^2);         % Betrag von L

%% Graphische Ausgabbe
figure()
for i=1:NPts
    clf
%   axis ([-v1,v1,-v2,v2,0,v3])
    axis ([-1,1,-1,1,0,1.5])
    h(1)=line([0,omega1(i)],[0,omega2(i)],[0,omega3],'color',Colors(4,:),...
             'LineStyle','-.','linewidth', 1.5);                % omega 
    h(2)=line([0,L1(i)],[0,L2(i)],[0,L3],'color', Colors(2,:),...
             'LineStyle','-','linewidth', 1.5);                 % L 
    h(3)=line([0,0],[0,0],[0,omega3],'color', Colors(11,:),...
             'linewidth',1.5); % Figurenachse
    grid on   
    pause(0.05)
end
hold on
h(4)=plot3 (omega1,omega2,omega3*(tv+0.01)./(tv+0.01),...
            'LineStyle','--','color', Colors(4,:)); 
% omega  Rastpolkegel 
h(5)=plot3 (L1,L2,L3*(tv+0.01)./(tv+0.01),...
            'LineStyle','-','color', Colors(2,:)); 
% L Nutationskegel
hh=legend(h,'\omega','L','\omega_3','Rastpolkegel','Nutationskegel',...
         'location','best');
legend box off
set(hh,'FontSize',14);
str1=cat(2,'Kräftefreier Kreisel, J_{12}=',num2str(J1,3),...
           'kgm^2, J_3=',num2str(J3,3),'kgm^2',', \omega=',num2str(omega,3)...
           ,'rad/s, L=',num2str(L,3),'kgm^2/s');
str1=cat(2,'Kräftefreier Kreisel, J_{12}=',num2str(J1,3),...
           'kgm^2, J_3=',num2str(J3,3),'kgm^2',' \omega_3=',num2str(omega3,3)...
           ,'rad/s'); 
title(str1,'FontSize',14, 'FontWeight','normal');
xlabel('\it x','FontSize',14');
ylabel('\it y','FontSize',14), 
zlabel('\it z','FontSize',14)
str2=cat(2,'\alpha_{RPK}=',...
            num2str(rad2deg(alpha_PK),3),'^o, \alpha_{PK}= ',...
            num2str(rad2deg(alpha_RPK),3),'^o, \theta_N=',...
            num2str(rad2deg(theta_N),3),'^o');
str3=cat(2,'\Omega=',num2str(omegaNu,3),'rad/s, \Omega´=',...
            num2str(OmegaS,3),'rad/s');
text (-0.8,1,0.4,str2,'FontSize',14)
text (-0.8,1,0.1,str3,'FontSize',14)
set(gca, 'FontSize',14);
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------

