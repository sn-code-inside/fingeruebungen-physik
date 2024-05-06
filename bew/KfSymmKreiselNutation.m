% -------------------------------------------------------------------------
% .m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Nutation des kräftefreine symmetrischen Kreisels
% 
% Beispielrechnungen und Simulation zur Nutation des kräftefreien
% symmetrischen Kreisels.
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

%%
% Variablen: Traegheitsmomente
J12 = 1.0;
J3  = 1.5;
% Integrationszeit
% Komponenten der Rotationsachse im koerperfesten System zu t=0
omega1 = 0.2;
omega2 = 0.2;
omega3 = 1.0;


%%
% Berechnen/Setzen der Anfangsbedingungen
phi0       = 0.;
omega_perp = sqrt(omega1^2 + omega2^2);
thetaN     = atan(J12*omega_perp/(J3*omega3));
deltaAlpha = atan2(omega2,omega1);
psi0       = -deltaAlpha+pi/2.;
Omega      = (J3-J12)/J12*omega3;
tmax = 2*pi/Omega;
tv=linspace(0.,tmax,50);
L  = sqrt((J12*omega1)^2+(J12*omega2)^2+(J3*omega3)^2);
OmegaS     = L/J12;


%% Numerische Berechnung
%
Y0      =[omega1;omega2;omega3;phi0;thetaN;psi0];
options = odeset('AbsTol',1.e-8,'RelTol',1.e-6); 
[t1,Y1] = ode45(@DGL,tv,Y0,options,J12,J12,J3);

%% Analytische Berechnung
%
theta = ones(length(tv),1);
theta = thetaN*theta;
psi   = -Omega*tv+psi0;
phi   = omega_perp*tv/sin(thetaN) + phi0;

%% Graphische Ausgabe

figure()

subplot(1,2,1)
plot(t1,Y1(:,1),'Linewidth',1,'Color',Colors(1,:),'LineStyle',Style{1});
axis([0 tmax -0.5 1.1]);
ylabel('\omega_{\it i} \rm in 1/s ','FontSize',14);
xlabel('{\it t} \rm in s ','FontSize',14);
hold on
plot(t1,Y1(:,2),'Linewidth',1,'Color',Colors(2,:),'LineStyle',Style{1});
plot(t1,Y1(:,3),'Linewidth',1,'Color',Colors(5,:),'LineStyle',Style{1});
legend('\omega_{1}','\omega_{2}','\omega_{3}','FontSize',14,...
       'Location','southwest');
legend box off;   
grid on;
set(gca,'FontSize',16)   
hold off;


% Euler-Winkel
subplot(1,2,2)
plot(t1,rad2deg(Y1(:,4)),'Linewidth',1,'Color',Colors(4,:),'LineStyle',Style{1});
axis([0 tmax -360 360]);
ylabel('\phi,\vartheta,\psi','FontSize',14);
xlabel('{\it t} \rm in s ','FontSize',14);
hold on
plot(t1,rad2deg(Y1(:,5)),'Linewidth',1,'Color','#A55621','LineStyle',Style{1});
plot(t1,rad2deg(Y1(:,6)),'Linewidth',1,'Color',Colors(3,:),'LineStyle',Style{1});
legend('\phi','\vartheta','\psi','FontSize',14,...
       'Location','southwest');
legend box off
grid on
set(gca,'FontSize',16)   

text(tmax*0.5,360-40,string(strcat(' \theta_N = ', ...
                        num2str(rad2deg(thetaN),'%4.1f°'))));
text(tmax*0.5,360-80,string(strcat(' \Omega = ', ...
                        num2str(Omega,'%4.1f 1/s'))));
text(tmax*0.5,360-120,string(strcat(' \Omega_S = ', ...
                        num2str(OmegaS,'%4.2f 1/s'))));
text(tmax*0.5,360-160,string(strcat(' J_3/J_{12} = ', ...
                        num2str(J3/J12 ,'%4.2f '))));

% % Zum Vergleich Analytische Lösung 
% % Winkel theta der Deutlichkeit wegen Faktor 10
% subplot(1,2,2)
% p1=plot(tv,rad2deg(phi),'Linewidth',1,'Color',Colors(4,:),'LineStyle',Style{1});
% hold on
% p2=plot(tv,rad2deg(theta)*10,'Linewidth',1,'Color',[],'LineStyle',Style{1});
% p3=plot(tv,rad2deg(psi),'Linewidth',1,  'Color',Colors(3,:),'LineStyle',Style{1});
% legend('\phi','\vartheta','\psi','FontSize',14,...
%        'Location','northeast');
% legend box off
% ylabel('\phi,\vartheta,\psi','FontSize',14);
% xlabel('{\it t} \rm in s ','FontSize',14);
% axis([0 tmax -360 360]);
% grid on
% set(gca,'FontSize',16)   

% print ("SymmKreiselNutation.pdf","-S652,425","-F:14");
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


%%
% Differentialgleichungssystem 
function dYdt = DGL(t,Y,J1,J2,J3)
    
    dYdt = [(J2-J3)/J1*Y(2)*Y(3)
		    (J3-J1)/J2*Y(3)*Y(1)
		    (J1-J2)/J3*Y(1)*Y(2)
            (Y(1)*sin(Y(6))+Y(2)*cos(Y(6)))/sin(Y(5))
            Y(1)*cos(Y(6))-Y(2)*sin(Y(6))
            Y(3)-cot(Y(5))*(Y(1)*sin(Y(6))+Y(2)*cos(Y(6)))];
end



% % Differentialgleichungssystem 
% function dYdt = DGL(t,Y,J1,J2,J3)
%     
%     dYdt = [(J2-J3)/J1*Y(2)*Y(3)
% 		 (J3-J1)/J2*Y(3)*Y(1)
% 		 (J1-J2)/J3*Y(1)*Y(2)
%      (Y(1)*sin(Y(6))+Y(2)*cos(Y(6)))/sin(Y(5))
%      Y(1)*cos(Y(6))-Y(2)*sin(Y(6))
%      Y(3)-cot(Y(5))*(Y(1)*sin(Y(6))+Y(2)*cos(Y(6)))];
% end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------

