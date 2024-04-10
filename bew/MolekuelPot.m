% -------------------------------------------------------------------------
% MolekuelPot.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Molekülpotentiale
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

%%

qb1=3/2;                    % Bindungslänge
qb2=1/2;                    % Bindungslänge
xmin=0.5;            
xmax=2*qb1;            
NPTS=100;                   % Anzahl Punkte
q=[xmin:2/NPTS:2*xmax]; 

% Fall 1 
U1=1./q.^3-1./q.^2;         % Potential 
F1=3./q.^4-2./q.^3;         % Kraft

% Fall 2 
U2=1./q.^12-1./q.^6;        % Potential 
F2=12./q.^13-6./q.^7;       % Kraft

vmin1=min(F1);              % Potentialminimum  
vmin2=min(F2);              % Potentialminimum
vmin=min(vmin1,vmin2);

%%
% Graphik
%
figure();
plot(q,U1,'Color',Colors(3,:),'LineWidth',2,'LineStyle',Style(1));
hold on
plot(q,F1,'Color',Colors(3,:),'LineWidth',2,'LineStyle',Style(5));
plot(q,U2,'Color',Colors(4,:),'LineWidth',2,'LineStyle',Style(1));
plot(q,F2,'Color',Colors(4,:),'LineWidth',2,'LineStyle',Style(5));
line([0,xmax],[0,0],'Color','k','LineStyle',Style(4));
axis([xmin, xmax,1.5*vmin,-1.5*vmin])
ylabel('\it U(q), F(q) \rm','FontSize',14);
xlabel('\it q \rm in \it q\rm_0','FontSize',14);
h=legend('U_1','F_1','U_2','F_2'); set(h,'FontSize',14);
legend box off;
grid on
set(gca,'Fontsize', 16);


%%
%
% Programmteil für die Schwingungslösung des Zweiatom-Potentials 
tmax=2.5;           % Maximalzeit in tau0
qb1=3/2;            % Gleichgewichtsposition
qb2=2^(1/6);        % Gleichgewichtsposition
qi1=0.1;            % Initialwert im Abstand von Gleichgewichtsposition 
qi2=0.1;            % Initialwert im Abstand von Gleichgewichtsposition 
AB1=[qb1+qi1;0.0];  % AB, Anfangsgeschwindigkeit = 0
AB2=[qb2+qi2;0.0];  % AB, Anfangsgeschwindigkeit = 0
% MATLAB's ODE45
opt=odeset('AbsTol',1.e-9,'RelTol',1.e-6);   
[t1,q1]=ode45(@Potmolec1,[0.0,tmax],AB1,opt);% numerical solution
[t2,q2]=ode45(@Potmolec2,[0.0,tmax],AB2,opt);% numerical solution
% HO Näherung
omega1 = sqrt(32/81);
tau01  = 2*pi/omega1;
qH1=qb1+qi1*cos(2*pi*t1);
omega2 = sqrt(72/2^(1/3));
tau02  = 2*pi/omega2;
qH2=qb2+qi2*cos(2*pi*t2);
figure();
axis([0 max(t1) min(qH1) max(qH1)])
subplot(2,1,1)
plot(t1, q1(:,1),'Color',Colors(3,:),'LineWidth',2,'LineStyle',Style(1));
hold on
plot(t1, qH1,'Color',Colors(3,:),'LineWidth',2,'LineStyle',Style(3));
line([0,max(t1)],[qb1,qb1],'color',Colors(3,:),'LineWidth',1,...
     'LineStyle',Style(4));
h=legend('Numer. Lösung 1','Harmon. Osz. 1','Gleichgewichtslage 1');
set(h,'FontSize',14)
xlabel('Zeit (\tau_0)','FontSize',14)
ylabel('Auslenkung \it{q} \rm in (\it{q}\rm_0)','FontSize',14)
legend box off;
grid on
set(gca,'Fontsize', 16);
subplot(2,1,2)
plot(t2, q2(:,1),'Color',Colors(4,:),'LineWidth',2,'LineStyle',Style(1));
hold on
plot(t2, qH2,'Color',Colors(4,:),'LineWidth',2,'LineStyle',Style(3));
line([0,max(t1)],[qb2,qb2],'color',Colors(4,:)','LineWidth',1,...
     'LineStyle',Style(4));
axis([0 max(t1) min(qH2) max(qH2)])
h=legend('Numer. Lösung 2','Harmon. Osz. 2','Gleichgewichtslage 2');
set(h,'FontSize',14)
xlabel('Zeit (\tau_0)','FontSize',14)
ylabel('Auslenkung \it{q} \rm in (\it{q}\rm_0)','FontSize',14)
legend box off;
grid on
set(gca,'Fontsize', 16);


% Berechnungen Helium Molekül
kB = 1.38e-23;
% q0 = 4.0e-09;
q0 = 2.5e-10;
% epstilde = 10*kB;
epstilde = 2.5e-22;
m_eff = 0.5*6.646e-27;
hp = 6.62607015e-34;
hq  = hp/2/pi;
Emin= epstilde;
Tmin = Emin/kB;

omega0 = sqrt(72*epstilde/q0^2/m_eff/2^(1/3));
TP = 2*pi/omega0;



EHO = hq*omega0/2;
THO = EHO/kB;

fprintf('\n Berechnungen zum Helium-Molekül \n');
fprintf('\n eps   in J   : %6.2e', epstilde);
fprintf('\n q0    in m   : %6.2e', q0);
fprintf('\n omega in 1/s : %6.2e', omega0);
fprintf('\n TP    in s   : %6.2e', TP);
fprintf('\n Emin  in J   : %6.2e', Emin);
fprintf('\n EH0   in J   : %6.2e', EHO);
fprintf('\n T in K       : %6.1f', Tmin);
fprintf('\n THO in K     : %6.1f', THO);
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------




%%
%
function dY1dt = Potmolec1( t, q, flag)
    % q(1)-Position, q(2)-Geschwindigkeit
    dY1dt = [ q(2); 81*pi^2*(3./(q(1).^4)-2./(q(1).^3))/8];
end

function dY2dt = Potmolec2( t, q, flag)
    % q(1)-Position, q(2)-Geschwindigkeit
    dY2dt = [ q(2); 4*2^(1/3)*pi^2*(12./(q(1).^13)-6./(q(1).^7))/18];
end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------
