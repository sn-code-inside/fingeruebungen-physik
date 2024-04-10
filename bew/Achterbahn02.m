% -------------------------------------------------------------------------
% Achterbahn02.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet die zeitliche Entwicklung  eines Massenpunktes auf
% einer vorgegebenen Kurve (1+tanh(x)) unter Einfluss der Gravitation und
% Berücksichtigung Reibung. 
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

% Parameter
g  = 9.81;                               % Schwerebeschleunigung

% AB
% Objekt erreicht sicher Plateau
H  = 5;                                  % Höhe Berg
x0   = -H;
hsx0 = hp(x0,H);
sicher  = 1.025;
dx0  = sicher*sqrt(2*g*H/(1+hsx0^2));     %inkl. Sicherheitspuffer 
z0   = bahn(x0,H);
dz0  = hsx0*dx0;
AB = [x0, dx0, z0, dz0];

NPoints =1000;
tspan  = linspace(0,4*sqrt(2*H/g)/sicher,NPoints);

% ODE 45
opt=odeset('AbsTol',1.e-10,'RelTol',1.e-8);           
[t,Y]=ode45(@(t,Y)DGL(t,Y,H,g),tspan,AB,opt); 
x  = Y(:,1);
vx = Y(:,2);
z  = Y(:,3);
vz = Y(:,4);
T  = 0.5*(vx.^2+vz.^2);
U  = g*z;

% Ableitungen
hs  = hp(x,H);
h2s = h2p(x,H);
fac = 1./(1+hs.^2);

% Beschleunigungen
ax = -(h2s.*vx.^2+g).*hs.*fac;
az = h2s.*vx.^2 + hs.*ax ;
% Normierung für graphische Darstellung
agx = ax*H/10/g;
agz = az*H/10/g;
Nz  = (az+g);
Ngz = Nz*H/10/g;
Ngx = agx;
N   = sqrt(Nz.^2+ax.^2)/g;
v   = sqrt(vz.^2+vx.^2);
theta1 = acosd((ax.*vx+Nz.*vz)./(N.^2+v.^2));
theta2 = atand(hs);


%% Bewegung auf Bahnkurve und Kräfte (Simulation)
figure()
hold on;
h=title('Stufe');
set(h,'FontSize',14,'FontWeight','normal'); 
axis([x0, -x0, -0.2*H, 1.2*H]);
ylabel('Höhe h in m','FontSize',14)
xlabel('x in m','FontSize',14)
set(gca,'FontSize',16);
grid on
plot(x,z,'color', Colors(7,:),'LineWidth',2);
for k = 1:NPoints/25:length(t)
        p1 = plot(x(k),z(k),'s','color', Colors(4,:),'LineWidth',2);
%       Geschwindigkeit
%       p2 = line([x(k),x(k)+vx(k)/5],[z(k),z(k)+vz(k)/5],...
%                     'color', Colors(3,:),'LineWidth',2);
%       Beschleunigung
        p2 = line([x(k),x(k)+agx(k)],[z(k),z(k)+agz(k)],...
                        'color', Colors(4,:),'LineWidth',2);
%       Normalkraft
        p3 = line([x(k),x(k)+agx(k)],[z(k),z(k)+Ngz(k)],...
                        'color', Colors(9,:),'LineWidth',2);
    pause(0.1);
    p1. Visible = 'off';
%     p2. Visible = 'off';
%     p3. Visible = 'off';
    grid on
end
legend([p2,p3],'Beschleunigung', 'Normalkraft',...
               'location','southeast','numcolumns',1);
legend box off
set(gca,'FontSize',16);




%% Beschleunigungen

figure()
h=title('Stufe');
set(h,'FontSize',14,'FontWeight','normal'); 
% yyaxis left
hold on;
% ylabel('Geschwindigkeit in m/s','FontSize',14)
% plot(t,vx,'color', Colors(3,:),'LineWidth',2);
% plot(t,vz,'color', Colors(4,:),'LineWidth',2);
ylabel('Normalkraft in m x g','FontSize',14)
% p(1)=plot(t,ax.^2+Nx/g,'color', Colors(3,:),'LineWidth',2);
% p(2)=plot(t,az/g,'color', Colors(4,:),'LineWidth',2);
p(1)=plot(t,N,'color', Colors(3,:),'LineWidth',2);
p(2)=plot(t,z,'color', Colors(7,:),'LineWidth',2);
axis([0, 2, -2, 10]);
xlabel('t in s','FontSize',14)
grid on
yyaxis right
% hold on
% plot(t,x,'color', Colors(7,:),'LineWidth',2);
axis([0, 2, -20, 100]);
p(3)=plot(t,theta1,'color', Colors(4,:),'LineWidth',2);
p(4)=plot(t,theta2,'color', Colors(6,:),'LineWidth',2, 'LineStyle', Style(1));
grid on
ylabel('Winkel in °','FontSize',14)
set(gca,'FontSize',16);
legend(p,'Normalkraft', 'z(t)', '\alpha(t)', '\theta(t)', ...
               'location','south','numcolumns',2);
legend box off


%% Energie
figure()  %Energie
Eg = (T(1)+U(1))/100;
plot(x,T/Eg, 'color', Colors(2,:),'LineWidth',2);
hold on
plot(x,U/Eg,'color', Colors(3,:),'LineWidth',2);
plot(x,(T+U)/Eg,'color', Colors(4,:),'LineWidth',2);
axis([x0, -x0, -20, 120]);
ylabel('Energie in % ','FontSize',14)
xlabel('x in m','FontSize',14)
grid on
legend('kin. Energie', 'pot. Energie','Gesamtenergie', ...
               'location','east','numcolumns',1);
legend box off
h=title('Stufe');
set(h,'FontSize',14,'FontWeight','normal'); 
set(gca,'FontSize',16);

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


%%
% Funktionen

function hb = bahn(x,H)
    hb  = 0.5*H*(1+tanh(x));
end

function hprime = hp(x,H)
    hprime = 0.5*H./(cosh(x)).^2;
end

function h2prime = h2p(x,H)
    h2prime = -H*tanh(x)./(cosh(x)).^2;
end

%Differentialgleichungssystem
function  dY  =  DGL(t,Y, H, g)
dY     =  zeros(4,1);  %  Es  muss  ein  Spaltenvektor  zurückgegeben  werden 
hs     =  0.5*H./(cosh(Y(1))).^2;
hss    =  -H*tanh(Y(1))./(cosh(Y(1))).^2;
dY(1)  =  Y(2);
dY(2)  = -(hss*Y(2)^2+g)*hs/(1+hs^2); 
dY(3)  =  Y(4);
dY(4)  = (hss*Y(2)^2-g*hs^2)/(1+hs^2); 
end

% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------



