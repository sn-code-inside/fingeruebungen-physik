% -------------------------------------------------------------------------
% Achterbahn03.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet die zeitliche Entwicklung eines Massenpunktes auf
% einer parametrisierten Kurve (Trisectrix) unter Einfluss der Gravitation.
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

% Anfangsbedingungen, Parameter
b  = 10;                             % Parameter für Trisectrix
H  = 3*b;                            % Ungefähre Höhe Trisectrix
g  = 9.81;                           % Schwerebeschleunigung


NPoints = 1000;
u     = linspace(-2*pi,2*pi,NPoints);
% Grenzfall Looping mit befestigten Wagen
tspan = linspace(0,6*pi,NPoints);   
% Grenzfall Looping mit freiem, nach unten hängendem Wagen
tspan = linspace(0,3*pi,NPoints);

u0   = -2*pi;

% Grenzfall Looping mit befestigtem Wagen 
du0  = 0.804725*pi;  

% Grenzfall Looping mit freiem, nach unten hängendem Wagen
du0  = 0.8625*pi;      
x0   = trisectx(b,u0);
dx0  = du0*xp(b,u0);
z0   = trisectz(b,u0);
dz0  = du0*zp(b,u0);

AB = [u0, du0];
opt=odeset('AbsTol',1.e-9,'RelTol',1.e-7);           
[t,Y]=ode45(@(t,Y)DGL(t,Y,b,g),tspan,AB,opt); 
u  = Y(:,1);
du = Y(:,2);
x  = trisectx(b,u);
z  = trisectz(b,u);
vx = xp(b,u).*du;
vz = zp(b,u).*du;
T  = 0.5*(vx.^2+vz.^2);
U  = g*z;

% Ableitungen
xs = xp(b,u);
x2s= x2p(b,u);
zs = zp(b,u);
z2s= z2p(b,u);

ddu= -(du.^2.*(xs.*x2s + zs.*z2s)-g*zs)./(xs.^2+zs.^2);

% Beschleunigungen
ax = x2s.*du.^2 + xs.*ddu;
az = z2s.*du.^2 + zs.*ddu;

% Normierung für graphische Darstellung
agx = ax*H/10/g;
agz = az*H/10/g;
Nz  = (az+g);
Ngz = Nz*H/10/g;
Ngx = agx;
N   = sqrt(Nz.^2+ax.^2)/g;
v   = sqrt(vz.^2+vx.^2);
theta1 = acosd((ax.*vx+Nz.*vz)./(N.^2+v.^2));
theta2 = atand(zs./xs);

%% Bewegung auf Bahnkurve und Kräfte (Simulation)
figure()
hold on;
h=title('Trisectrix');
set(h,'FontSize',14,'FontWeight','normal'); 
axis([x0, -x0, z0, 1.2*max(z)]);
ylabel('Höhe h in m','FontSize',14)
xlabel('x in m','FontSize',14)
set(gca,'FontSize',16);
grid on
plot(x,z,'color', Colors(7,:),'LineWidth',2);
for k = 1:NPoints/25:length(t)
        p1 = plot(x(k),z(k),'s','color', Colors(4,:),'LineWidth',2);
% %       Geschwindigkeit
%         p2 = line([x(k),x(k)+vx(k)/5],[z(k),z(k)+vz(k)/5],...
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
               'location','northeast','numcolumns',1);
legend box off
set(gca,'FontSize',16);

%% Energie
figure()  % Energie
if sign(U(1))< 0
    fac = -1;
else
    fac = 1;
end
Eg = (T(1) + U(1) + fac*U(1))/100;   
plot(x,T/Eg, 'color', Colors(2,:),'LineWidth',2);
hold on
plot(x,(U+fac*U(1))/Eg,'color', Colors(3,:),'LineWidth',2);
plot(x,(T+U+fac*U(1))/Eg,'color', Colors(4,:),'LineWidth',2);
axis([x0, -x0, -20, 120]);
ylabel('Energie in % ','FontSize',14)
xlabel('x in m','FontSize',14)
grid on
legend('kin. Energie', 'pot. Energie','Gesamtenergie', ...
               'location','east','numcolumns',1);
legend box off
h=title('Trisectrix');
set(h,'FontSize',14,'FontWeight','normal'); 
set(gca,'FontSize',16);


%% Beschleunigungen

figure()
h=title('Trisectrix');
set(h,'FontSize',14,'FontWeight','normal'); 
% yyaxis left
hold on;
% ylabel('Geschwindigkeit in m/s','FontSize',14)
% plot(t,vx,'color', Colors(3,:),'LineWidth',2);
% plot(t,vz,'color', Colors(4,:),'LineWidth',2);
ylabel('Normalkraft in m x g','FontSize',14)
p(1)=plot(x,N,'color', Colors(3,:),'LineWidth',2);
line ([-x0,x0], [1,1]);
line ([-x0,x0], [2,2]);
axis([x0, -x0, -2, 8]);
xlabel('x in m','FontSize',14)
grid on

yyaxis right
% hold on
% plot(t,x,'color', Colors(7,:),'LineWidth',2);
axis([x0, -x0, -10, 40]);
p(2)=plot(x,sqrt(vz.^2+vx.^2),'color', Colors(7,:),'LineWidth',2);
p(3)=plot(x,z,'color', Colors(7,:),'LineWidth',1);
line ([-x0,x0], [0,0]);
grid on
ylabel('Höhe in m','FontSize',14)
set(gca,'FontSize',16);
legend(p,'Normalkraft', 'v(x)', 'z(x)', ...
         'location','south','numcolumns',2);
legend box off

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


%%
% Funktionen

function xw = trisectx(b,u)
    xw  = b*u.*(u.^2-3)./(u.^2+1);
end
function zw = trisectz(b,u)
    zw  = - b.*(u.^2-3)./(u.^2+1);
end

function xprime = xp(b,u)
    xprime = b*(u.^4+6*u.^2-3)./(u.^2+1).^2;
end

function zprime = zp(b,u)
    zprime = 8*b*u./(u.^2+1).^2;
end

function x2prime = x2p(b,u)
    NN = (u.^2+1);
    x2prime = 6*b*u.*(1-(7/3*u.^2-3)./NN + (4/3*u.^2.*(u.^2-3))./NN.^2)./NN;
end

function z2prime = z2p(b,u)
    NN = (u.^2+1);
    z2prime = 2*b*(1-(5*u.^2-3)./NN + (4*u.^2.*(u.^2-3))./NN.^2)./NN;
end


% Differentialgleichungssystem
function  dY  =  DGL(t,Y,b,g)
dY  =  zeros(2,1);  %  Es  muss  ein  Spaltenvektor  zurückgegeben  werden 
xs     = xp(b,Y(1));
zs     = zp(b,Y(1));
x2s    = x2p(b,Y(1));
z2s    = z2p(b,Y(1));
dY(1)  =  Y(2);
dY(2)  = -(Y(2).^2.*(xs.*x2s + zs.*z2s)-g*zs)./(xs.^2+zs.^2); 
end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------


