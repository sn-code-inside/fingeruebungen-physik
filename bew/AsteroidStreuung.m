% -------------------------------------------------------------------------
% AsteroidStreuung.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Asteroiden-Streuung 
% 
% Programm berechnet Trajektorien für Asteroiden in der Erdbahn
% Analogie zur Rutherford-Streuung ohne Rückstoß
% Am 15. Februar 2013 passierte der ca. 45 m große Asteroid (367943) 
% Duende in einer Entfernung von knapp 28.000 km die Erde, also noch 
% unterhalb der Umlaufbahn der geostationären Satelliten.
% Streuteilchen: Asteroid
% Streuzentrum:  Erde 
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];
Marker = ['o','d','o','s','+'];

% Variablen, Konstanten 
G   = 6.67259e-11;              % Gravitationskonstante in m^3/kg/s^2
mE  = 5.9720e24 ;               % Erdmasse in kg
GE  = mE*G;                     % spez. Gravitationskonstante in m^3/s^2
mM  = 7.3483e22 ;               % Masse Mond in kg
GM  = mM*G;                     % spez. Gravitationskonstante in m^3/s^2
DA  = 25;                       % Durchmesser Asteroid in m
rhoA= 7000;                     % Dichte in kg/m^3
mA  = rhoA*(4/3)*pi*(DA/2)^3;   % Asteroidenmasse in kg
RZ = 0.0; 
vinf = 4000;                    % Geschwindigkeit relativ zur Erde m/s
En  = mA*vinf^2/2;              % Energie in Ws 
b   = 25e6;                     % Stoßparameter in m
x10 = 10*b;                     % Initiale x Position --> unendlich
y10 = b;                        % Initiale y Position = b
% initiales dxdt basierend auf Energy En Richtung Ziel, umgekehrtes Vorzeichen zu x10
d1x0=-sign(x10)*sqrt(2/mA*(En+G/x10)); 
d1y0=0.0;                       % Null y-Geschwindigkeit 
AB=[x10,d1x0,y10,d1y0];         % AB
tmax = 50*b/vinf;
tv = linspace(0,tmax,1000);
RME     = 3e8;                  % Abstand Erde Mond in m
phi0    = deg2rad(310);         % Phasenwinkel Mond Erde
omegaM  = sqrt((GE+GM)/RME^3);  % Umlauffrequenz Mond

[tout,Yout]   = asteroid_path(b, GE,  0*GM, RME, phi0,  AB, tv);
[tout1,Yout1] = asteroid_path(b, GE,  1*GM, RME, phi0,  AB, tv);
[tout2,Yout2] = asteroid_path(b, GE, 10*GM, RME, phi0,  AB, tv);

xM = RME*cos(omegaM*tout+phi0);
yM = RME*sin(omegaM*tout+phi0);
epsilon = b*vinf^2/G/mE;        % Exzentrizität
phi_asy = acosd(-1/epsilon);    % Asymptotenöffnungswinkel

[nr,nc]=size(Yout);             % Zeilen und spalten
L=15;                           % Betrachtungsfenster in units of b

figure(1)
hold on
p(1)=plot(0,0,'o','MarkerFaceColor', Colors(3,:),'MarkerSize',9); % Erde
p(2)=plot(x10/b,y10/b,'o','MarkerFaceColor',Colors(9,:),...
    'MarkerEdgeColor',Colors(9,:),'MarkerSize',6); % Asteroid
p(3)=plot(Yout(:,1)/b,Yout(:,3)/b,'Color', Colors(9,:),...
     'LineWidth',1,'LineStyle',Style(2)); % Asteroid Trajektorie ohne Mond
p(4)=plot(Yout1(:,1)/b,Yout1(:,3)/b,'Color', Colors(8,:),...
     'LineWidth',1,'LineStyle',Style(3)); % Asteroid Trajektorie mit Mond
p(5)=plot(Yout2(:,1)/b,Yout2(:,3)/b,'Color', Colors(8,:),...
     'LineWidth',1,'LineStyle',Style(1)); % Trajektorie 10x Mondmasse 
p(6)=plot(xM/b,yM/b,'Color', Colors(10,:),...
     'LineWidth',2,'LineStyle',Style(3)); % Mondtrajektorie
axis ([-L/3 L -L L/3])

for ki=1:5:round(nr)           % nicht alle Punkte dargestellt
    p1(1)=plot(xM(ki)/b,yM(ki)/b,'o', 'MarkerFaceColor', Colors(10,:),...
                  'MarkerSize',8); % Mond
    p1(2)=plot(Yout(ki,1)/b,Yout(ki,3)/b,'o', 'MarkerFaceColor', Colors(9,:),...
                'MarkerSize',7); % Asteroid
    p1(3)=plot(Yout2(ki,1)/b,Yout2(ki,3)/b,'o','MarkerFaceColor',Colors(8,:),...
        'MarkerEdgeColor',Colors(8,:),'MarkerSize',6); % Asteroid
    pause(0.05)
    p1(1).Visible='off';
    p1(2).Visible='off';
    p1(3).Visible='off';
end
 

xlabel('\it x \rm (Vielfache von \it b)','FontSize',13), ...
       ylabel('\it y \rm (Vielfache von \it b)','FontSize',14)

% Eingangsasymptote
line([-L/3,L],[1,1],'Color',Colors(4,:),'LineStyle',Style(2),...
     'LineWidth',1);
% Ausgangsasymptote
x=Yout(:,1);y=Yout(:,3);               
lx=length(x);                       
sa=(y(lx)-y(lx-1))/(x(lx)-x(lx-1)); % Ausgangsasymptote Anstieg 
xa=-L/3:L;
ya=y(lx)/b+sa*(xa-x(lx)/b);
p(7)=plot(xa,ya,'Color',Colors(4,:),'LineStyle',Style(2),...
     'LineWidth',1);                % Ausgangsasyptote
line([-L/3,L],[0,0],'Color',Colors(15,:),'LineStyle',Style(3),...
     'LineWidth',1);
line([0,0],[-L,L/3],'Color',Colors(15,:),'LineStyle',Style(3),...
     'LineWidth',1);
% Berechne r(t) und rmin
r=sqrt(x.^2+y.^2); rmin=min(r);     
h=legend(p(1:7),' Erde',' Asteroid',' Trajektorie ohne Mond', ' Trajektorie mit Mond',...
' Trajektorie 10x Mondmasse',' Mond ',' Asymptoten','location','bestoutside');
legend box off
set(h,'FontSize',14)

str1(1,:) = string(cat(2,' \it b       \rm  = ',num2str(b/1000,'%4.0f'),' km  '));
str1(2,:) = cat(2,' \it v_{inf} \rm     = ',num2str(vinf/1000,'%3.0f'), ' km/s');
str1(3,:) = string(cat(2,' \it r_{min} \rm    = ',num2str(rmin/1000,'%3.0f'), ' km  '));
str1(4,:) = cat(2,' \phi_{asy}      = ',num2str(phi_asy,'%4.1f'), ' Â°');
str1(5,:) = cat(2,' \epsilon          = ',num2str(epsilon,'%4.3f'), ' ');
for k = 1:5 
   text(L*1.05,-L+(5-k)*2, str1(k,:),'FontSize',12,'FontWeight','normal');
end
title('Numerische Simulation Asteroid','FontWeight','normal','FontSize',12);
grid on;
axis square;
set(gca,'FontSize',14)
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


%% Funktionen
function [tout,Yout] = asteroid_path(b, GE, GM, RME, phi0, AB, tv)
%    opt = odeset('AbsTol',1.e-9,'RelTol',1.e-8); 
%    [tout,Yout] = ode45(@(t,Y)dgl_ruther(t,Y,GE),tv, AB, opt);
    omegaM  = sqrt((GE+GM)/RME^3);  % Umlauffrequenz Mond
    opt = odeset('AbsTol',1.e-9,'RelTol',1.e-7,'Events',@events);
    [tout,Yout,te,ye,ie] = ode23(@(t,Y)dgl_ruther(t, Y, GE, GM, RME, omegaM, phi0),tv,AB,opt);
    function [value,isterminal,direction] = events(t,Y)
            % Locate the time when height passes through zero in a decreasing direction
            % and stop integration.  Here we use a nested function to avoid
            % passing the additional parameter P1 as an input argument.
            r=sqrt(Y(1).^2 + Y(3).^2);
            value = (Y(3)<-5*b) && (r > 10*b)  ;  % detect distance 
            isterminal = 1;   % Stop Integration
            direction  = 1;   % negative Richtung
    end
end

% -------------------------------------------------------------------------
function dY = dgl_ruther(t, Y, GE, GM, RME, omegaM, phi0)
    %  Berechnung im Schwerpunktsystem, Erde als Schwerpunkt angenommen mE>>mM
    %  Es muss ein Spaltenvektor zurückgegeben werden 
    dY     =  zeros(4,1); 
    r  = sqrt(Y(1).^2 + Y(3).^2);
    xM = RME*cos(omegaM*t+phi0);
    yM = RME*sin(omegaM*t+phi0);
    rM = sqrt((Y(1)-xM).^2+(Y(3)-yM).^2); 
    dY(1)  =  Y(2);
    dY(2)  = -GE*Y(1)./r.^3-GM*(Y(1)-xM)./rM.^3;
    dY(3)  =  Y(4);
    dY(4)  = -GE*Y(3)./r.^3-GM*(Y(3)-yM)./rM.^3;
end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------

