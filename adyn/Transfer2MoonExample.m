% -------------------------------------------------------------------------
% Transfer2MoonExample.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Berechnet verschiedene Parameter u. Trajektorien für einen Flug zum Mond
% Parameter r0, v0, gamma0 und lambda1 werden vorgegeben.
% 
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
% Style = ["-", "-.", ":", "--", ":"];


%% Initialisierung
% alle Daten in kg, m und s 
G   = 6.67430e-11;                 % G in m^3 /kg /s^2
ME  = 5.977e24;                    % Masse der Erde in kg
MM  = 7.348e22;                    % Masse Mond in kg
RE  = 6.378e6;                     % Erdradius in m
RM  = 1.7374e6;                    % Mondradius in m
rEM = 384400e3;                    % Abstand Erde-Mond in m
omegaM = 2.649e-6;                 % Umlauffrequenz Mond um Erde in rad/s
vM     = omegaM*rEM;               % Bahngeschwindigkeit Mond in m/s
RSOIM  = rEM*(MM/ME)^(2/5);        % SOI-Radius Mond
muE    = G*ME;
muM    = G*MM;


%% Beispiel

% Ausgangsparameter
h0      = 300000;
r0      = h0+RE;
v0      = 10.87*1000;
gamma0  = 0;
lambda1 = deg2rad(30);

% Energie Drehimpuls @ TLI
H0      = v0^2/2-muE/r0;
L0      = v0*r0;

% r, v  @ TLI 
r1      = sqrt(RSOIM^2+rEM^2-2*RSOIM*rEM*cos(lambda1));
v1      = sqrt(2*(H0+muE/r1));

% FPA @ SOI Mond
gamma1  = acos(L0/r1/v1);
gamma1d = rad2deg(gamma1);

% Phasenwinkel Moon @ SOI
phi1    = asin(RSOIM*sin(lambda1)/r1);
phi1d   = rad2deg(phi1);

% TOF
% Ellipsen-Parameter

p0      = L0.^2/muE;
a0      = -muE./H0/2;
exz0    = sqrt(1-p0./a0);
cosups0 = 1;
ups0    = acos(cosups0);
ups0d   = acosd(cosups0);
cosE0   = (exz0+cosups0)./(1+exz0.*cosups0);
E0      = acos(cosE0);
E0d     = acosd(cosE0);
sinE0   = sin(E0);

cosups1 = (p0-r1)./exz0/r1;
ups1    = acos(cosups1);
ups1d   = acosd(cosups1);
cosE1   = (exz0+cosups1)./(1+exz0.*cosups1);
E1      = acos(cosE1);
E1d     = acosd(cosE1);
sinE1   = sin(E1);

fac     = sqrt(a0.^3/muE);
TOF     = fac*(E1- E0 -exz0*(sinE1-sinE0));
TOFh    = TOF/3600;

% Phasenwinkel @ TLI
phi0    = ups1-ups0-phi1-omegaM*TOF;
phi0d   = rad2deg(phi0);

% Übergang Mond KOS
r2      = RSOIM;
v2      = sqrt(v1^2+vM^2-2*v1*vM*cos(gamma1-phi1)); 
eps2    = asin(vM*cos(lambda1)/v2 - v1*cos(lambda1+phi1-gamma1)/v2);
eps2d   = rad2deg(eps2);
H2      = v2^2/2-muM/r2;
% Falls H2 > 0 dann Hyperbelbahn und auch exz2 >1
L2      = v2*r2*sin(eps2);
p2      = L2.^2/muM;
exz2    = sqrt(1+2*H2*L2^2./muM^2);
rPM     = p2/(1+exz2);
vPM     = sqrt(2*(H2+muM/rPM));
hPM     = rPM-RM;

EH2      = acosh((1+r2*(exz2-1)/rPM)/exz2);
a2B      = r2/(exz2*cosh(EH2)-1);
cosups2 = a2B*(exz2-cosh(EH2))/r2;
ups2    = acos(cosups2);
ups2d   = acosd(cosups2);

% Flugzeit um Mond (innerhalb SOI)
TOFM  = 2*sqrt(a2B^3/muM)*(exz2*sinh(EH2)-EH2);
TOFMh = TOFM/3600;


%% Print Ausgabe

fprintf('\n')
fprintf('Vorgaben und Abflugparameter @TLI \n')
fprintf('| h0 (km) | gamma0(°) | v0(km/s) | lambda1(°) |\n');
fprintf('|   %5.1f |  %5.1f    |  %5.1f   | %5.1f      |\n',...
           h0/1000, rad2deg(gamma0), v0/1000, rad2deg(lambda1));
fprintf('\n')
fprintf('Abflugparameter @TLI \n')
fprintf('| a0 (km)  | p0 (km)  |    e0    |  phi0(°)  |\n');
fprintf('| %5.2e | %5.2e |  %6.4f  |  %7.2f  |\n',...
           a0/1000, p0/1000, exz0, phi0d);
fprintf('\n')
fprintf('Parameter auf Transfer zum Mond \n')
fprintf('| v1 (km/s) | v2 (km/s)| ups1(°) | phi1(°) | TOF (h)  |\n');
fprintf('|  %7.3f  | %7.3f  |  %5.1f  |  %5.2f  | %7.2f  |\n',...
           v1/1000, v2/1000, ups0d,phi1d,TOFh);
fprintf('\n')
fprintf('Parameter am Mond \n')
fprintf('| r2 (km)   | eps2(°)  | p2 (km) |    e2   | hPM (km) | TOF-SOI (h) |\n');
fprintf('| %5.2e  | %7.2f  |  %5.1f | %7.4f |  %5.1f  |  %7.2f    |\n ',...
           r2/1000, eps2d, p2/1000,exz2,hPM/1000,TOFMh);
fprintf('\n')



%% Graphik-Ausgabe

ups     = linspace(0,2*pi,3601);
upsM    = linspace(-pi/4,pi/6,3601);
EH      = linspace(-EH2,EH2,1001);
upsH    = 2*atan(sqrt((exz2+1)/(exz2-1))*tan(EH/2));


cups    = cos(ups);
sups    = sin(ups);
cupsM   = cos(upsM);
supsM   = sin(upsM);
domegt  = omegaM*TOF;
upsP    = -(domegt+phi0);
cupsP   = cos(ups-upsP);
supsP   = sin(ups-upsP);


figure('Name','Transferphase')
hold on
hp(1)=plot(RE*cups,RE*sups,'color',Colors(3,:),'linewidth',2);    %Erde
plot(RM*cups+rEM*cos(-domegt),RM*sups+rEM*sin(-domegt), ...
    'color', Colors(5,:),'linewidth',1);                     %Mond @ TLI
plot(RSOIM*cups+rEM*cos(-domegt),RSOIM*sups+rEM*sin(-domegt),...
    'color', Colors(5,:),'linewidth',1);                     %SOI  @ TLI
plot(RM*cups+rEM,RM*sups,'color',Colors(5,:));               %Mond @ arrival
hp(2)=plot(RSOIM*cups+rEM,RSOIM*sups, ...
    'color', Colors(5,:),'linewidth',1);                     %SOI  @ arival
hp(3)=plot(rEM*cos(upsM),rEM*sin(upsM),'color',Colors(5,:)); %Mondbahn

line([0 rEM],[0 0],'linestyle',':');                %Verbindungslinie
line([0 r1*cos(phi1)],[0 r1*sin(phi1)]);            %Phasenwinkel @ SOI
line([0 rEM*cos(domegt)],[0 -rEM*sin(domegt)]);     %Position Mond @ TLI
line([0 5*r0*cos(upsP)],...
     [0 5*r0*sin(upsP)]);                           %Phasenwinkel @ TLI
hp(4)=plot(p0*cups./(1+exz0*cupsP),p0*sups./(1+exz0*cupsP), ...
     'color',Colors(2,:),'linewidth',2);            %Transferbahn
grid on
axis equal
hp2=legend(hp,'Erde','SOI Mond','Mondbahn','Transferbahn',...
    'location','southwest'); 
set(hp2,'FontSize',14,'FontWeight','normal'); 
legend box off;
xlabel('x in m','FontSize',14); 
ylabel('y in m','FontSize',14)
grid on
set(gca,'FontSize',16);



%%

figure('Name','Phase um Mond')
hold on

p2   = L2^2/muM;
exz2 = sqrt(1+2*H2*L2^2./muM^2);
rPM     = p2/(1+exz2);
vPM     = sqrt(2*(H2+muM/rPM));
hPM     = rPM-RM;

hp(1)=plot(RE*cups,RE*sups,'color',Colors(3,:),'linewidth',2);    %Erde
plot(RM*cups+rEM*cos(-domegt),RM*sups+rEM*sin(-domegt), ...
    'color', Colors(5,:),'linewidth',1);                     %Mond @ TLI
plot(RSOIM*cups+rEM*cos(-domegt),RSOIM*sups+rEM*sin(-domegt),...
    'color', Colors(5,:),'linewidth',1);                     %SOI  @ TLI
plot(RM*cups+rEM,RM*sups,'color',Colors(5,:));               %Mond @ arrival
hp(2)=plot(RSOIM*cups+rEM,RSOIM*sups, ...
    'color', Colors(5,:),'linewidth',1);                     %SOI  @ arival
hp(3)=plot(rEM*cos(upsM),rEM*sin(upsM),'color',Colors(5,:)); %Mondbahn

line([0 rEM],[0 0],'linestyle',':');                %Verbindungslinie
line([0 r1*cos(phi1)],[0 r1*sin(phi1)]);            %Phasenwinkel @ SOI
line([0 rEM*cos(domegt)],[0 -rEM*sin(domegt)]);     %Position Mond @ TLI
line([0 5*r0*cos(upsP)],...
     [0 5*r0*sin(upsP)]);                           %Phasenwinkel @ TLI
hp(4)=plot(p0*cups./(1+exz0*cupsP),p0*sups./(1+exz0*cupsP), ...
     'color',Colors(2,:),'linewidth',2);            %Transferbahn 

alpha = ups2+pi+lambda1;
RZ = [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)];
vecH  = [a2B*(exz2-cosh(EH)); a2B*sqrt(exz2^2-1).*sinh(EH)];
vecHS = RZ*vecH;
hp(5) = plot(vecHS(1,:)+rEM , vecHS(2,:),...
     'color',Colors(9,:),'linewidth',2);            %"Mondumlaufbahn" 
axis([rEM-1*RSOIM, rEM+4*RSOIM, -2*RSOIM, 2*RSOIM]);
grid on
axis equal
hp2=legend(hp,'Erde','SOI Mond','Mondbahn','Transferbahn','Mondumlauf', ...
    'location','southwest'); 
set(hp2,'FontSize',14,'FontWeight','normal');
legend box off;
xlabel('x in m','FontSize',14); 
ylabel('y in m','FontSize',14)
grid on
set(gca,'FontSize',16);



