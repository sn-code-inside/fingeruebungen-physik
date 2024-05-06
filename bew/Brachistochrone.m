% -------------------------------------------------------------------------
% Brachistochrone.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Berechnungen zur Brachistochrone und Tautochrone, sowie
% Lagrange-Funktion und % Lagrange Gleichungen des Zykloidenpendels 
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

% Definiert die verallgemeinerten Koordinaten und Parameter

syms q dq u du phi dphi A m dbeta g  x z 'real'


% %% Symbolische Lösung Lagrange Funktion
% 
% % Verallgemeinerte Koordinaten und Ableitungen
% q  = u;
% dq = du;
% % kinetische Energien     
% T_kin = m/2 *4*A*A*4*(du^2)
% % potentielle Energien
% U_pot = -m*g*2*A*(1-u^2);
% %Lagrange Funktion
% L = T_kin - U_pot;
% 
% fprintf("------------------------------------\nLagrange-Funktion L:");
% fprintf("\n------------------------------------\n");
% fprintf(1, 'T = %s \n', T_kin);
% fprintf(1, 'U = %s \n', U_pot);
% fprintf(1, 'L = T - U = %s \n', L);
% fprintf("\n");
% 
% %Externe Kräfte 
% Q = 0;
% 
% % Berechnet Euler-Lagrange-Gleichungen:
% EQ = EulerLagrange(q,dq,L,Q,true);
% 

%% Numerische Simulation verschiedener Kurven

zA = 0;
xA = 0;
zE = -10;
xE = 25;

g = 9.81;
% m=1;

% Freier Fall
TE0 = sqrt(2*abs(zE)/g);

% Kurve 1 : Gerade
beta = (zE/xE);
func = @(u, beta) sqrt(-(1+beta^2)./(beta*u));
TE1 = sqrt(1/2/g)*integral(@(u) func(u,beta),0,xE);
x1 = linspace(0,xE,100);
z1 = beta*x1;

% Kurve 2 : Parabel
beta2 = -(zE/xE/xE);
func = @(u, xE, zE, beta2) sqrt(-(1+4*beta2^2*(u-xE).^2)...
         ./(zE+beta2*(u-xE).^2));
TE2 = sqrt(1/2/g)*integral(@(u) func(u, xE, zE, beta2),0,xE);
x2 = x1;
z2 = beta2*(x2-xE).^2+zE;

% Kurve 3 : Knickfunktion
TE3 = TE0 + xE/sqrt(2*abs(zE)*g);

% Kurve 4 : Kreis R = s
s       = sqrt(xE^2+zE^2);
R       = s;
myfun = @(u,s,xE,zE) s^2-(xE-u)^2-(zE-sqrt(s^2-u^2))^2;  
fun = @(u) myfun(u, s, xE, zE);    % function of x alone
xM = fzero(fun,0.1);
zM = sqrt(s^2-xM^2);
func = @(u, xM, zM, R) sqrt(-(1+((u-xM).^2./(R^2-(u-xM).^2)))...
         ./(zM-sqrt(R^2-(u-xM).^2)));
TE4 = sqrt(1/2/g)*integral(@(u) func(u, xM, zM, R),0,xE);
x4 = x1;
z4 = -sqrt(R^2-(x4-xM).^2) + zM;

% Kurve 5 : Zykloide
s1      = xE/zE;
myfun = @(u,s) s1*(1-cos(u))-u+sin(u);  % parameterized function
fun   = @(u) myfun(u, s1);              % function of x alone
phiE  = fzero(fun,-2*pi);
A     = zE/(1-cos(phiE));
phi  = linspace(0,phiE,100);
x5   = A*(phi-sin(phi));
z5   = A*(1-cos(phi));
TE5 = -sqrt(abs(A)/g)*phiE;


fprintf('\nTE0 : %5.4f s (Freier Fall)\n', TE0);
fprintf('TE1 : %5.4f s (Gerade)\n', TE1);
fprintf('TE2 : %5.4f s (Parabel)\n', TE2);
fprintf('TE3 : %5.4f s (Knickfunktion)\n', TE3);
fprintf('TE4 : %5.4f s (Kreis)\n', TE4);
fprintf('TE5 : %5.4f s (Zykloide)\n', TE5);

figure();
hp(1)=plot(x1,z1,'Color',Colors(1,:),'LineWidth',1,'LineStyle',Style(1));
hold on
hp(2)=plot(x2,z2,'Color',Colors(2,:),'LineWidth',1,'LineStyle',Style(2));
hp(3)=line([0 0],[0 zE],'Color',Colors(3,:),'LineWidth',1,'LineStyle',Style(1));
hp(3)=line([0 xE],[zE zE],'Color',Colors(3,:),'LineWidth',1,'LineStyle',Style(1));
hp(4)=plot(x4,z4,'Color',Colors(4,:),'LineWidth',1,'LineStyle',Style(3));
hp(5)=plot(x5,z5,'Color',Colors(5,:),'LineWidth',1,'LineStyle',Style(1));
h=legend(hp(1:5),strcat(num2str(TE1,4),' s - Gerade'),...
                 strcat(num2str(TE2,4),' s - Parabel'),...
                 strcat(num2str(TE3,4),' s - Knick'),...
                 strcat(num2str(TE4,4),' s - Kreis'),...
                 strcat(num2str(TE5,4),' s - Zykloide'),'NumColumns',1);
set(h,'FontSize',14)
axis([xA xE 1.2*zE zA]);
xlabel('\it{x} \rm in m','FontSize',14)
ylabel('\it{z} \rm in m','FontSize',14)
set(gca,'FontSize',14)

grid on;
legend box off;

%%
% Darstellung Zykloide


phi = linspace(0,4*pi,1000);
R   = 1;
A   = R;

x   = A*(phi-sin(phi));
y   = A*(1-cos(phi));

figure()
plot(x, y,'Color',Colors(4,:),'LineWidth',2,'LineStyle',Style(1));
axis equal;
axis([-R, 4*pi*R+R, 0, 4*R]);
hold on;
for k = 1:8
  gamma(k) = (k-1)*pi/3;  
  xM = R*gamma(k);
  xi = linspace(0,2*pi,100);
  hp(k)= plot(R*cos(xi)+xM, R*sin(xi)+R,...
      'Color', Colors(mod(k,6)+1,:),'LineWidth',1);
  lgdstr(k) = string(strcat(num2str(rad2deg(gamma(k)),3),'°'));
end
grid on
h=legend(hp, lgdstr,'NumColumns',4,'location','north');
legend box off
set(h,'FontSize',16)
xlabel('\it{x} ','FontSize',14)
ylabel('\it{y} ','FontSize',14)
set(gca,'FontSize',16)
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------









