% -------------------------------------------------------------------------
% Waschmaschine02.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet die zeitliche Entwicklung 
% der Vibrationen in der Auslaufphase
% einer Trommelwaschmaschine. 
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

% Parameter
R = 0.5; 
L2 = 1; 
M = 10; 
m = 1; 
g = 9.81; 
k = 4000;     
b= m*R/(M+m); 
Omega  = 75;
omega0 = sqrt(k/(M+m))*2;
stromega0 = strcat('\omega_0 =  ',num2str(omega0,' % 4.1f'),' 1/s');
eta    = 0.25*omega0;
streta = strcat('\eta =  ',num2str(eta,' % 4.1f'),' 1/s');
beta    = 0.01*omega0;
strbeta = strcat('\beta =  ',num2str(beta,' % 4.1f'),' 1/s');
tmax   = 5;               % Zeitspanne
tstep  = 0.002;             
tv=(0.0:tstep:tmax); N=length(tv);   % Simulationszeit, Zeit-Vektor
str2=cat(2,'\omega_0 = ',num2str(omega0,' % 4.1f'),' 1/s ' ,' \Omega = ',num2str(Omega,' % 4.1f'),' 1/s ');
str3=cat(2,'\eta   = ',num2str(eta,' % 4.1f'),' 1/s ',' \beta  = ',num2str(beta,' % 4.1f'),' 1/s');
str4=cat(2,str2,str3);

% Anfangswerte
r_stat = b*Omega^2/sqrt((omega0^2-Omega^2)^2+eta^2*Omega^2);
phi0  = atan(eta*Omega/(omega0^2-Omega^2));
dphi0 = Omega;
x0   = -r_stat*sin(phi0);            dx0   = Omega*r_stat*cos(phi0); 
z0   = +r_stat*cos(phi0)-g/omega0^2; dz0=Omega*r_stat*sin(phi0);     
AB=[x0;dx0;z0;dz0;phi0;dphi0];                 % AB Vektor


P = m*(M+m)*eta/(M*R*(M+2+m));
Q = m*(M+m)*omega0^2/(M*R*(M+2+m));
ddphi0 = 0;

% Runge-Kutta-Verfahren-MATLAB ODE45
opt=odeset('AbsTol',1.e-9,'RelTol',1.e-5);           
[t,w]=ode45(@(t,w)DGL(t,w,eta,beta,b,g,Q,P,omega0),tv,AB,opt); 


%%
% Position über die Zeit - Variation Eigenfrequenz
% w(1)=x, w(2)=dx, w(3)=z, w(4)=dz
figure()
omega0 = sqrt(k/(M+m))*1;
Q = m*(M+m)*omega0^2/(M*R*(M+2+m));
stromega1 = strcat('\omega_0 =  ',num2str(omega0,' % 4.1f'),' 1/s');
[t,w1]=ode45(@(t,w)DGL(t,w,eta,beta,b,g,Q,P,omega0),tv,AB,opt); 
omega0 = sqrt(k/(M+m))*3;
Q = m*(M+m)*omega0^2/(M*R*(M+2+m));
stromega2 = strcat('\omega_0 =  ',num2str(omega0,' % 4.1f'),' 1/s');
[t,w2]=ode45(@(t,w)DGL(t,w,eta,beta,b,g,Q,P,omega0),tv,AB,opt); 
x = w(:,1); z = w(:,3); r(1,:) = sqrt(x.^2 + z.^2);
x = w1(:,1); z = w1(:,3); r(2,:) = sqrt(x.^2 + z.^2);
x = w2(:,1); z = w2(:,3); r(3,:) = sqrt(x.^2 + z.^2);
subplot(2,1,1)
hold on
for ki = 1:3 
    h(ki)=plot(t, r(ki,:),'Color', Colors(ki+1,:));
end
xlabel ('t in s','FontSize',14),ylabel('r in m ','FontSize',14)
str1='Position Trommel - Variation Eigenfrequenz';
title(str1,'FontSize',12)
grid on;
hl=legend(h, stromega0, stromega1, stromega2,...
         'location','northwest'); 
set(hl,'FontSize',14)
legend box off

%%
% Position über die Zeit - Variation Dämpfung
% w(1)=x, w(2)=dx, w(3)=z, w(4)=dz
omega0 = sqrt(k/(M+m))*2;
Q = m*(M+m)*omega0^2/(M*R*(M+2+m));
beta    = 0.002*omega0; 
P = m*(M+m)*eta/(M*R*(M+2+m));
strbeta1 = strcat('\beta =  ',num2str(beta,' % 4.1f'),' 1/s');
[t,w1]=ode45(@(t,w)DGL(t,w,eta,beta,b,g,Q,P,omega0),tv,AB,opt); 
beta    = 0.05*omega0;
P = m*(M+m)*eta/(M*R*(M+2+m));
strbeta2 = strcat('\beta =  ',num2str(beta,' % 4.1f'),' 1/s');
[t,w2]=ode45(@(t,w)DGL(t,w,eta,beta,b,g,Q,P,omega0),tv,AB,opt); 
x = w(:,1); z = w(:,3); r(1,:) = sqrt(x.^2 + z.^2);
x = w1(:,1); z = w1(:,3); r(2,:) = sqrt(x.^2 + z.^2);
x = w2(:,1); z = w2(:,3); r(3,:) = sqrt(x.^2 + z.^2);
subplot(2,1,2)
hold on
for ki = 1:3 
    h(ki)=plot(t, r(ki,:),'Color', Colors(ki+1,:));
end
xlabel ('t in s','FontSize',14),ylabel('r in m ','FontSize',14)
str1='Position Trommel - Variation Dämpfung';
title(str1,'FontSize',12)
grid on;
hl=legend(h, strbeta, strbeta1, strbeta2); 
set(hl,'FontSize',14)
legend box off


%%
% Trajektorien

figure()
x = w(:,1);
z = w(:,3);
axis([-0.25,0.25,-0.25, 0.25])
for ki=1:5:N
     hold on
     strtime = strjoin({' t = ', num2str(t(ki), '% 6.2f s')});
     hs(1)=text(0.1, -0.2 ,strtime);
     hs(2)=plot(x(ki),z(ki),'o',...
          'MarkerSize',6,...
          'MarkerEdgeColor',Colors(2,:),...
          'MarkerFaceColor',Colors(2,:)'); %@ Trommelmittelpunkt x,z
     hs(3)=plot(x(ki),z(ki),'+',...
          'MarkerSize',2,...
          'MarkerEdgeColor',Colors(2,:),...
          'MarkerFaceColor',Colors(2,:)'); %@ Trommelmittelpunkt x,z
     pause(.01)
     hs(1).Visible = 'off';
     hs(2).Visible = 'off';
 end
text(0.1, -0.2 ,strtime);
plot(w(:,1), w(:,3),'color', Colors(2,:),'LineWidth',1,'LineStyle',Style(3));
hold on;
xlabel ('x in m ','FontSize',14), ylabel('z in m','FontSize',14)
title(str1,'FontSize',12)
axis([-0.25,0.25,-0.25, 0.25])
yp = max(w(:,3));
text(-0.20,yp*(1+0.1),str2)
text(-0.20,yp*(1+0.25),str3)
grid on
axis square
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


% DGL System
function dYdt = DGL(t,w,eta,beta,b,g,Q,P,omega0)
    % w(1):x
    % w(2):dx
    % w(3):z
    % w(4):dz
    
    x    = w(1);
    dx   = w(2);
    z    = w(3);
    dz   = w(4);
    phi  = w(5);
    dphi = w(6);
    ddphi = -beta*dphi+beta*dphi+Q*(x*cos(phi)+z*sin(phi))+ ...
                             P*(dx*cos(phi)+z*sin(phi));
    dYdt = [dx;...
       -eta*dx-omega0^2*x-b*(ddphi*cos(phi)-dphi^2*sin(phi));...
       dz;...
       -eta*dz-omega0^2*z-g-b*(ddphi*sin(phi)+dphi^2*cos(phi));...
       dphi;
       -beta*dphi+Q*(x*cos(phi)+z*sin(phi))+P*(dx*cos(phi)+z*sin(phi))];
end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------

