% -------------------------------------------------------------------------
% Doppelpendel03.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Doppelpendel - Kleinwinkelnäherung
% 
% Programm berechnet die zeitliche Entwicklung des Doppelpendels für
% verschiedene Parameter (ohne Reibung) in der Kleinwinkelnäherung. 
%
%--------------------------------------------------------------------------
% Benutzt Teile des Programms von Javier Hasbun
% Copyright (c) 2008, Javier Hasbun
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
% * Redistributions of source code must retain the above copyright
%   notice, this list of conditions and the following disclaimer.
% * Redistributions in binary form must reproduce the above copyright
%   notice, this list of conditions and the following disclaimer in
%   the documentation and/or other materials provided with the distribution
% * Neither the name of the  nor the names
%   of its contributors may be used to endorse or promote products derived
%   from this software without specific prior written permission.
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

% Parameter
L1=1; L2=1; m1=1; m2=1; g=9.81;      
tau0=sqrt(g/(L1+L2));               % Zeitspanne
tmax=5*tau0; tstep=0.01;             
tv=(0.0:tstep:tmax);N=length(tv);   % Simulationszeit, Zeit-Vektor

% Anfangswerte
th10 =0.10;  th20 =0.0;             % AB Winkel in rad
th10d=0.10;  th20d=0.0;             % AB Geschwindigkeit in rad/s
AB=[th10;th10d;th20;th20d];         % AB Vektor

%%
% Runge-Kutta-Verfahren-MATLAB ODE45
opt=odeset('AbsTol',1.e-9,'RelTol',1.e-5);           
[t,w]=ode45(@(t,w)DGL(t,w,L1,L2,m1,m2,g),tv,AB,opt); 

% Winkel über die Zeit
% w(1)=theta1, w(2)=theta1_dot, w(3)=theta2, w(4)=theta2_dot
str1=cat(2,'Doppelpendel: \theta_{10}=',num2str(th10,3),...
        ', \theta_{20}=',num2str(th20,3),' rad, (d\theta_{1}/dt)_0=',...
        num2str(th10d,3),', (d\theta_{2}/dt)_0=',num2str(th20d,3), 'rad/s');
str2=cat(2,'L_1=',num2str(L1,3),', L_2=',num2str(L2,3),' m, ');
str3=cat(2,'m_1=',num2str(m1,3),', m_2=',num2str(m2,3),' kg');
str4=cat(2,str2,str3);

% Trajektorien graphisch
figure();

subplot(2,2,1)
plot(t,w(:,1),'color', Colors(2,:),'LineWidth',1,'LineStyle',Style(3));
hold on;
plot(t,w(:,3),'color', Colors(4,:),'LineWidth',1,'LineStyle',Style(5));
xlabel('\it{t} \rm','FontSize',14);
ylabel('\theta_1, \theta_2','FontSize',14)
%title(str1,'FontSize',12)
ym=min([w(:,1);w(:,3)]); yp=max([w(:,1);w(:,3)]);    
axis([0,tmax,ym*(1+0.5),yp*(1+0.5)])
% text(.1,yp*(1-0.1),str2)
% text(.1,ym*(1-0.1),str3)
h=legend('\theta_1','\theta_2'); set(h,'FontSize',14)
legend box off
grid on;

subplot(2,2,3), plot(w(:,1),w(:,2),'color', Colors(2,:),'LineWidth',2)
xlabel ('\theta_1','FontSize',14);
ylabel('d\theta_1/d\it{t}','FontSize',14)
ym=min(w(:,2)); yp=max(w(:,2));    
xm=min(w(:,1)); xp=max(w(:,1));    
axis([xm*(1+0.5),xp*(1+0.5),ym*(1+0.5),yp*(1+0.5)])

grid on;

% subplot(2,1,2), plot(w(:,3),w(:,4),'color', Colors(4,:),'LineWidth',2)
% xlabel ('\theta_2','FontSize',14),ylabel('d\theta_2/d\it{t}','FontSize',14)
% grid on;

%%
% Kleinwinkelnäherung analytisch

im=sqrt(-1);
th10  = 0.1;
th10d = 0.0;
omeg1 = sqrt(2+sqrt(2))*sqrt(g/L1);
omeg2 = sqrt(2-sqrt(2))*sqrt(g/L1);
M=[  1      1     1      1;
    -1     -1     1      1;
    +im*omeg1 -im*omeg1 +im*omeg2 -im*omeg2;
    -im*omeg1 +im*omeg1 +im*omeg2 -im*omeg2];
AB=[th10 0 th10d 0]';
Coef = M\AB;
A = Coef(1);
B = Coef(2);
C = Coef(3);
D = Coef(4);

lsg1 = A*exp(im*omeg1*t) + B*exp(-im*omeg1*t) + ...
       C*exp(im*omeg2*t) + D*exp(-im*omeg2*t);
lsg2 = -sqrt(2)*(A*exp(im*omeg1*t) + B*exp(-im*omeg1*t)) + ...
       sqrt(2)*(C*exp(im*omeg2*t) + D*exp(-im*omeg2*t));
dlsg1= im*omeg1*(A*exp(im*omeg1*t) - B*exp(-im*omeg1*t))+ ...
       im*omeg2*sqrt(2)*(C*exp(im*omeg2*t) - D*exp(-im*omeg2*t));
      

subplot(2,2,2)
plot(tv,lsg1,'color', Colors(2,:),'LineWidth',1,'LineStyle',Style(3));
hold on;
plot(tv,lsg2,'color', Colors(4,:),'LineWidth',1,'LineStyle',Style(5));
xlabel('\it{t} \rm','FontSize',14);
ylabel('\theta_1, \theta_2','FontSize',14)
% title(str1,'FontSize',12)
ym2=min([w(:,1);w(:,3)]); yp2=max([w(:,1);w(:,3)]);    
axis([0,tmax,ym2*(1+0.5),yp2*(1+0.5)])
text(.1,yp*(1-0.1),str2)
% text(.1,ym*(1-0.1),str3)
h=legend('\theta_1','\theta_2'); set(h,'FontSize',14)
legend box off
grid on;

subplot(2,2,4);
plot(lsg1,dlsg1,'color', Colors(2,:),'LineWidth',2)
xlabel ('\theta_1','FontSize',14);
ylabel('d\theta_1/d\it{t}','FontSize',14)
axis([xm*(1+0.5),xp*(1+0.5),ym*(1+0.5),yp*(1+0.5)])
grid on;
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------

% DGL System
function dYdt = DGL(t,w,L1,L2,m1,m2,g)
    % w(1):th1, w(2):th1_dot, w(3):th2, w(4):th2_dot
    tp = w(1)-w(3); cs=sin(tp); cc=cos(tp);
    ta = m2*L2*cc/(m1+m2)/L1;
    tb = (m2*L2*w(4).^2.*cs+(m1+m2)*g*sin(w(1)))/(m1+m2)/L1;
    tc = (L1*tb.*cc+L1*w(2).^2.*cs-g*sin(w(3)))./(L2-ta.*L1.*cc);
    dYdt = [w(2);-ta.*tc-tb;w(4);tc];
end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------

