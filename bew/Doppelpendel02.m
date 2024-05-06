% -------------------------------------------------------------------------
% Doppelpendel02.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet die zeitliche Entwicklung des Doppelpendels für
% verschiedene Parameter (mit/ohne Reibung). 
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
L1=1; L2=1; m1=4; m2=1; g=9.81;      
tau0=sqrt(g/(L1+L2));               %Zeitspanne
tmax=10*tau0; tstep=0.05;             
tv=(0.0:tstep:tmax);N=length(tv);   %Simulationszeit, Zeit-Vektor

% Anfangswerte
th10=3; th20=2;               %AB Winkel in rad
th10d=0.0; th20d=0.0;               %AB Geschwindigkeit in rad/s
AB=[th10;th10d;th20;th20d];         %AB Vektor

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

%%
% Trajektorien
figure();
% subplot(3,1,1)
% plot(t,w(:,1),'color', Colors(2,:),'LineWidth',1,'LineStyle',Style(3));
% hold on;
% plot(t,w(:,3),'color', Colors(4,:),'LineWidth',1,'LineStyle',Style(5));
% xlabel ('\it{t} \rm','FontSize',14), ylabel('\theta_1, \theta_2','FontSize',14)
% title(str1,'FontSize',12)
% ym=min([w(:,1);w(:,3)]); yp=max([w(:,1);w(:,3)]);    
% axis([0,tmax,ym*(1+0.5),yp*(1+0.5)])
% text(.1,yp*(1-0.1),str2)
% % text(.1,ym*(1-0.1),str3)
% h=legend('\theta_1','\theta_2'); set(h,'FontSize',14)
% legend box off

subplot(2,1,1), plot(w(:,1),w(:,2),'color', Colors(2,:),'LineWidth',2)
xlabel ('\theta_1','FontSize',14),ylabel('d\theta_1/d\it{t}','FontSize',14)
grid on;

subplot(2,1,2), plot(w(:,3),w(:,4),'color', Colors(4,:),'LineWidth',2)
xlabel ('\theta_2','FontSize',14),ylabel('d\theta_2/d\it{t}','FontSize',14)
grid on;

%%
% Simulation  
x1=L1*sin(w(:,1)); 
y1=L1*cos(w(:,1));
x2=x1+L2*sin(w(:,3)); 
y2=y1+L2*cos(w(:,3));
v=L1+L2; y1=v-y1; y2=v-y2;  
figure();
vx=max([v;x1;x2]);vy=max([v;y1;y2]);
axis([-vx,vx,0,vy]) 
 for ki=1:N
     clf
     axis([-vx,vx,0,vy])
     hold on
     strtime = strjoin({' t = ', num2str(t(ki), '% 6.2f s')});
     text(-vx+0.15, vy - 0.15,strtime);

     plot(0,v,'s',...
          'MarkerSize',7,...
          'MarkerEdgeColor','k',...
          'MarkerFaceColor','k');              %Aufhängepunkt
     line([0,x1(ki)],[v,y1(ki)],'color', Colors(2,:),'LineWidth',2); %Arm1
     h(1)=plot(x1(ki),y1(ki),'o',...
          'MarkerSize',6,...
          'MarkerEdgeColor',Colors(2,:),...
          'MarkerFaceColor',Colors(2,:)'); %m1 @ x1,y1
     line([x1(ki),x2(ki)],[y1(ki),y2(ki)],'color',Colors(4,:),...
          'LineWidth',2);  %Arm2 
     h(2)=plot(x2(ki),y2(ki),'o',...
          'MarkerSize',6,...
          'MarkerEdgeColor',Colors(4,:),...
          'MarkerFaceColor',Colors(4,:)); %m2 @ x2,y2
     pause(.025)
 end
 
h(3)=plot(x1,y1,'color', Colors(2,:),'LineWidth',1,'LineStyle',Style(3));                               
                                                             %m1 Spur
h(4)=plot(x2,y2,'color', Colors(4,:),'LineWidth',1,'LineStyle',Style(5));
                                                             %m2 Spur
h=legend(h,'m_1','m_2','m_1 Trajektorie', 'm_2 Trajektorie');
set(h,'FontSize',13)
xlabel('x_1, x_2','FontSize',13), ylabel('y_1, y_2','FontSize',13)
title(str1,'FontSize',12), 
text(-vx+0.15,0 + 0.25, str4)
legend box off
grid on;
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


% DGL System
function dYdt = DGL(t,w,L1,L2,m1,m2,g)
    %w(1):theta1, w(2):theta1_dot, w(3):theta2, w(4):theta2_dot
    tp = w(1)-w(3); cs=sin(tp); cc=cos(tp);
    ta = m2*L2*cc/(m1+m2)/L1;
    tb = (m2*L2*w(4).^2.*cs+(m1+m2)*g*sin(w(1)))/(m1+m2)/L1;
    tc = (L1*tb.*cc+L1*w(2).^2.*cs-g*sin(w(3)))./(L2-ta.*L1.*cc);
    dYdt = [w(2);-ta.*tc-tb;w(4);tc];
end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------

