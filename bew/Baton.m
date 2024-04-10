% -------------------------------------------------------------------------
% Baton.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Bewegungsgleichung Baton
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

% Konstanten:

   global g
   global l1
   global l2
   global l3
   global l4
   global m1
   global m2
g = 9.81;
l1 = input("l1 in m: ");
l2 = input("l2 in m: ");
l3 = input("l3 in m: ");
l4 = input("l4 in m: ");
m1 = input("m1 in kg: ");
m2 = input("m2 in kg: ");

% Anfangswinkel phi0:
phiadeg=input("phianfang in Grad: ");
global phi0
phi0=phiadeg/180*pi;

% Anfangswinkel theta0:
global theta0
theta0=phi0-pi/2;

% Anfangsiwnkel sigma0:
global sigma0
sigma0=pi - phi0;

% Anfangsbedingungen:
anf =[phi0, 0, theta0 ,0 ,sigma0, 0];

% Lösen der Differentialgleichung:
[t,L]=ode45(@f,[0 5],anf);

% Position des Projektils:
x1=-l1*sin(L(:,1))+l3*sin(L(:,1)-L(:,3));
y1=l1*cos(L(:,1))-l3*cos(L(:,3)-L(:,1));

% Position des GG:
x2=l2*sin(L(:,1))-l4*sin(L(:,1)+L(:,5));
y2=-l2*cos(L(:,1))+l4*cos(L(:,1)+L(:,5));

% Position des vorderen Ende des Arms:
x3=-l1*sin(L(:,1));
y3=l1*cos(L(:,1));

% Position des hinteren Ende des Arms:
x4=l2*sin(L(:,1));
y4=-l2*cos(L(:,1));

% Plot:
figure
subplot(2,1,1);
plot(t,L(:,1),"r-");
title("phi");
grid on

subplot(2,1,2);
plot(t,L(:,2),"r-");
title("phidot");
grid on

figure
plot(x1,y1,"-r", x2,y2,"-b", x3,y3,"-m", x4,y4,"-c", 0,0,"og");
title("???");
xlabel("x in m");
ylabel("y in m");
grid on
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------



% Differentialgleichung
function dy = f(t,y)
   global g
   global l1
   global l2
   global l3
   global l4
   global m1
   global m2
   
   dy = zeros(6,1);
   
   % ddotphi:
   dy(1) = y(2);
   dy(2) =  -((y(2)^2*l1^2*m1*cos(y(1))-g*l1*m1)*sin(y(1))*sin(y(1)-y(3))^2+((-y(2)^2*l1^2*m1*sin(y(1))^2+y(2)^2*l1^2*m1*cos(y(1))^2-g*l1*m1*cos(y(1)))*cos(y(1)-y(3))+(-y(4)^2+2*y(2)*y(4)-y(2)^2)*l1*l3*m1*cos(y(1)))*sin(y(1)-y(3))-y(2)^2*l1^2*m1*cos(y(1))*sin(y(1))*cos(y(1)-y(3))^2+(y(4)^2-2*y(2)*y(4)+y(2)^2)*l1*l3*m1*sin(y(1))*cos(y(1)-y(3))+(y(2)^2*l2^2*m2*cos(y(1))+g*l2*m2)*sin(y(1))*sin(y(5)+y(1))^2+((y(2)^2*l2^2*m2*sin(y(1))^2+y(2)^2*l2^2*m2*cos(y(1))^2+g*l2*m2*cos(y(1)))*cos(y(5)+y(1))+(-y(6)^2-2*y(2)*y(6)-y(2)^2)*l2*l4*m2*cos(y(1)))*sin(y(5)+y(1))+y(2)^2*l2^2*m2*cos(y(1))*sin(y(1))*cos(y(5)+y(1))^2+(y(6)^2+2*y(2)*y(6)-y(2)^2)*l2*l4*m2*sin(y(1))*cos(y(5)+y(1))+(g*l1*m1-g*l2*m2)*sin(y(1)))/(l1^2*m1*sin(y(1))^2*sin(y(1)-y(3))^2+2*l1^2*m1*cos(y(1))*sin(y(1))*cos(y(1)-y(3))*sin(y(1)-y(3))+l1^2*m1*cos(y(1))^2*cos(y(1)-y(3))^2+l2^2*m2*sin(y(1))^2*sin(y(5)+y(1))^2+2*l2^2*m2*cos(y(1))*sin(y(1))*cos(y(5)+y(1))*sin(y(5)+y(1))+l2^2*m2*cos(y(1))^2*cos(y(5)+y(1))^2-l2^2*m2-l1^2*m1) ;
   
   % ddottheta:
   dy(3) = y(4);
   dy(4) = -(((y(4)^2-2*y(2)*y(4)+2*y(2)^2)*l1^2*l3*m1*cos(y(1))-g*l1*l3*m1)*sin(y(1))*sin(y(1)-y(3))^2+(((-y(4)^2+2*y(2)*y(4)-2*y(2)^2)*l1^2*l3*m1*sin(y(1))^2+(y(4)^2-2*y(2)*y(4)+2*y(2)^2)*l1^2*l3*m1*cos(y(1))^2-g*l1*l3*m1*cos(y(1)))*cos(y(1)-y(3))+(-g*l2^2-g*l1*l2)*m2*sin(y(1))^2*sin(y(5)+y(1))^2+(((y(2)^2*l1*l2^2*m2*cos(y(1))^2+(-2*g*l2^2-g*l1*l2)*m2*cos(y(1)))*sin(y(1))-y(2)^2*l1*l2^2*m2*sin(y(1))^3)*cos(y(5)+y(1))+(y(6)^2+2*y(2)*y(6)+y(2)^2)*l1*l2*l4*m2*cos(y(1))*sin(y(1)))*sin(y(5)+y(1))+(-y(2)^2*l1*l2^2*m2*cos(y(1))*sin(y(1))^2+y(2)^2*l1*l2^2*m2*cos(y(1))^3-g*l2^2*m2*cos(y(1))^2)*cos(y(5)+y(1))^2+(-y(6)^2-2*y(2)*y(6)+y(2)^2)*l1*l2*l4*m2*sin(y(1))^2*cos(y(5)+y(1))+(g*l1*l2*m2-g*l1^2*m1)*sin(y(1))^2+(((-y(4)^2+2*y(2)*y(4)-y(2)^2)*l1*l3^2-y(2)^2*l1^3)*m1-y(2)^2*l1*l2^2*m2)*cos(y(1))+g*l2^2*m2+g*l1^2*m1)*sin(y(1)-y(3))+(-y(4)^2+2*y(2)*y(4)-2*y(2)^2)*l1^2*l3*m1*cos(y(1))*sin(y(1))*cos(y(1)-y(3))^2+(((-y(2)^2*l1*l2^2*m2*cos(y(1))^2-g*l1*l2*m2*cos(y(1)))*sin(y(1))-y(2)^2*l1*l2^2*m2*sin(y(1))^3)*sin(y(5)+y(1))^2+((-3*y(2)^2*l1*l2^2*m2*cos(y(1))*sin(y(1))^2-y(2)^2*l1*l2^2*m2*cos(y(1))^3-g*l1*l2*m2*cos(y(1))^2)*cos(y(5)+y(1))+(y(6)^2+2*y(2)*y(6)+y(2)^2)*l1*l2*l4*m2*cos(y(1))^2)*sin(y(5)+y(1))-2*y(2)^2*l1*l2^2*m2*cos(y(1))^2*sin(y(1))*cos(y(5)+y(1))^2+(-y(6)^2-2*y(2)*y(6)+y(2)^2)*l1*l2*l4*m2*cos(y(1))*sin(y(1))*cos(y(5)+y(1))+((g*l1*l2*m2-g*l1^2*m1)*cos(y(1))+y(2)^2*l1*l2^2*m2+((y(4)^2-2*y(2)*y(4)+y(2)^2)*l1*l3^2+y(2)^2*l1^3)*m1)*sin(y(1)))*cos(y(1)-y(3))+(y(2)^2*l2^2*l3*m2*cos(y(1))+g*l2*l3*m2)*sin(y(1))*sin(y(5)+y(1))^2+((y(2)^2*l2^2*l3*m2*sin(y(1))^2+y(2)^2*l2^2*l3*m2*cos(y(1))^2+g*l2*l3*m2*cos(y(1)))*cos(y(5)+y(1))+(-y(6)^2-2*y(2)*y(6)-y(2)^2)*l2*l3*l4*m2*cos(y(1)))*sin(y(5)+y(1))+y(2)^2*l2^2*l3*m2*cos(y(1))*sin(y(1))*cos(y(5)+y(1))^2+(y(6)^2+2*y(2)*y(6)-y(2)^2)*l2*l3*l4*m2*sin(y(1))*cos(y(5)+y(1))+(g*l1*l3*m1-g*l2*l3*m2)*sin(y(1)))/(l1^2*l3*m1*sin(y(1))^2*sin(y(1)-y(3))^2+2*l1^2*l3*m1*cos(y(1))*sin(y(1))*cos(y(1)-y(3))*sin(y(1)-y(3))+l1^2*l3*m1*cos(y(1))^2*cos(y(1)-y(3))^2+l2^2*l3*m2*sin(y(1))^2*sin(y(5)+y(1))^2+2*l2^2*l3*m2*cos(y(1))*sin(y(1))*cos(y(5)+y(1))*sin(y(5)+y(1))+l2^2*l3*m2*cos(y(1))^2*cos(y(5)+y(1))^2-l2^2*l3*m2-l1^2*l3*m1) ;
   
   % ddotsigma:
   dy(5) = y(6);
   dy(6) = (((g*l1*l2+g*l1^2)*m1*sin(y(1))^2*sin(y(5)+y(1))+(y(2)^2*l1^2*l2*m1*sin(y(1))^3+(g*l1*l2*m1*cos(y(1))-y(2)^2*l1^2*l2*m1*cos(y(1))^2)*sin(y(1)))*cos(y(5)+y(1))+(y(2)^2*l1^2*l4*m1*cos(y(1))-g*l1*l4*m1)*sin(y(1)))*sin(y(1)-y(3))^2+(((y(2)^2*l1^2*l2*m1*sin(y(1))^3+(y(2)^2*l1^2*l2*m1*cos(y(1))^2+(g*l1*l2+2*g*l1^2)*m1*cos(y(1)))*sin(y(1)))*sin(y(5)+y(1))+(3*y(2)^2*l1^2*l2*m1*cos(y(1))*sin(y(1))^2-y(2)^2*l1^2*l2*m1*cos(y(1))^3+g*l1*l2*m1*cos(y(1))^2)*cos(y(5)+y(1))-y(2)^2*l1^2*l4*m1*sin(y(1))^2+y(2)^2*l1^2*l4*m1*cos(y(1))^2-g*l1*l4*m1*cos(y(1)))*cos(y(1)-y(3))+(y(4)^2-2*y(2)*y(4)+y(2)^2)*l1*l2*l3*m1*cos(y(1))*sin(y(1))*sin(y(5)+y(1))+(y(4)^2-2*y(2)*y(4)+y(2)^2)*l1*l2*l3*m1*cos(y(1))^2*cos(y(5)+y(1))+(-y(4)^2+2*y(2)*y(4)-y(2)^2)*l1*l3*l4*m1*cos(y(1)))*sin(y(1)-y(3))+((y(2)^2*l1^2*l2*m1*cos(y(1))*sin(y(1))^2+y(2)^2*l1^2*l2*m1*cos(y(1))^3+g*l1^2*m1*cos(y(1))^2)*sin(y(5)+y(1))+2*y(2)^2*l1^2*l2*m1*cos(y(1))^2*sin(y(1))*cos(y(5)+y(1))-y(2)^2*l1^2*l4*m1*cos(y(1))*sin(y(1)))*cos(y(1)-y(3))^2+((-y(4)^2+2*y(2)*y(4)-y(2)^2)*l1*l2*l3*m1*sin(y(1))^2*sin(y(5)+y(1))+(-y(4)^2+2*y(2)*y(4)-y(2)^2)*l1*l2*l3*m1*cos(y(1))*sin(y(1))*cos(y(5)+y(1))+(y(4)^2-2*y(2)*y(4)+y(2)^2)*l1*l3*l4*m1*sin(y(1)))*cos(y(1)-y(3))+((y(6)^2+2*y(2)*y(6)+2*y(2)^2)*l2^2*l4*m2*cos(y(1))+g*l2*l4*m2)*sin(y(1))*sin(y(5)+y(1))^2+(((-y(6)^2-2*y(2)*y(6)+2*y(2)^2)*l2^2*l4*m2*sin(y(1))^2+(y(6)^2+2*y(2)*y(6)+2*y(2)^2)*l2^2*l4*m2*cos(y(1))^2+g*l2*l4*m2*cos(y(1)))*cos(y(5)+y(1))+(g*l2^2*m2-g*l1*l2*m1)*sin(y(1))^2+(((-y(6)^2-2*y(2)*y(6)-y(2)^2)*l2*l4^2-y(2)^2*l2^3)*m2-y(2)^2*l1^2*l2*m1)*cos(y(1))-g*l2^2*m2-g*l1^2*m1)*sin(y(5)+y(1))+(-y(6)^2-2*y(2)*y(6)+2*y(2)^2)*l2^2*l4*m2*cos(y(1))*sin(y(1))*cos(y(5)+y(1))^2+((g*l2^2*m2-g*l1*l2*m1)*cos(y(1))+((y(6)^2+2*y(2)*y(6)-y(2)^2)*l2*l4^2-y(2)^2*l2^3)*m2-y(2)^2*l1^2*l2*m1)*sin(y(1))*cos(y(5)+y(1))+(g*l1*l4*m1-g*l2*l4*m2)*sin(y(1)))/(l1^2*l4*m1*sin(y(1))^2*sin(y(1)-y(3))^2+2*l1^2*l4*m1*cos(y(1))*sin(y(1))*cos(y(1)-y(3))*sin(y(1)-y(3))+l1^2*l4*m1*cos(y(1))^2*cos(y(1)-y(3))^2+l2^2*l4*m2*sin(y(1))^2*sin(y(5)+y(1))^2+2*l2^2*l4*m2*cos(y(1))*sin(y(1))*cos(y(5)+y(1))*sin(y(5)+y(1))+l2^2*l4*m2*cos(y(1))^2*cos(y(5)+y(1))^2-l2^2*l4*m2-l1^2*l4*m1) ;
end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------
