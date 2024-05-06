% -------------------------------------------------------------------------
% WurfmitohneReibung.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet Auf- und Abstiegszeiten senkrechter
% Wurf mit und ohne Reibung
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];
Marker = ['o','d','o','s','+'];

%%
% Import Data
% alle Daten in kg, m und s

v0     = 50;
g      = 9.81;
% Begin Rechnungen Quadratisch
eta    = 0.01;
vfQ    = sqrt(g/eta);
AB_up    = [0 v0];
options1  = odeset('AbsTol',1.e-9,'RelTol',1.e-5,'Events',@events1);
options2  = odeset('AbsTol',1.e-9,'RelTol',1.e-5,'Events',@events2);
tspan   = linspace(0,20,1000);
[tupQ,yupQ,te_upQ,ye_upQ,ie] = ode23(@(t,y)[y(2);-g-eta*y(2)^2],...
                     tspan,AB_up,options1);
AB_down = [ye_upQ(1) 0];
[tdownQ,ydownQ,te_downQ,ye_downQ,ie] = ode23(@(t,y)[y(2);-g-eta*y(2)*abs(y(2))],...
                      tspan,AB_down,options2);

                  
% Begin Rechnungen Linear
eta    = 0.1;
vf     = g/eta;
[tupL,yupL,te_upL,ye_upL,ie] = ode23(@(t,y)[y(2);-g-eta*y(2)],...
                     tspan,AB_up,options1);
AB_down = [ye_upL(1) 0];
[tdownL,ydownL,te_downL,ye_downL,ie] = ode23(@(t,y)[y(2);-g+eta*y(2)],...
                      tspan,AB_down,options2);
      
% Begin Rechnungen ohne Reibung
tmax  = v0/g;    
tup   = linspace(0,tmax,1000);
hmax  = v0^2/2/g; 
yup   = -g*tup.^2/2+v0*tup;
vup   = -g*tup+v0;
ydown = -g*tup.^2/2+hmax;
vdown = -g*tup;              


LgdStr (1,:) = sprintf('%s','quadratisch up  ');
LgdStr (2,:) = sprintf('%s','quadratisch down');
LgdStr (3,:) = sprintf('%s','linear up       ');
LgdStr (4,:) = sprintf('%s','linear down     ');
LgdStr (5,:) = sprintf('%s','ohne R up       ');
LgdStr (6,:) = sprintf('%s','ohne R down     ');

               
% Graphische Ausgabe 
figure();
hold on;
h(1)=plot(tupQ,yupQ(:,1), 'Color', Colors(2,:),...
              'LineWidth',2);
h(2)=plot(te_upQ(1)+tdownQ, ydownQ(:,1), 'Color', Colors(3,:),...
              'LineWidth',2);
line([te_upQ te_upQ], [0 max(yupQ(:,1))],'Color', Colors(4,:),...
              'LineWidth',1,'lineStyle',Style(3));
line([0 te_upQ+te_downQ],[max(yupQ(:,1)) max(yupQ(:,1))],...
          'Color', Colors(4,:),'LineWidth',1,'lineStyle',Style(3));
   
h(3)=plot(tupL,yupL(:,1), 'Color', Colors(2,:),...
              'LineWidth',2,'lineStyle',Style(2));
h(4)=plot(te_upL(1)+tdownL, ydownL(:,1), 'Color', Colors(3,:),...
              'LineWidth',2,'lineStyle',Style(2));
line([te_upL te_upL], [0 max(yupL(:,1))],'Color', Colors(4,:),...
              'LineWidth',1,'lineStyle',Style(3));
line([0 te_upL+te_downL],[max(yupL(:,1)) max(yupL(:,1))],...
          'Color', Colors(4,:),'LineWidth',1,'lineStyle',Style(3));
      
h(5)=plot(tup,yup, 'Color', Colors(2,:),...
              'LineWidth',1);
h(6)=plot(tmax+tup, abs(ydown), 'Color', Colors(3,:),...
              'LineWidth',1);
line([tmax tmax], [0 hmax],'Color', Colors(4,:),...
              'LineWidth',1,'lineStyle',Style(3));
line([0 2*tmax],[hmax hmax],...
          'Color', Colors(4,:),'LineWidth',1,'lineStyle',Style(3));
      
grid on
header2 = sprintf('z(t)'); 
ttl = title(header2);
ttl.FontSize = 16; ttl.FontWeight = 'normal'; 
ylabel('Höhe in m');
xlabel('Zeit in s');
lgd=legend(h,LgdStr(1:6,:),'Location',...
    'bestoutside','NumColumns',1);
lgd.FontSize=16; legend boxoff; 
set(gca, 'Fontsize', 14, 'linewidth', 1);



% Graphische Ausgabe 
figure();
hold on;
h(1)=plot(tupQ,yupQ(:,2), 'Color', Colors(2,:),...
              'LineWidth',2);
h(2)=plot(te_upQ(1)+tdownQ, abs(ydownQ(:,2)), 'Color', Colors(2,:),...
              'LineWidth',2);
line([te_upQ te_upQ], [0 max(yupQ(:,2))],'Color', Colors(4,:),...
              'LineWidth',1,'lineStyle',Style(3));
line([0 te_upQ+te_downQ],[vfQ vfQ],...
          'Color', Colors(4,:),'LineWidth',1,'lineStyle',Style(3));
h(3)=plot(tupL,yupL(:,2), 'Color', Colors(2,:),...
              'LineWidth',2,'lineStyle',Style(2));
h(4)=plot(te_upL(1)+tdownL, abs(ydownL(:,2)), 'Color', Colors(3,:),...
              'LineWidth',2,'lineStyle',Style(2));
line([te_upL te_upL], [0 max(yupL(:,2))],'Color', Colors(4,:),...
              'LineWidth',1,'lineStyle',Style(3));
line([0 te_upL+te_downL],[max(yupL(:,2)) max(yupL(:,2))],...
          'Color', Colors(4,:),'LineWidth',1,'lineStyle',Style(3));    
h(5)=plot(tup,vup, 'Color', Colors(2,:),...
              'LineWidth',1);
h(6)=plot(tmax+tup, abs(vdown), 'Color', Colors(3,:),...
              'LineWidth',1);
line([tmax tmax], [0 v0],'Color', Colors(4,:),...
              'LineWidth',1,'lineStyle',Style(3));
line([0 2*tmax],[v0 v0],...
          'Color', Colors(4,:),'LineWidth',1,'lineStyle',Style(3));
      
grid on
header2 = sprintf('z(t)'); 
ttl = title(header2);
ttl.FontSize = 16; ttl.FontWeight = 'normal'; 
ylabel('Geschwindigkeit in m/s');
xlabel('Zeit in s');
lgd=legend(h,LgdStr(1:6,:),'Location',...
    'bestoutside','NumColumns',1);
lgd.FontSize=16; legend boxoff; 
set(gca, 'Fontsize', 14, 'linewidth', 1);
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


 function [value,isterminal,direction] = events1(t,y)
        % Locate the time when velocity passes through zero 
        % and stop integration.  Here we use a nested function to avoid
        % passing the additional parameter P1 as an input argument.
        value = y(2);     % detect velocity = 0
        isterminal = 1;   % stop the integration
        direction = 0;    % negative direction
end

 function [value,isterminal,direction] = events2(t,y)
        % Locate the time when height passes through zero 
        % and stop integration.  Here we use a nested function to avoid
        % passing the additional parameter P1 as an input argument.
        value = y(1);     % detect height = 0
        isterminal = 1;   % stop the integration
        direction = -1;   % negative direction
 end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------
