% -------------------------------------------------------------------------
% HarmonOszillator.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Beispielberechnungen zum Harmonischen Oszillator
% a) verschiedene Dämpfungen (Reibungskoeffizienten)
% b) verschiedene Anregungsfrequenzen und Dämpfungen(Reibungskoeffizienten)
% c) Frequenz-, Phasengang und Energieübertrag bei Anregung
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

% Parameter
m     = 0.2;                             %Masse
k     = 1;                               %rücktreibende Kraftkonstante
q0    = 1.0;                             %Auslemkung t = 0 
qdot0 = 5.0;                             %Anfangsgeschwindigkeit
omega0= sqrt(k/m);                       %Eigenfrequenz

tmax  = 8;                               %Maximalzeit
NPkt  = 500;                             %Anzahl Punkte
t     = linspace(0,tmax,NPkt);           %Zeitbereich


%%
%Teil a) 

%% 
% Fall ohne Dämpfung
eta10 = 0.00;                            %Reibungsterm 1
B     = sqrt(q0^2+qdot0^2/omega0^2);     %Konstante B  
theta = atan(omega0*q0/(qdot0));         %Winkel theta
q     = B.*sin(omega0*t+theta);          %Ungedämpfter Oszillator

%% 
% Fall mit schwacher Dämpfung
eta11  = 0.05;                            %Reibungsterm 1
gamma = eta11/2/m;                        %Kritischer Parameter gamma
Disc  = omega0^2-gamma^2;                 %Diskriminator: hier positiv
if Disc <= 0                         
    disp('keine schwache Dämpfung'); 
    return; 
end
omega = sqrt(Disc);                               %Aktuelle Frequenz
B     = sqrt(q0^2+(qdot0+gamma*q0)^2/omega^2);    %Konstante B 
theta = atan(omega*q0/(qdot0+gamma*q0));          %Winkel theta
q1    = B.*exp(-gamma*t).*sin(omega*t+theta);     %Schwach gedämpfter 
                                                  %Oszillator
q1e   = B.*exp(-gamma*t);                         %Einhüllende

%% 
% Fall mit starker Dämpfung
eta12 = 1.2;                             %Reibungsterm 1
fac   = eta12/2/m;
gamma1= fac - sqrt(fac^2-k/m);           %Kritischer Parameter gamma1
gamma2= fac + sqrt(fac^2-k/m);           %Kritischer Parameter gamma2
Disc  = fac^2-k/m;                       %Diskriminator: hier positiv
if Disc <= 0                         
    disp('keine starke Dämpfung'); 
    return; 
end
A     = q0 -(qdot0+gamma1*q0)/(gamma1-gamma2);    %Konstante A
B     = (qdot0+gamma1*q0)/(gamma1-gamma2);        %Konstante B
q2    = A.*exp(-gamma1*t) + B.*exp(-gamma2*t);    %Stark gedämpfter 
                                                  %Oszillator


%% 
% Fall mit kritischer Dämpfung
eta13 = 2*m*omega0;                      %Reibungsterm 1
gamma = eta13/2/m;                       %Kritischer Parameter gamma
A     = q0;                              %Konstante A
B     = (qdot0+gamma1*q0);               %Konstante B
q3    = A.*exp(-gamma1*t) + B.*t.*exp(-gamma2*t); %Kritisch gedämpfter 
                                                  %Oszillator

                                                  
%%
% Graphik

figure();
title('Harmonischer Oszillator');
plot(t,q,'Color',Colors(4,:),'Linewidth',1,'LineStyle',Style(1));
hold on
plot(t,q1, 'Color',Colors(2,:),'Linewidth',1,'LineStyle',Style(2));
plot(t,q1e,'Color',Colors(2,:),'Linewidth',1,'LineStyle',Style(1));
plot(t,q2,'Color',Colors(3,:),'Linewidth',1,'LineStyle',Style(4));
plot(t,q3,'Color',Colors(4,:),'Linewidth',1,'LineStyle',Style(5));
ylabel('Auslenkung  \it{q} \rm in m','FontSize',14);
xlabel('Zeit \it{t} \rm in s ','FontSize',14);
legend('ungedämpft','schwache Dämpfung','Einhüllende',...
       'überkritische Dämpfung', 'kritische Dämpfung','location',...
       'south','NumColumns',2);
legend box off
axis([0 8 -4 4 ])
grid on
set(gca,'Fontsize', 16);


%%
% Teil b)

%%
% Variable Anregungsfrequenz

omegaD = omega0*[0.5; 0.9; 1.1; 1.5];    %Anregungsfrequenz
eta    = [0.01; 0.05; 0.25; 0.75];       %Reibunsgkeoffizient
FA     = 0.5;                            %Amplitude der Anregung
phiA   = 0;                              %Phase der Anregung

tmax  = 50;                              %Maximalzeit
NPkt  = 500;                             %Anzahl Punkte
t     = linspace(0,tmax,NPkt);           %Zeitbereich

for index = 1:4
    indexo  = index;
    indexe  = 2;
    omegaA = omegaD(indexo);
    gamma  = eta(indexe)/2/m;                %Kritischer Parameter gamma
    Disc   = (2*gamma*omegaA).^2+(omega0^2-omegaA^2).^2; 
    A      = FA/m/sqrt(Disc);                %Amplitude durch ANregung
    NN     = (omega0^2-omegaA^2);            %Nenner
    if NN==0
        NN = 1.e-4; 
    end
    if(omegaA <= omega0)
      phi = atan(2*gamma*omegaA/NN);      %Phase zwischen Anregung und Osz. 
    else
      phi = pi + atan(2*gamma*omegaA/NN); %Shift um pi bei omegaA > omega0
    end

    delta = theta - phi;              %Phase der Lösung
    qp = A*cos(omegaA*t + delta);     %Spezielle Lösung der inhomogenen DGL    
    Disc=omega0^2-gamma^2;            %Positiv für schwache Dämfung
    if Disc <= 0                       
        disp('zu starke Dämpfung'); 
        return; 
    end
    omega1 = sqrt(Disc);                 
    theta  = atan(omega1*q0/(qdot0+gamma*q0));      %Winkel theta
    B      = sqrt(q0^2+(qdot0+gamma*q0)^2/omega1^2);%Konstante B 
    qh     = B*exp(-gamma*t).*sin(omega1*t+theta);  %Lösung homogene DGL
    q(index,:) = qp + qh;
    ttlstr(index,:)=string(strcat('\gamma =  ',num2str(gamma,3),...
                   '  \omega_A/\omega_0 = ',num2str(omegaA/omega0,3)));
end

figure();
% text(2,A*1.75,str,'FontSize',12,'Color','red');
for index = 1:4
    hold on;
    subplot(2,2,index)
    plot(t,q(index,:),'Color',Colors(index+1,:),'Linewidth',1,'LineStyle',...
         Style(1));
    title(ttlstr(index,:),'FontSize',12)
    ylabel('q(t) in m');
    xlabel('t in s','FontSize',8);
    axis([0 tmax 1.25*min(q(2,:)) 1.25*max(q(2,:))])
    grid on
    set(gca,'Fontsize', 16);
end

%%
% Variable Reibung

for index = 1:4
    indexo  = 2;
    indexe  = index;
    omegaA = omegaD(indexo);
    gamma  = eta(indexe)/2/m;                %Kritischer Parameter gamma
    Disc   = (2*gamma*omegaA).^2+(omega0^2-omegaA^2).^2; 
    A      = FA/m/sqrt(Disc);                %Amplitude durch ANregung
    NN     = (omega0^2-omegaA^2);            %Nenner
    if NN==0
        NN = 1.e-4; 
    end
    if(omegaA <= omega0)
      phi = atan(2*gamma*omegaA/NN);      %Phase zwischen Anregung und Osz. 
    else
      phi = pi + atan(2*gamma*omegaA/NN); %Shift um pi bei omegaA > omega0
    end

    delta = theta - phi;              %Phase der Lösung
    qp = A*cos(omegaA*t + delta);     %Spezielle Lösung der inhomogenen DGL    
    Disc=omega0^2-gamma^2;            %Positiv für schwache Dämfung
    if Disc <= 0                       
        disp('zu starke Dämpfung'); 
        return; 
    end
    omega1 = sqrt(Disc);                 
    theta  = atan(omega1*q0/(qdot0+gamma*q0));      %Winkel theta
    B      = sqrt(q0^2+(qdot0+gamma*q0)^2/omega1^2);%Konstante B 
    qh     = B*exp(-gamma*t).*sin(omega1*t+theta);  %Lösung homogene DGL
    q(index,:) = qp + qh;
    ttlstr(index,:)=string(strcat('\gamma =  ',num2str(gamma,3),...
                   '  \omega_A/\omega_0 = ',num2str(omegaA/omega0,3)));
end

figure();
% text(2,A*1.75,str,'FontSize',12,'Color','red');
for index = 1:4
    hold on;
    subplot(2,2,index)
    plot(t,q(index,:),'Color',Colors(index+1,:),'Linewidth',1,'LineStyle',...
         Style(1));
    title(ttlstr(index,:),'FontSize',12)
    ylabel('q(t) in m');
    xlabel('t in s','FontSize',8);
    axis([0 tmax 1.25*min(q(2,:)) 1.25*max(q(2,:))])
    grid on
    set(gca,'Fontsize', 16);
end


%%
% Teil c)

%%
% Variable Anregungsfrequenz

omegaMin = 0.001*omega0;
omegaMax = 2.5*omega0;
NPkt     = 500;
omegaA = linspace(omegaMin,omegaMax,NPkt);  
etamin = 0.02;                          %Minimum eta
etamax = 2*m*omega0/sqrt(2);            %Maximum eta
etastep=(etamax-etamin)/5;  

figure()
subplot(1,3,1)
jPlot = 0;
for eta=etamin:etastep:etamax              
    gamma = eta/2/m;                    %gamma
    jPlot = jPlot+1;
    Disc  = (2*gamma*omegaA).^2+(omega0^2-omegaA.^2).^2; 
    A     = FA/m./sqrt(Disc);           %Amplitude
    semilogy(omegaA/omega0,A,'Color',Colors(jPlot,:),'Linewidth',1,... 
             'LineStyle',Style(1));                   %Plot Amplitude
    hold on
    omega_res = sqrt(omega0^2-2*gamma^2);             %Resonanzfrequenz
    Amax      = FA/2/m/gamma/sqrt(omega0^2-gamma^2);  %Max @ Resonanz
    line([omega_res/omega0;omega_res/omega0],[Amax;Amax] ,...
       'Color',Colors(jPlot,:),'Marker','+','MarkerSize',8,'LineWidth',1);
    str=cat(2,'\gamma=',num2str(gamma,2));
    text(omega_res/omega0+0.01,1.1*Amax,str,'FontSize',10,...
        'Color',Colors(jPlot,:));
end
title('Amplitude über Anregungsfrequenz','FontSize',14)
ylabel('Amplitude \it A','FontSize',14);
xlabel('\omega_A/\omega_0','FontSize',14);
axis([0 2.5 0.1 20])
grid on
set(gca,'Fontsize', 16);

subplot(1,3,2)
jPlot = 0;
for eta=etamin:etastep:etamax              
   gamma = eta/2/m;                    %gamma
   jPlot = jPlot+1;
   NN  = (omega0^2-omegaA.^2);  
   for i=1:NPkt
     if(omegaA(i)<=omega0)
       phi(i)= atand(2*gamma*omegaA(i)/NN(i));       %Phasendifferenz    
     else
       phi(i)= 180+atand(2*gamma*omegaA(i)/NN(i));  %Phasenshift um pi 
     end
   end
   plot(omegaA/omega0,phi,'Color',Colors(jPlot,:),'Linewidth',1,... 
             'LineStyle',Style(1));                   %Plot Phase
   hold on
   str=cat(2,'\gamma=',num2str(gamma,2));
   text(1.05,(1+0.1)*phi(100)+5,str,'FontSize',10,...
        'Color',Colors(jPlot,:));
end
title('Phasendifferenz','FontSize',14)
ylabel('Phase \phi','FontSize',14);
xlabel('\omega_A/\omega_0','FontSize',14);
axis([0 2.5 0 180])
grid on
set(gca,'Fontsize', 16);

subplot(1,3,3)
jPlot = 0;
for eta=etamin:etastep:etamax              
   gamma = eta/2/m;                        %gamma
   jPlot = jPlot+1;
   Disc  = (2*gamma*omegaA).^2+(omega0^2-omegaA.^2).^2; 
   NN  = (omega0^2-omegaA.^2);  
   A     = FA/m./sqrt(Disc);               %Amplitude
   Power = 0.5*FA.*A.*omegaA.*sind(phi);   %Power
   semilogy(omegaA/omega0,Power,'Color',Colors(jPlot,:),'Linewidth',1,... 
             'LineStyle',Style(1));        %Plot Power
   hold on
   [Power,j]=max(Power);            %Power ist maximal
   str=cat(2,'\gamma=',num2str(gamma,2));
   text(1.05,Power + 0.02,str,'FontSize',10,...
        'Color',Colors(jPlot,:));  
end
title('Leistungstransfer ','FontSize',14)
ylabel('Leistung \it P ','FontSize',14);
xlabel('\omega_A/\omega_0','FontSize',14);
axis([0 2.5 0.01 10])
grid on
set(gca,'Fontsize', 16);
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------

