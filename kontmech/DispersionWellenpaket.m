% -------------------------------------------------------------------------
% DispersionWellenpaket.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik des Kontinuums" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm zeigt Dispersion eines Gauß-förmigen Wellenpakets.
% 
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

%% Parameter

%Parameterspace für Raumkoordinate x
xend   = 4000;
xsteps = 40000;
x      = linspace(-xend/10,xend*9/10,xsteps);
deltax = xend/xsteps;
x0     = 100;
b      = 15;                          % Anfangsbreite in m

%Parameterspace für Zeitkoordinate t
tend   = 200;
tsteps = 40000;
t = linspace(-tend,tend,tsteps);
deltat = tend/tsteps;
t0     = 0;                             % Anfangspulsbreite

vph0    = 10.0;                         % mittlere Phasengeschwindigkeit m/s
vd      = -5;                           % Dispersionsparameter m/s
lambda0 = 4.0;                          % Wellenlänge 1 in m
k0      = 2*pi/lambda0 ;                % Wellenzahl  in m/s
omega0  = k0*vph0;

%% Berechnung Ausgangs-Wellenpaket und Spektrum
psi_inp = exp(-((x-x0)/(b)).^2) .* exp(1i*k0*x);
% Berechnung Fouriertransformation und Wellenzahlen
PSI_inp=fft(psi_inp);
PSI_shift=fftshift(psi_inp);
kx = 1/deltax/xsteps;

% Berechnung Vektoren der Raumfrequenzen und Geschwindigkeiten 
for kz=1:xsteps
  kw(kz)    = 2*pi*kx*(kz-1);
  vph(kz)   = vph0 + vd*(abs(kw(kz))-k0)/k0;
  omega(kz) = -abs(vph(kz)*kw(kz)); % Damit nach rechts laufende Welle muss
                                    % omega > 0 sein. 
  if  kw(kz) > 0 
           omega(kz) = -omega(kz);
  end
end
%% Kontrollbild (non shifted)

figure('name','Kontrollbild')
h2 = title("v_{ph}(k) und \omega(k)");
xlabel("k/(2\pi) (1/m)")

yyaxis left
plot(kw/2/pi,vph)
ylabel("v_{ph}(k) in m/s")
% axis([0 5 -5 15]);

yyaxis right
plot(kw/2/pi,omega)
ylabel("\omega(k) in 1/s")

grid on
set(h2, 'FontSize',14,'FontName','Times','FontWeight','normal');
set(gca,'FontSize',14,'FontName','Times');


%% Wellenbild zu zwei verschiedenen Zeiten
figure('name','Wellenpaket zu zwei verschiedenen Zeiten')
tmess=150;
tmess_str = strcat('t = ',num2str(tmess,5));
tmess_str = strcat(tmess_str,' s');
hold on
plot(x,real(ifft(PSI_inp.*exp(-1i*omega*150))));
plot(x,real(psi_inp));
axis([0 1200 -1 1]);
grid on
xlabel("x (m)")
ylabel("\psi")
h2 = title(strcat('Wellenpaket in x-Domain @t=0 s und @ ',tmess_str));
set(h2, 'FontSize',14,'FontName','Times', 'FontWeight','normal');
set(gca,'FontSize',14,'FontName','Times');

figure('name','Wellenpaket Spektrum')
hold on
subplot(2,1,1)
plot(kx*(-xsteps/2:xsteps/2-1),abs(PSI_shift))
title("fft Spectrum")
xlabel("k/(2\pi) (1/m)")
ylabel("Betrag ")
axis([0 2 -inf inf]);
grid on
h2 = title(strcat('Powerspectrum @',tmess_str)); 
set(h2, 'FontSize',14,'FontName','Times','FontWeight','normal');
set(gca,'FontSize',14,'FontName','Times');

subplot(2,1,2)
hold on
plot(kx*(-xsteps/2:xsteps/2-1),omega*tmess/2/pi)
title("Spectrum des Filters")
xlabel("k/(2\pi) (1/m)")
ylabel("Phase \omega(k)t in 2\pi  ")
axis([0 2 -inf inf]);
axis([0 2 0 2e5]);
grid on
h2 = title(strcat('Phase des Filters @',tmess_str)); 
set(h2, 'FontSize',14,'FontName','Times','FontWeight','normal');
set(gca,'FontSize',14,'FontName','Times');

%% Simulation der Ausbreitung der Wellenfunktion bis zum Zeitpunkt tmess
% tmess =2000;
% 
figure('name','Simulation Ausbreitung Wellenpaket');
hold on
grid on
axis([0 1200 -1 1]);
xlabel("x (m)")
ylabel("Amplitude")
plot(x,real(psi_inp));
izend = 31;
tstr(izend,8) ="        ";
for iz = 1:izend
   t(iz) = (iz-1)*5; % Zeit in s
   tempstr = strcat(num2str(t(iz), 5),' s');
   tstr(iz,:) = strcat('t = ',tempstr);
   psi_out(iz,:) = real(ifft(PSI_inp(:).*exp(-1i*omega(:)*t(iz))));
end
hp2 = text(1000,0.95,tstr(1));
pause(0.1)
for iz = 2:izend
    hp2.Visible = 'off';
    hp2 = text(1000,0.95,tstr(iz));
    hold on
    hp1 = plot(x(:),psi_out(iz,:),'color',Colors(2,:));
    pause(0.2)
    hp1.Visible = 'off';
end
hp = plot(x(:),psi_out(izend,:),'color',Colors(2,:));

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------

