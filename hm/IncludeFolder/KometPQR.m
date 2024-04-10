% -------------------------------------------------------------------------
% KometPQR.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Berechnet die heliozentrisch-ekliptikalen Koordinaten der Kometen über
% Gaussche Vektoren
%
% Eingabe:
%   T       Zeit als Julianisches Datum 
%   GM      Produkt aus Gravitationskonstante und Sonnenmasse AE^3*d°-2
%   BaPaK   Bahnparameter der Kometen
%   k       Nr des Kometen
% Ausgabe:
%  v        Geschwindigkeit des Kometen in kartesischen Koordinaten
%  xyz      xyz Koordinaten des Kometen in kartesischen Koordinaten
%  Time     Zeitstuetzstellen in JD
%  Name     Name des Kometen
% -------------------------------------------------------------------------

function KometData=KometPQR(GM, T, BaPaK,k)
    pi2=2*pi;
    %   T0n      Periheldurchgang in Julianischem Datum 
    %   qPn      Periheldistanz in AE
    %   ePn      Exzentrizität
    ePn=BaPaK.eP(k);
    qPn=BaPaK.qP(k);
    T0n=BaPaK.T0(k);
    Tn=(T-2451545.0);
    T0n=(T0n-2451545.0);
    delta=1-ePn;
    invax=delta/qPn;
    tau=sqrt(GM)*(Tn-T0n);
    % Mittlere Anomalie 
    M  = tau*sqrt(invax*invax*invax);
    % Bahnkoordinaten aus  exzent. Anomalie für Ellipsenbahn
    aPn=1/invax;
    kf=sqrt(GM/aPn);
    Ecc=EAnom(M,ePn);
    cosE=cos(Ecc);
    sinE=sin(Ecc);
    fac=sqrt((1-ePn)*(1+ePn));
    rho=1-ePn.*cosE;
    x=aPn*(cosE-ePn); 
    y=aPn*fac.*sinE; 
    nu=2*atan(sqrt((1+ePn)/(1-ePn)).*tan(Ecc/2));
    R=aPn.*rho;
    vx=-kf.*sinE./rho;
    vy=kf*cosE*fac./rho;
    Temp.x_ell   = [x;y;R*0]; 
    Temp.v_ell   = [vx; vy; 0.0*vx];

    % Bahnkoordinaten aus Stumpff-Funktionen für Parabelbahn
    fac =0.5*ePn;
    E2=100;
    E20=0;
    iter=0;
    kf=sqrt(GM/(qPn*(1+ePn)));
    while abs(E2-E20)>eps
      iter=iter+1;
      E20 = E2;
      A = 1.5*sqrt(fac/(qPn*qPn*qPn)).*tau;  
      B = (sqrt(A.*A+1.0)+A).^(1/3);
      u  = B - 1.0./B;  
      u2 = u.*u;  
      E2 = u2*(1.0-ePn)./fac;
      [c1,c2,c3]=Stumpff(E2); 
      fac = 3.0*ePn.*c3;
    end
    R  = qPn*(1 + u2.*c2*ePn./fac );
    x=qPn*(1.0-u2.*c2./fac);
    y=qPn*sqrt((1.0+ePn)./fac).*u.*c1;
    z=0*x;
    vx=-kf*y./R;
    vy=kf*(x./R+ePn); 
    Temp.x_par   = [x;y;z]; 
    Temp.v_par   = [vx; vy; 0.0*vx];

    % Optimieren der Bahn (entweder Ellipse oder Parabel) 
    for n=1:length(M)
      if (abs(M(n)) > 0.2) 
           Komet(k).xyz(:,n)= Temp.x_ell(:,n);
           Komet(k).v(:,n)  = Temp.v_ell(:,n);
      else
           Komet(k).xyz(:,n)= Temp.x_par(:,n);
           Komet(k).v(:,n)  = Temp.v_par(:,n);
      end
    end
    
    % Winkeldistanz vom Fruehlingspunkt bis zum aufsteigenden Knoten
    % in der Ekliptik in Bogenmaß
    OmegaPn=deg2rad(BaPaK.OmegaP(k));
    % Perihelargument in Bogenmaß
    omegaPn=deg2rad(BaPaK.omegaP(k));
    % Inklination 
    iPn=deg2rad(BaPaK.iP(k));
    
    % PQR Matrixmultiplikation
    yin=[Komet(k).xyz(1,:);Komet(k).xyz(2,:);Komet(k).xyz(3,:)];
    yout=mtimes(PQR(-OmegaPn,-iPn,-omegaPn),yin);
    Komet(k).xyz=yout;
    
    % heliozentrische Geschwindigkeitsvektoren des Kometen
    yin=[Komet(k).v(1);Komet(k).v(2);Komet(k).v(3)];
    yout=mtimes(PQR(-OmegaPn,-iPn,-omegaPn),yin);
    Komet(k).v=yout;
    
    Komet(k).Name = BaPaK.Name(k);
    Komet(k).Time = T;
    KometData=Komet(k);
end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------