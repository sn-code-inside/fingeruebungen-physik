function PlotCircle (xM,yM,rC,Col,LW)
    u=linspace(0,360,360);
    nx=zeros(360);
    ny=zeros(360);
    Nmax = 100; 
    for k=1:Nmax
        nx= rC*k*cosd(u)/Nmax+xM;
        ny= rC*k*sind(u)/Nmax+yM;
        plot(nx,ny,'LineWidth',LW,'Color',Col);
    end
end

