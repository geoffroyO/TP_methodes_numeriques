
dx=Lx/Nx;
dy=Ly/Ny;

exec("dif-conv-f.sce")

maillage_x=linspace(0,(Nx-1)/Nx*Lx,Nx)';
maillage_y=linspace(0,(Ny-1)/Ny*Ly,Ny)';

cx=zeros(Ny,Nx); //composante x de la vitesse de convection
cy=zeros(Ny,Nx); //composante y de la vitesse de convection

phi=zeros(Ny,Nx);   //fonction à calculer
phi_i=zeros(Ny,Nx); //condtion initiale

for i=1:Nx
    for j=1:Ny
        phi(i,j)=phi_0((j-1)*dx,(i-1)*dy); // Initialisation de phi
        cy(i,j)=conv((j-1)*dx,(i-1)*dy)(1); // Initialisation de cx
        cx(i,j)=conv((j-1)*dx,(i-1)*dy)(2); // Initialisation de cy
    end
end

dt=min(calcul_dt(cx,dx),calcul_dt(cy,dy))
Nt=floor(Tf/dt)
phi_i=phi

for k=1:Nt
    phi = solveur_2D(phi, cx, cy, Nx, Ny, nu, dt, dx, dy); // On itère le solveur 2D
end
