
dx=Lx/Nx;
dy=Ly/Ny;

exec("dif-conv-f.sce")

maillage_x=linspace(0,(Nx-1)/Nx*Lx,Nx)';
maillage_y=linspace(0,(Ny-1)/Ny*Ly,Ny)';

cx=zeros(Ny,Nx); //composante x de la vitesse de convection
cy=zeros(Ny,Nx); //composante y de la vitesse de convection

phi=zeros(Ny,Nx);   //fonction à calculer
phi_i=zeros(Ny,Nx); //condtion initiale

//------------------------------------------
//TODO remplir les tableaux cx cy phi phi_i
//------------------------------------------

dt=min(calcul_dt(cx,dx),calcul_dt(cy,dy));
Nt=floor(Tf/dt);
for k=1:Nt
    //-----------------
    // TODO
    //----------------
end
