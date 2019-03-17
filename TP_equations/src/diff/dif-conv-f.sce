
function dt=calcul_dt(u,dx)
  dt=dx/max(abs(u));
endfunction

function vort_s=solveur_1D(vort_p,u,m,kappa,dt,dx)
  M = zeros(Nx,Nx) // Initialisation des matrices M et N
  N = zeros(Nx,Nx)
  for i=1:Nx
    N(i,i) = 1 + 2*kappa*dt/(dx^2)
    M(i,i) = 1 - (u(i)*dt/dx)^2
    if i<Nx then
        N(i,i+1) = -kappa*dt/dx^2
        M(i,i+1) =-u(i)*(dt/(2*dx)) + ((u(i)*dt/dx)^2)/2
    end
    if  i>1 then
        N(i,i-1) = -kappa*dt/dx^2
        M(i,i-1) = u(i)*(dt/(2*dx)) + ((u(i)*dt/dx)^2)/2
    end    
end
  N(1,Nx) = -kappa*dt/dx^2
  M(1,Nx) = u(1)*(dt/(2*dx)) + ((u(1)*dt/dx)^2)/2
  N(Nx,1) = -kappa*dt/dx^2
  M(Nx,1) = -u(Nx)*(dt/(2*dx)) + ((u(Nx)*dt/dx)^2)/2
  vort_s = umfpack(sparse(N),"\",M*vort_p') // Résolution du système linéaire
endfunction


function vort_s=solveur_2D(vort_p, ux, uy, Nx, Ny, kappa, dt, dx, dy) // Solveur 2D avec méthode de splitting
  vort_s = vort_p;
  for k=1:Nx
    vort_s(k,:)=solveur_1D(phi(k,:), ux(k,:), Nx, kappa, dt, dx); 
  end
  for k=1:Ny
    vort_s(:,k) =solveur_1D(phi(:,k)', uy(:,k), Ny, kappa, dt, dy)
  end
endfunction
  

