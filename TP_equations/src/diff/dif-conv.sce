Nx=200;
dx=1/Nx;
Nt=200;
dt=1/Nt;
kappa=0.01;

function y=phi_0(x) // Condition initiale
    if (0 <= x) & (x < 0.25) then 
        y = 0
    elseif (0.25 <= x) & (x < 0.375) then
        y = 2*(x-0.25)
    elseif (0.375 <= x) & (x < 0.5) then
        y = 2*(0.5 - x)
    elseif (0.5 <= x) & (x <= 1) then
        y =0
    end
endfunction

function y=conv(x) // Vitesse de convection 
    y = 0.4*(x-0.25)
endfunction

phi=zeros(Nx,1); // Initialisation de de phi
for n=1:(Nx)
    phi(n) = phi_0((n-1)*dx)
end
phi_i = phi
// Initialisation de M et N
N=zeros(Nx, Nx)
M=zeros(Nx, Nx)
N(1,1) = 1+2*kappa*(dt/dx^2)
N(1,2) = -kappa*dt/dx^2
N(1,Nx) = -kappa*dt/dx^2
N(Nx,1) = -kappa*dt/dx^2
N(Nx,Nx-1) = -kappa*dt/dx^2
N(Nx,Nx) = 1+2*kappa*(dt/dx^2)
for k=2:(Nx-1)
    N(k,k) = 1+2*kappa*(dt/dx^2);
    N(k,k-1) = -kappa*dt/dx^2;
    N(k,k+1) = -kappa*dt/dx^2;
end
M(1,1)=1-2*conv(0)^2*dt^2/(2*dx^2);
M(1,2)=conv(0)*dt/(2*dx)*(conv(0)*dt/dx - 1);
M(1,Nx)=conv(0)*dt/(2*dx)*(conv(0)*dt/dx + 1);
M(Nx,1)=conv((Nx-1)*dx)*dt/(2*dx)*(conv((Nx-1)*dx)*dt/dx - 1);
M(Nx,Nx-1)=conv((Nx-1)*dx)*dt/(2*dx)*(conv((Nx-1)*dx)*dt/dx + 1);
M(Nx,Nx)=1-2*conv((Nx-1)*dx)^2*dt^2/(2*dx^2);
for i=2:(Nx-1)
    M(i,i)=1-2*conv((i-1)*dx)^2*dt^2/(2*dx^2);
    M(i,i+1)=conv((i-1)*dx)*dt/(2*dx)*(conv((i-1)*dx)*dt/dx - 1);
    M(i,i-1)=conv((i-1)*dx)*dt/(2*dx)*(conv((i-1)*dx)*dt/dx + 1);
end

fin=Nt;

// Fonctions pour la méthode de Cholesky

function [A]=cholesky_fact(A)
    // Factorisation choleski de la matrice A
      [m,n]=size(A);
      if (m~=n) then
        print(%io(2), "error, not a square matrix");
      else
          A(1,1) = sqrt(A(1,1))
          for p = 1 : n
              
              for i = 1 : n
                  if i ~= p then
                    somme = 0
                      for k = 1 : p-1
                          somme = somme + A(i,k)*A(p,k)
                      end
                      A(i,p) = (A(i, p) - somme)/A(p,p)
                   end
              end
              if p+1<=n then
                somme = 0
                for k = 1 : p
                    somme = somme + A(p+1, k)**2
                end
                A(p+1,p+1) = (A(p+1, p+1) - somme)**(1/2)
              end
          end
      end
    endfunction
    
    function [y]=up_sweep_cholesky(A,x)
      [m,n]=size(A);
      if (m~=n) then
        print(%io(2), "error, not a square matrix");
      else
        y = zeros(n, 1)
        for i = n : -1 : 1
            S = 0
            for k  = n : -1 : (i+1)
                S = S + A(i,k)*y(k,1)
            end
            y(i,1) = (x(i,1)-S)/A(i,i)
            
        end
      end
    endfunction
    
    function [y]=down_sweep_cholesky(A,x)
      [m,n]=size(A);
      if (m~=n) then
        print(%io(2), "error, not a square matrix");
      else
        y = zeros(n, 1)
        for i = 1 : n
            S = 0
            for k  = 1 : (i-1)
                S = S + A(i,k)*y(k,1)
            end
            y(i,1) = (x(i,1)-S)/A(i,i)
            
        end
     end
    endfunction

// Calcul des phi par itération
T = cholesky_fact(N)
for i=1:Nt 
    y = down_sweep_cholesky(T, M*phi);
    phi = up_sweep_cholesky(T', y)
end

// On plot les résultats
X =0:1/Nx:(Nx-1)/Nx;
scf;
plot(X,[phi_i phi]);
xtitle('Resolution de lequation', 'x', '$\phi(x)$');
legend(['$\phi(t=0)$';'$\phi \quad  au \,temps \,final$']);

