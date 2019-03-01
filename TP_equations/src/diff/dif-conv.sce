exec("my_cholesky.sce")
Nx=200;
dx=1/Nx;
Nt=200;
dt=1/Nt;
kappa=0.0001;

function y=phi_0(x)
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

function y=conv(x)
    y = 0.4*(x-0.25)
endfunction

phi=zeros(Nx,1);
for n=1:(Nx)
    phi(n) = phi_0((n-1)*dx)
end

N=zeros(Nx, Nx)
M=zeros(Nx, Nx)
N(1,1) = 1+2*kappa*(dt/dx**2)
N(1,2) = -kappa*dt/dx**2
N(1,Nx) = -kappa*dt/dx**2
N(Nx,1) = -kappa*dt/dx**2
N(Nx,Nx-1) = -kappa*dt/dx**2
N(Nx,Nx) = 1+2*kappa*(dt/dx**2)
for k=2:(Nx-1)
    N(k,k) = 1+2*kappa*(dt/dx**2);
    N(k,k-1) = -kappa*dt/dx**2;
    N(k,k+1) = -kappa*dt/dx**2;
end
M(1,1) = 1-2*conv(0)**2;
M(1,2) = -conv(0)*dt/(2*dx) + conv(0)**2*dt**2/(2*dx**2)
M(1,Nx) = conv(0)*(dt/2*dx) + conv(0)**2*(dt**2/2*dx**2)
M(Nx,Nx) = 1-2*conv((Nx-1)*dx)**2;
M(Nx,1) = -conv((Nx-1)*dx)*dt/(2*dx) + conv((Nx-1)*dx)**2*dt**2/(2*dx**2)
M(Nx,Nx-1) = conv((Nx-1)*dx)*(dt/2*dx) + conv((Nx-1)*dx)**2*(dt**2/2*dx**2)
for k=2:(Nx-1)
    M(k,k) = 1-2*conv(((k-1)*dx))**2;
    M(k,k+1) = -conv((k-1)*dx)*dt/(2*dx) + conv((k-1)*dx)**2*dt**2/(2*dx**2);
    M(k,k-1) = conv((k-1)*dx)*(dt/2*dx) + conv((k-1)*dx)**2*(dt**2/2*dx**2);
end
fin=Nt;
exec("my_cholesky.sce")

for i=1:fin
    phi = my_cholesky(N,M*phi);
end
