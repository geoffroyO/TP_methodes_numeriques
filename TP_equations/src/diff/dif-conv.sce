Nx=200;
Nt=200;
kappa=0.0001;
exec("my_cholesky.sce")

function y=phi_0(x)
//--------------------
//TODO condition initiale
//--------------------
endfunction

function y=conv(x)
//--------------------
//TODO fonction de convection
//--------------------
endfunction


//--------------------
//TODO initialiser phi et assembler les matrices M et N
//-------------------
fin=Nt;
for i=1:fin
    phi=my_choleski(N,M*phi);
end
scf;
plot(maillage,[phi_i phi]);
