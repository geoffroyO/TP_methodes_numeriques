
// Retourne la fréquence d'échantillonage de la transformée de Fourier discrète
function [freqs]=fftfreq(N, L)
    k = modulo(N, 2)
    if k == 0 then
        freqs = 2*%i*%pi * cat(2, linspace(0, N/2 -1, N/2 ), linspace(-N/2, -1, N/2))/L; // Cas N pair
    else  
        freqs = 2*%i*%pi * cat(2, linspace(0, (N-1)/2, (N+1)/2), linspace(-(N-1)/2, -1, (N-1)/2))/L; // Cas N impair
    end
endfunction


// Résolution de l'équation de Poisson en dimension 2 en utilisant la FFT
//    laplacien(psi) = f
// Entrée: f de taille (Ny,Nx) sur une domaine de taille (Ly,Lx)
// Sortie: psi, solution de l'équation
function [psi]=poisson_2d(f, Nx, Ny, Lx, Ly)
    kx = fftfreq(Nx, Lx)
    ky = fftfreq(Ny, Ly)
    psi_hat = zeros(Ny,Nx)
    f_hat = fft(f, "nonsymmetric")
    for p=1:Ny
        for q=1:Nx
            if p==1 & q==1 then 
                    psi_hat(1,1) = 0; // Cas ou kx et ky sont simultanément nuls
            else 
                psi_hat(p,q) = f_hat(p,q)/(kx(q)**2 + ky(p)**2)
            end
        end
    end
    psi = real(ifft(psi_hat, "nonsymmetric")) // Transformée de Fourier inverse réelle
endfunction
// Résolution de l'équation de Poisson avec rot en dimension 2 en utilisant la FFT
//    laplacien(Ux) = -dW/dy
//    laplacien(Uy) = +dW/dx
// Entrée: champs de vorticité W de taille (Ny,Nx) sur un domaine de taille (Ly,Lx)
// Sortie: Ux et Uy, vitesses solution des équations
function [Ux,Uy]=poisson_curl_2d(W, Nx, Ny, Lx, Ly)
    kx = fftfreq(Nx, Lx)
    ky = fftfreq(Ny, Ly)
    u_x_hat = zeros(Ny,Nx)
    u_y_hat = zeros(Ny,Nx)
    w_hat = fft(W, "nonsymmetric")
    for p=1:Ny
        for q=1:Nx
            if p==1 & q==1 then 
                    u_x_hat(1,1) = 0; // Cas où kx et ky sont simultanément nuls
                    u_y_hat(1,1) = 0;
            else 
                u_x_hat(p,q) = (-ky(p)/(kx(q)**2 + ky(p)**2))*w_hat(p,q); // Calcul des coefficients de Fourier
                u_y_hat(p,q) = (kx(q)/(kx(q)**2 + ky(p)**2))*w_hat(p,q);
            end
        end
    end
    Ux = real(ifft(u_x_hat, "nonsymmetric"))
    Uy = real(ifft(u_y_hat, "nonsymmetric"))
endfunction

