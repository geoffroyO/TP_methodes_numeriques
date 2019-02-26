
// Retourne la fréquence d'échantillonage de la transformée de Fourier discrète
function [freqs]=fftfreq(N, L)
    i = modulo(N, 2)
    if i == 0 then
        freqs = (2*%i*%pi/L) * cat(1, linspace(0, N/2 -1, N/2 ), linspace(-N/2, -1, N/2));
    else  
        freqs = (2*%i*%pi/L) * cat(1, linspace(0, (N-1)/2, (N+1)/2), linspace(-(N-1)/2, -1, (N-1)/2));
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
            psi_hat(p,q) = f_hat(p,q)/(kx(q)**2 + ky(p)**2)
        end
    end
    psi = real(ifft(psi_hat, "nonsymmetric"))
endfunction

// Résolution de l'équation de Poisson avec rot en dimension 2 en utilisant la FFT
//    laplacien(Ux) = -dW/dy
//    laplacien(Uy) = +dW/dx
// Entrée: champs de vorticité W de taille (Ny,Nx) sur un domaine de taille (Ly,Lx)
// Sortie: Ux et Uy, vitesses solution des équations
function [Ux,Uy]=poisson_curl_2d(W, Nx, Ny, Lx, Ly)
    // TODO: Calculer Ux et Uy à partir de la vorticité par FFT avec l'option 'nonsymmetric'
endfunction

