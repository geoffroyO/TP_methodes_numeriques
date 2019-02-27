// Initialise W(x,y)
function [W]=init_field(y,x)
    W = 8*%pi**2*cos(2*%pi*x)*cos(2*%pi*y)
endfunction

// Solution de référence du problème 
function [Refx]=solution_fieldx(y,x)
    Refx = -2*%pi*cos(2*%pi*x)*sin(2*%pi*y)
endfunction

function [Refy]=solution_fieldy(y,x)
    Refy = 2*%pi*sin(2*%pi*x)*cos(2*%pi*y) 
endfunction

function plot_error(W, Refx, Ux, Refy, Uy, Errx, Erry, Y, X)  // Je sais pas mettre de légendes
    fig = gcf()
    subplot(341)
    plot3d(Y, X, W)
    subplot(342)
    plot3d(Y, X, Refx)
    subplot(343)
    plot3d(Y, X, Ux)
    subplot(344)
    plot3d(Y, X, Errx)
    subplot(345)
    plot3d(Y, X, Refy)
    subplot(346)
    plot3d(Y, X, Uy)
    subplot(347)
    plot3d(Y, X, Erry)
    xs2png(fig, "poisson_curl_error.png")
endfunction

// Fonction de test pour le solveur de Poisson
function test_poisson_curl(Lx, Ly, Nx, Ny)
    printf("::Testing poisson operator::")
    printf("\n  Domain size:    [%0.2f, %0.2f]", Lx, Ly)
    printf("\n  Discretization: [%i, %i]", Nx, Ny)
    
    // X[i] = i*dx avec dx = Lx/Nx et i=0..Nx-1
    // Y[i] = j*dy avec dy = Ly/Ny et j=0..Ny-1
    X = linspace(0.0, Lx*(Nx-1)/Nx, Nx)
    Y = linspace(0.0, Ly*(Ny-1)/Ny, Ny)
    
    printf("\n\n  Initializing field W(x,y).")
    W   = feval(Y, X, init_field)
    
    printf("\n  Initializing reference solution Refx, Refy.")
    Refx = feval(Y, X, solution_fieldx)
    Refy = feval(Y, X, solution_fieldy)
    dir  = get_absolute_file_path("test_poisson_curl.sce")
    file = dir+"poisson.sce" 
    printf("\n\n  Loading poisson_curl_2d function from file %s%s%s.", char(39), file, char(39))
    exec(file, -1)

    printf("\n\n  Computing Poisson solution Ux, Uy.")
   [Ux, Uy] = poisson_curl_2d(W, Nx, Ny, Lx, Ly)

    printf("\n  Computing error |Ux-Refx|(x,y).")
    Errx = abs(Ux-Refx)

    printf("\n  Computing error |Uy-Refy|(x,y).")
    Erry = abs(Uy-Refy)
    
    file = pwd()+"/poisson_curl_error.png"
    printf("\n\n  Plotting everything to %s%s%s.", char(39), file, char(39))
    plot_error(W, Refx, Ux, Refy, Uy, Errx, Erry, Y, X)
    
    printf("\n\n")
    mErrx = max(Errx)
    mErry = max(Erry)
    max_error = 1e-12
      
    if (mErrx > max_error) | (mErry > max_error) then
        printf(" TEST FAILURE ")
        exit(1)
    else
        printf(" TEST SUCCESS.\n")
        exit(0)
    end
endfunction


// Taille du domaine
Lx = 1.0
Ly = 1.0

// Discretisation du domaine
Nx = 64
Ny = 32

test_poisson_curl(Lx, Ly, Nx, Ny)