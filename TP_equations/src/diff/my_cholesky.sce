
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

function [U]=my_cholesky(N,S)
    T = cholesky_fact(N)
    y = down_sweep_cholesky(T,S)
    U = up_sweep_cholesky(T.',y)
endfunction
