
function [T]=cholesky_fact(A) // factorisation de cholesky, donne la matrice 
  [n, m] = size(A)            // triangulaire inferieur
  T = zeros(n, m)
  T(1,1) = sqrt(A(1,1))
  for j=2:n
    T(j,1) = A(1,j)/T(1,1)
  end;
  for i=2:(n-1)
    S = 0
    for k=1:(i-1)
      S = S + T(i,k)**2
    end;
    T(i,i) = sqrt( A(i,i) - S)
    for j=(i+1):n
      S1 = 0
      for k=1:(i-1)
        S1 = S1 + T(i,k)*T(j,k)
      end;
      T(j,i)=(A(i,j)-S1)/T(i,i)
    end;
  end;
  S2 = 0
  for i=1:(n-1)
    S2 = S2 + T(n,i)**2
  end;
  T(n,n) =sqrt(A(n,n) - S2)
endfunction

function [y]=up_sweep_cholesky(A,x)
  [m,n]=size(A);
  if (m~=n) then
    print(%io(2), "error, not a square matrix");
  else
//-------------------------
//TODO
//------------------------
  end
endfunction

function [y]=down_sweep_cholesky(A,x)
  [m,n]=size(A);
  if (m~=n) then
    print(%io(2), "error, not a square matrix");
  else
//-----------------------
//TODO
//----------------------
 end
endfunction

function [U]=my_cholesky(N,S)
//---------------
//TODO
//--------------
endfunction
