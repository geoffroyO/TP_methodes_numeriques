
function [A]=cholesky_fact(A)
//-----------------------
//TODO factorise en place
//----------------------
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
