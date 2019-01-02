function [time,iteration] = SOR(A,n,b,x0,tol,itMax,omega)
% Successive Over-Relaxation Method for solving linear systems
% Inputs: 
% matrix A, dimension of matrix n, vector b, initial solution x0, tolerance tol, max iterations allowed itMax,
% and adjustable weight parameter omega
% Outputs:
% time and iterations to converge to solution within certain tolerance
it = 0;
xOld = x0;
tic
while it < itMax
    if (it == 500000 || it == 750000 || it == 1000000)
        display(it)
    end
    xNew(1) = (b(1) - A(1,2:n)*xOld(2:n))/A(1,1);
    for i = 2:n-1
       xNew(i) = (b(i) - A(i,1:i-1)*transpose(xNew(1:i-1)) ...
           -A(i,i+1:n)*xOld(i+1:n))/A(i,i);
    end
    xNew(n) = (b(n) - A(n,1:n-1)*transpose(xNew(1:n-1)))/A(n,n);
    xNewSOR = omega.*transpose(xNew) + (1-omega).*xOld;
    delta = xNewSOR - xOld;
    deltaNorm = norm(delta);
    xNorm = norm(xNewSOR);
    if (deltaNorm < tol*xNorm)
        display(['SOR solution converged in ' num2str(it) ' iterations']);
        iteration = it;
        xApprox = xNewSOR;
        %display(xApprox)
        time = toc;
        return
    end
    xOld = xNewSOR;
    it = it+1;
end
error(['Convergence in ' num2str(itMax) ' iterations failed'])
