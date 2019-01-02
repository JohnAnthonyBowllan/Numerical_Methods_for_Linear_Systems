function [xApprox,time] = CGmethod(A,b,x0,tol,itMax)
% Conjugate Gradient Method for solving linear systems
% Inputs: 
% matrix A, vector b, initial solution x0, tolerance tol, and max iterations allowed itMax
% Outputs:
% approximate solution xApprox and time
xOld = x0;
rOld = b-A*xOld;
dOld = rOld;
bNorm = norm(b);
tic % start timer
for k = 1:itMax
    % calculate alpha
    alpha = dot(dOld,rOld)/dot(dOld,A*dOld); 
    % update x and r
    xNew = xOld  + alpha.*dOld;
    rNew = b-A*xNew;
    % calculate convergence criteria
    delta = xNew-xOld;
    deltaNorm = sqrt(delta'*delta);
    xNorm = sqrt(xNew'*xNew);
    resNorm = sqrt(rNew'*rNew);
    % test for convergence
    if (resNorm < tol * bNorm)
        display(['Solution Converged in ' num2str(k) ' iterations'])
        xApprox = xNew;
        time = toc % stop timer and save time
        return % exit function file
    end
    % update descent direction 
    beta = (rNew'*rNew)/(rOld'*rOld);
    dOld = rNew + beta*dOld; %dNew 
    % update r and x
    rOld = rNew;
    xOld = xNew;
end
error(['Convergence in ' num2str(itMax) ' iterations failed'])
