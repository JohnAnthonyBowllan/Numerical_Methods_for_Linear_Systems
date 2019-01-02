function [xApprox,time,iteration,xCoord,yCoord] = steepestDescent(A,b,x0,tol,itMax)
% Steepest Descent Method for solving linear systems
% Inputs: 
% matrix A, vector b, initial solution x0, tolerance tol, and max iterations allowed itMax
% Outputs:
% approximate solution xApprox and time
 

xOld = x0;
rOld = b - A*xOld;
iteration = [];
xCoord = [];
yCoord = [];
tic % start timer
for k = 1:itMax
    % Calculate alpha, xNew, rnew, and necessary norms and suppress;
    alpha = dot (rOld,rOld)/dot(rOld,A*rOld);
    xNew = xOld + alpha*rOld;
    iteration = [iteration; k];
    xCoord = [xCoord; xNew(1)];
    yCoord = [yCoord; xNew(2)];
    rNew = b - A*xNew;
    resNorm = norm(rNew); 
    if (resNorm < tol)
        display(['Solution Converged in ' num2str(k) ' iterations'])
        xApprox = xNew;
        time = toc % stop timer and save time
        return % exit function file
    end
    rOld = rNew;
    xOld = xNew;
end
error(['Convergence in ' num2str(itMax) ' iterations failed'])
