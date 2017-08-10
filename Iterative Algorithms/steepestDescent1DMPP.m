%% John Bowllan Steepest Descent 1DMPP function file

function [xApprox,time] = steepestDescent1DMPP(Ax,fVec,uOld,tol,itMax2,n)
% This function, called by the helper file, carries out the steepest 
% descent algorithm. We must input a matrix A, vector b, an initial
% guess x0 for the solution, a certain tolerance tol, and a potential 
% cap on the number of allowed iterations itMax. This function will 
% output the number of iterations it took the algorithm to find a solution
% to the system within a certain tolerance and the time it took the algorithm
% to run. 
 

xOld = uOld;
rOld = fVec - Ax;

% iteration, xCoord, yCoord are solely for constructing a table at the end


tic % start timer
for k = 1:itMax2
    % Calculate alpha, xNew, rnew, and necessary norms and suppress;
    ArOld = [0; -rOld(1:n-1)] + 2*rOld  - [rOld(2:n); 0];
    alpha = dot(rOld,rOld)/dot(rOld,ArOld);
    xNew = xOld + alpha*rOld;
    rNew = fVec - ([0; -xNew(1:n-1)] + 2*xNew  - [xNew(2:n); 0]);
    resNorm = norm(rNew); % ADD MORE CONVERGENCE CRITERIA
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