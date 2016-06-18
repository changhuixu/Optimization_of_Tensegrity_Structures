function [x,fval,xs,Ds] = ...
    newton_tr(x,Delta,Delta_max,eta,tol,trace)
% function [x,fval,xs,Ds] = ...
%	newton_tr(x,Delta,Delta_max,eta,tol,trace)
% Performs Newton's method from x.
% Returns x and func(x) as fval, and xs, the list of x values
% for each iteration.
%
%   Delta  -- original trust region size
%   Delta_max -- maximum trust region size
%   eta    -- quality control parameter for taking steps
%   tol    -- stopping tolerance: stop if ||grad f(x)||_2 < tol
%   trace  -- optional; if present and non-zero
%	then print out information about process

global mu;
mu = 100;
number_of_points = length(x);

% Check input data
if Delta > Delta_max || Delta <= 0 || eta < 0 || eta > 0.25
	fprintf('newton_tr: Error in data (eta, Delta, Delta_max)\n');
	fval = energy(x);
	return
end

% Initialization
xs = [];
Ds = [];

while 1	% forever do
%     if Delta < tol
%         fprintf('\n****** Delta value is too small. Bail out.\n\n');
%         return
%     end

    if trace ~= 0
        xs = [xs x];
        Ds = [Ds Delta];
    end
    % Set up model function
    [lambda,gLagrange] = findLambda(x);
    fval = energy(x);
    gval = gradientE(x);
    H = hessianLagrange(x,lambda);
   
    if trace ~= 0
        fprintf('newton_tr: function value = %g\n', fval);
    end
    %norm_GL= norm(gLagrange);
    %gLagrange = gLagrange;
    [m,n] = size(gval);
    if m < n
        gval = gval';
    end
    if trace ~= 0
        %fprintf('newton_tr: ||g|| = %g\n', norm(gval,2));
        fprintf('newton_tr: ||gLagrange|| = %g\n', norm(gLagrange));
        fprintf('newton_tr: ||C(x*)|| = %g\n', norm(constraintE(x)));
    end
    if norm(gLagrange,2) < tol && norm(constraintE(x)) < tol
        fprintf('*****************************************************\n');
        fprintf('newton_tr: ||gLagrange|| = %g\n', norm(gLagrange));
        fprintf('newton_tr: ||C(x*)|| = %g\n', norm(constraintE(x)));
        fprintf('newton_tr: minimum eigenvalue of N\''*HessL*N = %g\n',...
            min(eig(NHN)));
        return          % loop-ending condition
    end
  
    if trace ~= 0
        fprintf('newton_tr: Delta = %g\n', Delta);
    end

    % find phat
    [phat,N] = findPhat(x);
    % phat = phat;
    % norm_phat = norm(phat)
    NHN = N.'*H*N;
	if trace ~= 0
        fprintf('newton_tr: ||phat|| = %g\n', norm(phat));
        fprintf('newton_tr: minimum eigenvalue of N\''*HessL*N = %g\n',...
            min(eig(NHN)));
	end
    
    % compute the solution (p) for model function.
	if norm(phat) > Delta
        p = Delta*phat/norm(phat);
        at_bdry = 1;
        pvec = reshape(p,[3,number_of_points]);
        d = zeros(3*number_of_points,1);
        dvec = reshape(d,[3,number_of_points]);
    else
        [q,at_bdry] = findQ(phat,N,H,gval,Delta,0);
        p = phat + N*q;
        pvec = reshape(p,[3,number_of_points]);
        d = findD(x+pvec);
        dvec = reshape(d,[3,number_of_points]);

        if trace ~= 0
            if at_bdry
                fprintf('newton_tr: p at trust-region boundary\n');
            else
                fprintf('newton_tr: p interior to trust region\n');
            end
        end
	end

  % Compute "quality" parameter
    temp = norm(constraintE(x)+gradientC(x)*p)-norm(constraintE(x));
    denominator = gval'*p+0.5*p'*H*p+mu*temp;
    numerator = meritFunc(x+pvec+dvec)-meritFunc(x);
% 	if abs(denominator) < tol
%          denominator = -tol;
% 	end
    rho = numerator/denominator;

	if trace ~= 0
        fprintf('newton_tr: Quality parameter rho = %g\n', rho);
	end
	if rho < 0.25
        Delta = 0.25*norm(p+d,2);
    else
        if rho > 0.75 && at_bdry
            Delta = min(2*Delta,Delta_max);
        end
	end % if.rho < 0.25.else
    if rho > eta
        if trace ~= 0
            fprintf('newton_tr: take step\n');
        end
        x = x + pvec + dvec;
    end % if rho > eta
    %pause
end % while

end
