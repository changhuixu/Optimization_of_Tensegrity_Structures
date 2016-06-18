function [p,at_bdry] = dogleg_tr(g,H,Delta,trace)
% function [p,at_bdry] = dogleg_tr(g,H,Delta,trace)
%
% Uses the dogleg method to get an approximate solution
% for the trust region problem
%
%   min g'*p + (1/2)*p'H*p
% subject to
%   norm(p) <= Delta
%
% See Nocedal & Wright "Numerical Optimization", for example.

% Compute the Cauchy point
gHg = g'*H*g;
gnorm = norm(g);
if ( gHg > 0  && gnorm^3 < Delta*gHg )
  pC = -gnorm^2/gHg*g;
  at_bdry = 0;
else
  pC = -Delta*g/gnorm;
  at_bdry = 1;
end

% First check if H is positive semi-definite
[~,k] = chol(H);
if ( k > 0 ) % not positive definite
  % use Cauchy point
  p = pC;
  if trace > 0
    fprintf('dogleg_tr: H not positive definite; Using Cauchy point\n');
  end
else
  % Use dogleg with "Newton" method
  if ( at_bdry )
    % If the Cauchy point is already at the boundary
    % then return that.
    p = pC;
    if trace > 0
      fprintf('dogleg_tr: Using Cauchy point\n');
    end
  else
    % Else move from the Cauchy point toward the Newton point
    % until we hit the boundary of the trust region.
    pN = -H\g;
    if norm(pN) < Delta
      p = pN;
      at_bdry = 0;
      if trace > 0
	fprintf('dogleg_tr: Using full Newton step\n');
      end
    else
      % Solve ||pC+theta*(pN-pC)||^2 = Delta^2
      % This gives a quadratic equation a*theta^2+b*theta+c == 0;
      % pick positive root <= 1 using stable formula
      at_bdry = 1;
      pdiff = pN-pC;
      a = pdiff'*pdiff;
      b = 2*pC'*pdiff;
      c = pC'*pC - Delta^2;
      if ( a == 0 )
	theta = min(-b/c, 1);
      else
	theta = min(-2*c/(b+sqrt(b^2-4*a*c)), 1);
	theta = max(theta, 0);  % ensure theta >= 0
      end
      p = pC + theta*pdiff;
      if trace > 0
      fprintf('dogleg_tr: Using partial Newton step: theta = %g\n', theta);
      end
    end
  end
end

end
