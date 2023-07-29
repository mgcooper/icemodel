function Tsfc = SEBSOLVE(Tair,Qsi,Qli,ea,albedo,De,Pa,wspd,cv_air,...
   emiss,SB,Tf,Qc,xTsfc,chi,roL,scoef,fopts,liqflag,metiter)
%SEBSOLVE solve the surface energy balance for the skin temperature

debug = false;
dflag = true;

% gather terms in the SEB equation. Note that if icond_flag == 0, Qc = 0
AAA = cv_air * De;                           % [W m-2 K-1]
CCC = 0.622 / Pa;                            % [Pa-1] = [m3 J-1]
EEE = chi*(1.0-albedo)*Qsi + emiss*Qli + Qc; % [W m-2]
FFF = roL * De;                              % [W m-2]

% these replace the block above
B1 = scoef(2)/(Tair*wspd^2);
B2 = scoef(3)/(sqrt(Tair)*wspd);

Sfnc = @STABLEFN;
Vfnc = @VAPOR;
fSEB = @(Tsfc) EEE - emiss*SB.*Tsfc.^4 + ...
   AAA.*Sfnc(Tair,Tsfc,wspd,scoef).*(Tair-Tsfc) + ...
   FFF.*CCC.*Sfnc(Tair,Tsfc,wspd,scoef) .* ...
   (ea-Vfnc(Tsfc,Tf,liqflag));

% find the region with zero crossing
[a,b] = fzero_guess_to_bounds(fSEB, xTsfc, xTsfc-50, xTsfc+50);

if isnan(a)
   %try
      % this might return NaN during spinup (doesn't catch), see check at end,
      % but don't use for now, need method to deal with spinup when xTsfc-Tsfc
      % is very large and leads to bad solution
   %   Tsfc = fzero(fSEB,xTsfc,fopts);
   %catch
      % Tsfc = Tair;
      Tsfc = min(Tair,Tf);
   %end
else
   try
      Tsfc = fzero_brent(fSEB, a, b, fopts.TolX);
   catch
      % fallback to fsolve
      try
         Tsfc = fzero(fSEB,[a b],fopts);
      catch
         % fallback to secant method (or fsolve? )
         try
            Tsfc = fzero_secant(fSEB, a, b, fopts.TolX);
         catch
            % Tsfc = Tair; % set to air temperature if all methods fail
            Tsfc = min(Tair,Tf);
         end
      end
   end
end

% This should only occur during spinup, and only if built-in fzero is used
if isnan(Tsfc)
   Tsfc = min(Tair,Tf);
   % Tsfc = Tair;
end


function [root, iter] = fzero_secant(f, a, b, tol,varargin)
iter = 0;
max_iter = 2000;
fa = f(a);
fb = f(b);
while iter < max_iter
   iter = iter + 1;
   c = b - fb * (b - a) / (fb - fa);
   fc = f(c);
   if abs(fc) < tol || abs(c - b) < tol
      break;
   end
   a = b; fa = fb; b = c; fb = fc;
end
root = c;

function b = fzero_brent(f, a, b, tol, varargin)
%FZERO_BRENT  Find a root of a univariate function within a given interval
%             using Brent's method
%
% x = fzero_brent(f,a,b,t)
% finds x in the interval [a,b] satisfying |x-y| <= t/2 where f(y) = 0.
% f(a) and f(b) must have opposite signs.
%
% ... = fzero_brent(f,lb,ub,t,x,...) passes additional inputs directly
% to f.
%
% This function is compatible with MATLAB's code generation -- so long f is
% similarly compatible. Many root-finding problems can be solved with
% fzero_brent by writing another function which calls fzero_brent inside a
% for loop, and this can be made fast by using code generation on that
% wrapper function. Note that f will only be called with a scalar input as
% its first argument; codegen knows this, and might strip out unnecessary
% code from the function definition underlying f.
% --- Example:
% Simple bisection between bounds of -0.5 and 1.5 would fail to find the
% root. Starting with a guess of .85 and expanding outwards finds a root.
% This example is shown graphically on the MATLAB File Exchange.
%
% c = poly([0 1]); % a polynomial with roots at 0 and 1
% f = @(x) polyval(c,x); % encapsulate extra parameters
% [a,b] = fzero_guess_to_bounds(f, .85, -.5, 1.5);
% root = fzero_brent(f, a, b, .05);
%
%  Discussion:
%
%    The interval [A,B] must be a change of sign interval for F.
%    That is, F(A) and F(B) must be of opposite signs.  Then
%    assuming that F is continuous implies the existence of at least
%    one value C between A and B for which F(C) = 0.
%
%    The location of the zero is determined to within an accuracy
%    of 6 * EPS * abs ( C ) + 2 * T, where EPS is the machine epsilon.
%
%    Thanks to Thomas Secretin for pointing out a transcription error in the
%    setting of the value of P, 11 February 2013.
%
%    Additional parameters given by varargin are passed directly to F.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    14 September 2019
%    30 June 2020 - Geoff Stanley
%
%  Author:
%
%    Original FORTRAN77 version by Richard Brent.
%    MATLAB version by John Burkardt.
%    Minor changes (passing extra arguments) by Geoff Stanley.
%
%  Reference:
%
%    Richard Brent,
%    Algorithms for Minimization Without Derivatives,
%    Dover, 2002,
%    ISBN: 0-486-41998-3,
%    LC: QA402.5.B74.
%
%  Parameters:
%
%    Input, real A, B, the endpoints of the change of sign interval.
%
%    Input, real T, a positive error tolerance.
%
%    Input, real value = F ( x ), the name of a user-supplied
%    function which evaluates the function whose zero is being sought.
%
%    Output, real VALUE, the estimated value of a zero of
%    the function F.
%

% Matt Cooper, for compatibility with fzero calling syntax
% a = ab(1);
% b = ab(2);

fa = f(a, varargin{:});
fb = f(b, varargin{:});

c = a;
fc = fa;
e = b - a;
d = e;

while ( true )

   if ( abs ( fc ) < abs ( fb ) )

      a = b;
      b = c;
      c = a;
      fa = fb;
      fb = fc;
      fc = fa;

   end

   tol = 2.0 * eps * abs ( b ) + tol;
   m = 0.5 * ( c - b );

   if ( abs ( m ) <= tol || fb == 0.0 )
      break
   end

   if ( abs ( e ) < tol || abs ( fa ) <= abs ( fb ) )

      e = m;
      d = e;

   else

      s = fb / fa;

      if ( a == c )

         p = 2.0 * m * s;
         q = 1.0 - s;

      else

         q = fa / fc;
         r = fb / fc;
         p = s * ( 2.0 * m * q * ( q - r ) - ( b - a ) * ( r - 1.0 ) );
         q = ( q - 1.0 ) * ( r - 1.0 ) * ( s - 1.0 );

      end

      if ( 0.0 < p )
         q = - q;
      else
         p = - p;
      end

      s = e;
      e = d;

      if ( 2.0 * p < 3.0 * m * q - abs ( tol * q ) && p < abs ( 0.5 * s * q ) )
         d = p / q;
      else
         e = m;
         d = e;
      end

   end

   a = b;
   fa = fb;

   if ( tol < abs ( d ) )
      b = b + d;
   elseif ( 0.0 < m )
      b = b + tol;
   else
      b = b - tol;
   end

   fb = f(b, varargin{:});

   if ( ( 0.0 < fb && 0.0 < fc ) || ( fb <= 0.0 && fc <= 0.0 ) )
      c = a;
      fc = fa;
      e = b - a;
      d = e;
   end

end

function [a, b] = fzero_guess_to_bounds(f, x, A, B, varargin)
%FZERO_GUESS_TO_BOUNDS  Search for a sign change bounding a zero of a
%                       univariate function, expanding geometrically
%                       outward from an initial guess.
%
%
% [a, b] = fzero_guess_to_bounds(f, x)
% finds a < b such that f(a) and f(b) have different sign*, meaning a
% solution exists within the interval [a,b].  The bounds a,b are expanded
% outward in geometric progression from an initial guess for the root of f
% at x. If f evaluates to NaN at any point during the search, then a = nan
% and b = nan are immediately returned.  If the function is genuinely
% single-signed, or even if it is not but its values of opposite sign are
% skipped over, it is possible to enter an infinite loop.  Calling the
% function in this form is therefore not recommended unless you know the
% function will not result in such an infinite loop.
%
% [a, b] = fzero_guess_to_bounds(f, x, A, B)
% as above, but limits [a,b] to remain inside the subset [A, B].  If x is
% outside of [A, B], it is immediately moved into this range. If no
% sign-change is found within [A, B], then a = nan and b = nan are
% returned.  Note, as above, it is possible that a sign-change is skipped
% over as f is only evaluated at finitely many x values.
%
% [a,b] = fzero_guess_to_bounds(f, x, A, B, ...)
% passes all additional arguments to the function f.
%
% * Note: for computational speed, herein the "sign" of 0 is considered the
% same as the sign of a negative number.
%
% This function is compatible with MATLAB's code generation.
%
%
% --- Input:
%   f       : handle to a function that accepts a real scalar as its first
%             input and returns a real scalar
%   x       : initial scalar guess for a root of f
%   A       : scalar lower bound
%   B       : scalar upper bound
%   varargin: All additional inputs are passed directly to f
%
%
% --- Output:
%   a : lower bound for interval containing a root, scalar
%   b : upper bound for interval containing a root, scalar
%
%
% --- Acknowledgements:
% Expansion from initial guess inspired by MATLAB's fzero.m.
%
%
% Author    : Geoff Stanley
% Email     : geoffstanley@gmail.com
% Version   : 1.0
% History   : 01/07/2020 - initial release
%           : 22/07/2020 - fix infinite loop in bounded case, arising from machine precision rounding


% Matt Cooper, for compatibility with fzero syntax NEED TO REVISIT WHEN I HAVE
% MORE TIME TO FIX THE PARSING
% A = AB(1);
% B = AB(2);

% Geometrically expand from the guess x, until a sign change is found
sqrttwo = 1.414213562373095;


if nargin >= 4 && isscalar(A) && isscalar(B)

   % Handle bad inputs
   if isnan(A) || isnan(B) || isnan(x)
      a = nan;
      b = nan;
      return
   end

   fa = f(A, varargin{:});
   if fa == 0
      a = A;
      b = A;
      return
   end

   fb = f(B, varargin{:});
   if fb == 0
      a = B;
      b = B;
      return
   end

   x = min(max(x, A), B);

   % bounds are given
   dxp = (B - x) / 50;
   dxm = (x - A) / 50;

   % Set a = x, except when x is so close to A that machine roundoff makes dxm identically 0,
   % which would lead to an infinite loop below.  In this case, set a = A.
   if dxm == 0
      a = A;
   else
      a = x;
   end
   fapos = f(a, varargin{:}) > 0;

   % Similarly, set b = x, except for machine precision problems.
   if dxp == 0
      b = B;
      fbpos = f(b, varargin{:}) > 0;
   else
      b = x;
      if dxm == 0
         fbpos = fapos; % since a = b = x
      else
         fbpos = f(b, varargin{:}) > 0;
      end
   end

   while true
      if a > A
         % Move a left, and test for a sign change
         dxm = sqrttwo * dxm;
         a = max(x - dxm, A);
         fapos = f(a, varargin{:}) > 0;
         if xor(fapos, fbpos) % fa and fb have different signs
            return
         end
      elseif b == B % also a == A, so cannot expand anymore
         if xor(fapos, fbpos) % one last test for sign change
            return
         else % no sign change found
            a = nan;
            b = nan;
            return
         end
      end

      if b < B
         % Move b right, and test for a sign change
         dxp = sqrttwo * dxp;
         b = min(x + dxp, B);
         fbpos = f(b, varargin{:}) > 0;
         if xor(fapos, fbpos) % fa and fb have different signs
            return
         end
      elseif a == A % also b == B, so cannot expand anymore
         if xor(fapos, fbpos) % one last test for sign change
            return
         else % no sign change found
            a = nan;
            b = nan;
            return
         end
      end
   end

else

   % Handle bad inputs
   if isnan(x)
      a = nan;
      b = nan;
      return
   end

   % no bounds given
   if x ~= 0
      dx = abs(x) / 50;
   else
      dx = 1/50;
   end

   a = x;
   b = x;
   fb = f(b, varargin{:});
   fa = fb; %#ok<NASGU> % needed for codegen
   while true

      dx = sqrttwo * dx;

      % Move a left, and test for a sign change
      a = x - dx;
      fa = f(a, varargin{:});
      if isnan(fa)  % Outside the valid bounds of the function.  Give up.
         a = nan;
         b = nan;
         return
      elseif xor(fa > 0, fb > 0) % fa and fb have different signs
         return
      end

      % Move b right, and test for a sign change
      b = x + dx;
      fb = f(b, varargin{:});
      if isnan(fb)  % Outside the valid bounds of the function.  Give up.
         a = nan;
         b = nan;
         return
      elseif xor(fa > 0, fb > 0) % fa and fb have different signs
         return
      end
   end
end