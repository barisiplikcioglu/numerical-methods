function [res,err] = tanh_sinh_quad(func,lim1,lim2,tol)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% 
% A basic application of the tanh-sinh (double exponential) quadrature on MATLAB
%
% Author: Suat Barış İplikçioğlu, MSc.
%         Department of Electrical and Electronics Engineering
%         Koç University, Istanbul, Turkey
% 
% INPUTS: 
% * func:   Function to be integrated. If the function has singularities in 
%           the integration domain, the integral should be split and these
%           should be in the integral boundaries.
% * lim1:   Lower bound of the integral.
% * lim2:   Upper bound of the integral.
% * tol:    Tolerance in terms of number of significant digits (10^-tol).
%
% OUTPUT:
% * res:    Quadrature result
% * err:    Relative error in between last two iterations
%
% NOTES:
% * This function applies the double exponential numerical integration 
%   scheme introduced by H. Takahasi and M. Mori in 1974.
% * Integration step sizes ('h') and number of summed terms are 
%   are approximated at each step. The routine terminates when the 
%   difference between two successive step size iterations are less
%   than 10^(-tol).
% * Singularities in the integration domain should be at the bounds. 
%
% SOURCES:
% *  H. Takahasi and M. Mori, “Double exponential formulas for numerical
%    integration,” Publications of RIMS, Kyoto University, vol. 9 (1974), 
%    pg. 721–741.
% *  D. H. Bailey, "Tanh-Sinh High-Precision Quadrature" (2006)
% *  D. H. Bailey. K Jeyabalan. X. S. Li. "A comparison of three 
%    high-precision quadrature schemes," Experiment. Math. 14 (3) (2005), 
%    p. 317-329
% *  Truncation formula by: https://www.hpmuseum.org/forum/thread-7690.html
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
tol=ceil(tol);
if((tol<1))
    error('"tol" should be a positive integer.');
end

bf1=(lim1+lim2)/2;
bf2=(-lim1+lim2)/2;

hx=2;
k=0;
itmax=16;
told=10^(-tol/2);

while(1)
   h=hx^-k;
   ri=@(t) tanh(pi/2*sinh(h.*t));
   wi=@(t) pi/2.*cosh(h.*t)./cosh(pi./2.*sinh(h.*t)).^2;
   
   N=abs(ceil(1/h.*acosh(1/pi.*log(2*10^tol.*min([1,bf2])))));
   
   fvx=-N:1:N;
   
   iappr=bf2.*h.*sum(func(bf1+bf2.*ri(fvx)).*wi(fvx));
   if(k==0)
       tolerr=100;
   else
       tolerr=abs(itemp-iappr);
   end
   
   itemp=iappr;
   k=k+1;
   
   if( (tolerr<told) || k>itmax )
       break;
   end 
end

res=iappr;
err=tolerr;

end
