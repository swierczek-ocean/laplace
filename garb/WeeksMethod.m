function varargout = WeeksMethod(FLaplace,Timevec,RelErrorTol,SearchSwitch)
%WEEKSMETHOD A Weeks' method numerical inverse Laplace transform of F(s) to f(t) 
%
%  Use:
%  Two Output Formats are avaiable:
%
%  Version 1.0 Format)
%  If SearchSwitch<2 (sigma,b) search
%  [Invertf,sigmaP,rhoP,RelTotalError,AbsTotalError,AbsTruncateError,AbsRoundoffError,ToleranceMetFlag,LaguerreCoef] =...
%  WeeksMethod(FLaplace,Timevec,RelErrorTol,SearchSwitch)
%  else (alpha,rho) search
%  [Invertf,sigmaP,rhoP,RelTotalError,AbsTotalError,AbsTruncateError,AbsRoundoffError,ToleranceMetFlag,LaguerreCoef] =...
%  WeeksMethod(FLaplace,Timevec,RelErrorTol,SearchSwitch)
%  end
%
%  
%  Version 2.0 Format) 
%  [Invertf,OutParamS,OutErrorS,LaguerreCoef]=
%  WeeksMethod(FLaplace,Timevec,RelErrorTol,SearchSwitch)
%  
%  Input:
%  Of the 4 inputs listed, the first 2 are required.  The last 2 are optional.
%
%  At least the Laplace space function and time(s) at which to invert must be provided:
%  1)FLaplace = a character string expression for the Laplace transform space
%  function F(s) in terms of 's'  
%  2) Timevec = times where f(t) are to be computed (vector)
%  
%  Optional are the search methods parameters, 
%  note that options 0,1,2 are searches over (sigma,b) parameters
%  while 3,4 are searches over (alpha,rho) parameters
%  and 5 is a 1D search over only alpha with (rho=1, equivalently b=1/2)    
%  
%  RelErrorTol (optional) = an input total relative error tolerance for each time [default 5 percent]
%  
%  SearchSwitch (optional) = the search method [default 0]
%   0 - CPU local fminbnd search over sigma and b
%   1 - CPU global search over sigma and b
%   2 - CPU local fmisearch 2D
%   3 - CPU local fminbnd search over alpha and rho
%   4 - CPU global search over alpha and rho
%   5 - CPU global 1D search over only alpha with adaptive integration
%
%  Output: 
%  (all outputs have a length equal to that of the time input vector (Timevec)
%
%  Format version 1.0)
%  Invertf = numerically computed f(t)
%
%  sigmaP = automatically determined sigma or the user defined value if an input 
%  bP = automatically determined beta or the user defined value if an input
%
%  RelTotalError = the estimated total relative error
%  AbsTotalError = the estimated total absolute error
%  AbsTruncateError = the estimated absolute trunction error
%  AbsRoundoffError = the estimated absolute round-off error
%  ToleranceMetFlag = this flag indicates if the relative error tolerance was met: 0-fail, 1-success
%
%  LaguerreCoef = the Laguerre polynomial expansion coefficients
%
%  Format version 2.0)
%  Invertf = numerically computed f(t)
%  
%  OutParamS = array of structures with elements 
%    (the parameters are either automatically determined or the user defined values)
%    sigma = (rho/2) - alpha;
%    b     = rho/2
%    alpha = b - sigma; 
%    rho   = 2b;
%    MoebiusMP = Moebius transform parameters (determinant is NOT scaled to 1)
%              = [1 -sigma-b; 1 -sigma+b]
%              = [1 alpha-rho; 1 rho]
%    Moebiusdet = Moebius transform determinant
%               = 2b
%               = rho
%    Moebiustrace = Moebius transform trace
%                 = 1-sigma+b
%                 = 1+alpha
%
%  OutErrorS = array of sturctures with elements
%    RelTotalError = the estimated total relative error
%    AbsTotalError = the estimated total absolute error
%    AbsTruncateError = the estimated absolute trunction error
%    AbsRoundoffError = the estimated absolute round-off error
%    ToleranceMetFlag = this flag indicates if the relative error tolerance was met: 0-fail, 1-success
%
%  LaguerreCoef = the Laguerre polynomial expansion coefficients
%
%  Comment:
%  The Weeks method is typically expressed in terms of expansion parameters (sigma,b).
%
%  This corresponds to a Moebius transform with of the form
%  w=(s-sigma-b)/(s-sigma+b) with the following Laguerre polynomial expansion expresion for F(s):
%  s = \sigma - b\frac{w+1}{w-1} = \sigma -b + \frac{2b}{1-w}
%  \sum_{n=0}^{\infty} a_{n}w^{n} = \frac{2b}{1-w}F(s)
%
%  The time domain function f(t) is then:
%  f(t) = \exp(\sigma t)\sum_{n=0}^{\infty}a_{n} exp{(-bt)} L_{n}(2bt)
%
%  where the Laguerre polynomials are defined by:
%  L_{n}(x) = \frac{e^{x}}{n!} \frac{d^{n}}{dx^{n}}(x^{n} e^{-x})
%
%  The Moebius transform however can be also be expressed in a form with a
%  more intutitive interpretation as a composition of transfrom from the
%  s to the w plane.
%  w = \beta + \rho*e^{i\theta)}/(s+\alpha)
%  where \alpha,\beta,\rho, and \theta are all real.
%  Or for the Weeks method, 
%  w = 1 + -2b/(s+(b-\sigma)
%  \beta = 1
%  \rho = 2b
%  \alpha = b-\sigma
%  \theta = \pi
%
%  In terms of alpha and rho, then
%  \sum_{n=0}^{\infty} a_{n}w^{n} = \frac{\rho}{1-w}F(\frac{\rho}{1-w} - \alpha)
%  f(t) = \exp(-\alpha t)\exp((\rho/2)t)\sum_{n=0}^{\infty}a_{n} exp{(-(\rho/2)t)} L_{n}(\rho t)
%  f(t) = \exp(-\alpha t)\sum_{n=0}^{\infty}a_{n}L_{n}(\rho t)
%  L_{n}(x) = \frac{e^{x}}{n!} \frac{d^{n}}{dx^{n}}(x^{n} e^{-x})
%
%  Note that although the (\sigma,b) and (\alpha,\rho) representations are
%  algebrically equivalently, numerical optimization over the two parameter
%  spaces typically yields different approximations. This can be seen as another 
%  manifestation of the well known sensitivity of numerical Laplace transform inversion.
%
%  Finally, this code already returns the Moebius transform coefficients in the standard form
%  W=(m11*Z+m12)/(m21*Z+m22)
%  Using the fact that 3 points in the Z and W planes are sufficient to 
%  determine m11,m12,m21,and m22, and Zj->Wj as the 3 known mapped points, then 
%  m11 = w3*(w2-w1)*(z2-z3) - w1*(w2-w3)*(z2-z1)
%  m12 = w1*z3*(w2-w3)*(z2-z1) - w3*z1*(w2-w1)*(z2-z3)
%  m21 = (w2-w1)*(z2-z3) - (w2-w3)*(z2-z1)
%  m22 = z3*(w2-w3)*(z2-z1) - z1*(w2-w1)*(z2-z3)
%  For example, z1=0->w1=i, z2=1->w2=2, z3=-1->w3=4
%  m11 = (16-6i)
%  m12 = 2i
%  m21 = (6-2i)
%  m22 = 2
%  or, if normalized by sqrt(m11*m22 - m12*m21) = sqrt(84-72i)
%  m11 \approx 1.6246 - 0.0072i
%  m12 \approx -0.0660 + 0.1783i
%  m21 \approx 0.6010 + 0.0196i
%  m22 \approx 0.1783 + 0.0660i
%
%  References: 
%  Weeks' method:
%  1) Algorithms for Parameter Selection in the Weeks Method for Inverting the Laplace Transform, 
%     J. A. C. Weideman, SIAM Journal on Scientific Computing, vol. 21, pp. 111-128, 1999. 
%     http://dip.sun.ac.za/~weideman/research/weeks.html
%  2) Application of Weeks method for the numerical inversion of
%     the Laplace transform to the matrix exponential, P. Kano, M. Brio, and J. Moloney,
%     Communications in mathematical sciences", vol. 3, no. 3, pp. 335-372, 2005.
%     http://math.arizona.edu/~brio/WEEKS_METHOD_PAGE/pkanoWeeksMethod.html
%  3) Software for an Implementation of Weeks' Method for the Inverse Laplace Transform Problem
%     S. Garbow, G. Giunta, and J. N. Lyness, A. Murli, ACM, 1988.
%     http://www.utd.edu/~cantrell/ee6481/lectures/Laplace/Garbow-88.pdf
%  4) Moebius Transformations Revealed
%     Douglas N. Arnold and Jonathan Rogness
%     Notices of the AMS, Vol. 55, No. 11, p. 1226-1231, Nov. 2008.
%     http://www.ima.umn.edu/~arnold/moebius/
%  5) The Geometry of Moebius Transformations
%     John Olsen, Spring 2010 Notes,
%     www.johno.dk/mathematics/moebius.pdf
% 
%  Author: 
%  Patrick Kano, Moysey Brio - 2011, 2016
%  Modification Date [M/D/Y]:
%  03/31/2011 - Version 1.0
%  06/10/2016 - Version 2.0 with adapative alpha integration

%% Initialize Output Vectors
Ntimes = length(Timevec);
Invertf = zeros(1,Ntimes);
sigmaP = zeros(1,Ntimes);
bP = zeros(1,Ntimes);
alphaP = zeros(1,Ntimes);
rhoP = zeros(1,Ntimes);
RelTotalError = zeros(1,Ntimes);
AbsTotalError = zeros(1,Ntimes);
AbsTruncateError = zeros(1,Ntimes);
AbsRoundoffError = zeros(1,Ntimes);
ToleranceMetFlag = zeros(1,Ntimes); %0-failed to meet tolerance, 1-succeeded in meet tolerance
LaguerreCoef = cell(1,Ntimes); %this can be a large amount of space 

%Define default (sigma, b) space search parameters
sigmin = 0;
sigmax = 20;
tols = [1.0 0.5 0.25]; %grid resolution

bmax = 40;
tolb = [1.0 0.5 0.25]; %grid resolution

%Define default (rho, alpha) space search parameters
alphamax = 5.0;
alphamin = -5.0; %-1.0*alphamax;
tolalpha = [0.1 0.05 0.025]; %grid resolution

rhomax = 2*alphamax;
tolrho = [0.1 0.05 0.025]; %grid resolution

%Define the default number of Laguerre expansion coefficients
%NLag = [32 64 128 256 512]; %very large number of coefficients is fine for
%FFT but will not work with the adaptive integration
NLag = 8;


switch nargin
 case 2
  SearchSwitch = 0; %default local fminbnd search
  RelErrorTol = 0.05*ones(1,Ntimes);
 case 3
  SearchSwitch = 0; %default local fminbnd search
 case 4
  %do nothing, the user has specified all input parameters to the wrapper
 otherwise
  error('WeeksMethod:nargincheck','Incorrect number of inputs.  At least 2 are required.');
end %switch 

%%%Main Search
for tidx=1:Ntimes
 %The results will increase with accuracy until the tolerance is met   
 for ntolidx=1:length(tols)
  for nLagidx=1:length(NLag)
   %tempindices = [tidx,ntolidx,nLagidx]
   
   if(SearchSwitch==0 || SearchSwitch==1) %sigma,b searches
   
    [testInvertf,testsigmaP,testbP,testRelTotal,testAbsTotal,testAbsTruncate,testAbsRoundoff,testLaguerreCoef]=...  
    wfnWeeksCoreSigmab(FLaplace,Timevec(tidx),NLag(nLagidx),sigmin,sigmax,bmax,SearchSwitch,tols(ntolidx),tolb(ntolidx));

    Invertf(tidx)          = testInvertf;
    RelTotalError(tidx)    = testRelTotal;
    AbsTotalError(tidx)    = testAbsTotal;
    AbsTruncateError(tidx) = testAbsTruncate;
    AbsRoundoffError(tidx) = testAbsRoundoff;
    LaguerreCoef{tidx}     = testLaguerreCoef;

    sigmaP(tidx) = testsigmaP;
    bP(tidx)     = testbP;
    
    alphaP(tidx) = bP(tidx) - sigmaP(tidx);
    rhoP(tidx)   = 2.0*bP(tidx);
        
    if testRelTotal<RelErrorTol(tidx)
     ToleranceMetFlag(tidx) = 1;
     break;
    else %if(RelErrorTol(tidx)<testRelTotal) %update
     ToleranceMetFlag(tidx) = 0;
    end
    
   elseif(SearchSwitch==2 || SearchSwitch==3 || SearchSwitch==4) %SearchSwitch alpha, rho search
   
    [testInvertf,testalphaP,testrhoP,testRelTotal,testAbsTotal,testAbsTruncate,testAbsRoundoff,testLaguerreCoef]=...  
    wfnWeeksCoreAlphaRho(FLaplace,Timevec(tidx),NLag(nLagidx),alphamin,alphamax,rhomax,SearchSwitch,tolalpha(ntolidx),tolrho(ntolidx));

    Invertf(tidx)          = testInvertf;    
    RelTotalError(tidx)    = testRelTotal;
    AbsTotalError(tidx)    = testAbsTotal;
    AbsTruncateError(tidx) = testAbsTruncate;
    AbsRoundoffError(tidx) = testAbsRoundoff;
    LaguerreCoef{tidx}     = testLaguerreCoef;
    
    alphaP(tidx) = testalphaP;
    rhoP(tidx)   = testrhoP;

    sigmaP(tidx) = 0.5*rhoP(tidx) - alphaP(tidx);
    bP(tidx)     = 0.5*rhoP(tidx);
    
    if testRelTotal<RelErrorTol(tidx)
     ToleranceMetFlag(tidx) = 1;
     break;
    else %if(RelErrorTol(tidx)<testRelTotal) %update
     ToleranceMetFlag(tidx) = 0;
    end
   
   elseif(SearchSwitch==5) %search alpha with adaptive integration
    %Unlike the other methods, time is not considered in the optimization
    %function used for determining the optimal parameter
       
    disp('Selected Method: search alpha with adaptive integration');
    if(tidx==1)
     alphamin = -50.0;
     alphamax = 5%0.5;
     tolalpha = 0.01;
    
     Fofs = @(s)eval(FLaplace); %convert from a string to a function handle
     [testInvertf,testalphaP,testrhoP,testRelTotal,testAbsTotal,testAbsTruncate,testAbsRoundoff,testLaguerreCoef]=...  
     wfnWeeksCoreAdaptiveIntegrate(Fofs,Timevec(tidx),NLag(nLagidx),alphamin,alphamax,tolalpha);
    else
     %just use the values 
     [testInvertf,testalphaP,testrhoP,testRelTotal,testAbsTotal,testAbsTruncate,testAbsRoundoff,testLaguerreCoef]=...  
     wfnWeeksCoreAdaptiveIntegrate(Fofs,Timevec(tidx),NLag(nLagidx),testalphaP);
    end
    
    Invertf(tidx)          = testInvertf;    
    RelTotalError(tidx)    = testRelTotal;
    AbsTotalError(tidx)    = testAbsTotal;
    AbsTruncateError(tidx) = testAbsTruncate;
    AbsRoundoffError(tidx) = testAbsRoundoff;
    LaguerreCoef{tidx}     = testLaguerreCoef;
    
    alphaP(tidx) = testalphaP;
    rhoP(tidx)   = testrhoP;

    sigmaP(tidx) = 0.5*rhoP(tidx) - alphaP(tidx);
    bP(tidx)     = 0.5*rhoP(tidx);
    
    if testRelTotal<RelErrorTol(tidx)
     ToleranceMetFlag(tidx) = 1;
     break;
    else %if(RelErrorTol(tidx)<testRelTotal) %update
     ToleranceMetFlag(tidx) = 0;
    end
       
   else
    error([mfilename,':SearchSwitch'],'Incorrect choice of the search switch.');
    
   end %SearchSwitch 
  end %nLagidx
  
  if(ToleranceMetFlag(tidx)==1), break; end %out to next time index
 end %ntolidx
end %tidx    

%% Summarize the outputs into the appropriate format
if nargout==4 
%  
%  Version 2.0 Format) 
%  [Invertf,OutParamS,OutErrorS,LaguerreCoef]=
%  WeeksMethod(FLaplace,Timevec,RelErrorTol,SearchSwitch)

 varargout{1} = Invertf;
 
 for tidx=1:Ntimes
  
  OutParamS(tidx).sigmaP = sigmaP(tidx);
  OutParamS(tidx).bP     = bP(tidx);

  OutParamS(tidx).alpha = alphaP(tidx); %bP(tidx)-sigmaP(tidx);
  OutParamS(tidx).rho   = rhoP(tidx); %2.0*bP(tidx);  

  %Zvec = [0 1 -1];
  %Wvec(1) = (Zvec(1) + sigmaP(tidx) - bP(tidx))/(Zvec(1) + sigmaP(tidx) + bP(tidx));
  %Wvec(2) = (Zvec(2) + sigmaP(tidx) - bP(tidx))/(Zvec(2) + sigmaP(tidx) + bP(tidx));
  %Wvec(3) = (Zvec(3) + sigmaP(tidx) - bP(tidx))/(Zvec(3) + sigmaP(tidx) + bP(tidx));  
   
  % MoebiusMP = Moebius transform parameters (determinant is NOT scaled to 1)
%              = [1 -sigma-b; 1 -sigma+b]
%              = [1 alpha-rho; 1 rho]
  OutParamS(tidx).MoebiusMP(1,1) = 1;
  OutParamS(tidx).MoebiusMP(1,2) = -sigmaP(tidx) - bP(tidx);
  OutParamS(tidx).MoebiusMP(2,1) = 1;
  OutParamS(tidx).MoebiusMP(2,2) = -sigmaP(tidx) + bP(tidx);    

  % Moebiusdet = Moebius transform determinant
  %            = 2b
  %            = rho  
  OutParamS(tidx).Moebiusdet = 2*bP(tidx);
 
  % Moebiustrace = Moebius transform trace
  %              = 1+sigma+b
  %              = 1+alpha 
  OutParamS(tidx).Moebiustrace = 1+sigmaP(tidx)+bP(tidx);
 end %tidx
 varargout{2} = OutParamS;
 
 OutErrorS.RelTotalError    = RelTotalError;
 OutErrorS.AbsTotalError    = AbsTotalError;
 OutErrorS.AbsTruncateError = AbsTruncateError;
 OutErrorS.AbsRoundoffError = AbsRoundoffError;
 OutErrorS.ToleranceMetFlag = ToleranceMetFlag;
 varargout{3} = OutErrorS;
 
 varargout{4} = LaguerreCoef; 
elseif nargout==1 %just return Invertf
 varargout{1} = Invertf;
else
%  Version 1.0 Legacy Format with only sigma,b searches
%  [Invertf,sigmaP,bP,RelTotalError,AbsTotalError,AbsTruncateError,AbsRoundoffError,ToleranceMetFlag,LaguerreCoef] =...
%  WeeksMethod(FLaplace,Timevec,RelErrorTol,SearchSwitch)

%  or (depending on the search switch)
%  [Invertf,alphaP,rhoP,RelTotalError,AbsTotalError,AbsTruncateError,AbsRoundoffError,ToleranceMetFlag,LaguerreCoef] =...
%  WeeksMethod(FLaplace,Timevec,RelErrorTol,SearchSwitch)

  varargout{1} = Invertf;
  if(SearchSwitch<2)
  varargout{2} = sigmaP;
  varargout{3} = bP;
  else
  varargout{2} = alphaP;
  varargout{3} = rhoP;     
  end
  varargout{4} = RelTotalError;
  varargout{5} = AbsTotalError;
  varargout{6} = AbsTruncateError;
  varargout{7} = AbsRoundoffError;
  varargout{8} = ToleranceMetFlag;
  varargout{9} = LaguerreCoef;
  % error([mfilename,':nargout'],'The user may select 1, 4, or 9 outputs');
end 
end %function definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
