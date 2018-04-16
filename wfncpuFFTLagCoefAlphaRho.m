function LaguerreCoef = wfncpuFFTLagCoefAlphaRho(FLaplace, Ncoeff, alphaP, rhoP)
%WFNCPUFFTLAGCOEFALPHARHO Laguerre polynomial coefficients {a_{n}} computation
%  The function utilizes the midpoint version of the FFT 
%  to compute the Laguerre polynomial expansion coefficients from 
%  the sampled Laplace space function.    
%  This function uses standard MATLAB/CPU computations.  
%
%  Use:
%  LaguerreCoeff = wfncpuFFTLagCoefAlphaRho(FLaplace, Ncoeff, alphaP, rhoP)
%  
%  Input: 
%  FLaplace = a symbolic expression for the Laplace transform space function F(s)
%  Ncoeff = number of Laguerre expansion coefficients, 
%  alphaP = Weeks alpha parameter
%  rhoP = Weeks rho parameter
%
%  Output:
%  LaguerreCoeff = a vector of Laguerre polynomial expansion coefficients  
%  LaguerreCoeff is a standard MATLAB variable.  
%  
%  Comment:
%  By using an FFT, the summation can be performed in O(Nlog(N)) operations 
%  as opposed O(N*N) required by the direct application of the midpoint rule.
%  The coefficients corresponding to negative indices which result from the FFT 
%  are not needed for the summation and thus are neglected. 
%  To use the FFT, the number of F(s) samples is twice Ncoeff.
%
%  Author: 
%  Patrick Kano, Moysey Brio - 2016
%
%  Modification Date [M/D/Y]:
%  03/04/2016 - Version 1.0

FFTSamples = zeros(1,2*Ncoeff,'double'); %Twice the samples as the number of coefficients 

%
jidxvec = -Ncoeff:(Ncoeff-1);

Wtemp = exp(1i*(jidxvec+1/2)*pi/Ncoeff);

s = (rhoP./(1-Wtemp)) - alphaP;

Gval = eval(FLaplace); %FLaplace is an expresion in terms of s

FFTSamples(jidxvec+Ncoeff+1) = (rhoP./(1-Wtemp)).*Gval;
 
%Note the order: FFTSamples(1,1:2*N)
TempCoef = fftshift(fft(fftshift(FFTSamples)))/(2*Ncoeff); 

%Use only part of the TempCoef
LaguerreCoef = TempCoef(Ncoeff+1:2*Ncoeff).*exp(-1i*(0:Ncoeff-1)*pi/(2*Ncoeff));

end %function definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
