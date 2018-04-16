function LaguerreCoef = wfncpuFFTLagCoefSigmab(FLaplace, Ncoeff, sigmaP, bP)
%WFNCPUFFTLAGCOEFSIGMAB Laguerre polynomial coefficients {a_{n}} computation
%  The function utilizes the midpoint version of the FFT 
%  to compute the Laguerre polynomial expansion coefficients from 
%  the sampled Laplace space function.    
%  This function uses standard MATLAB/CPU computations.  
%
%  Use:
%  LaguerreCoeff = wfncpuFFTLagCoefSigmab(FLaplace, Ncoeff, sigmaP, bP)
%  
%  Input: 
%  FLaplace = a symbolic expression for the Laplace transform space function F(s)
%  Ncoeff = number of Laguerre expansion coefficients, 
%  sigmaP = Weeks sigma parameter
%  bP = Weeks b parameter
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
%  Patrick Kano, Moysey Brio - 2011
%
%  Modification Date [M/D/Y]:
%  03/30/2011 - Version 1.0
%  03/04/2016 - Just a name change

FFTSamples = zeros(1,2*Ncoeff,'double'); %Twice the samples as the number of coefficients 

%
jidxvec = -Ncoeff:(Ncoeff-1);

Wtemp = exp(1i*(jidxvec+1/2)*pi/Ncoeff);

s = sigmaP-bP+(2*bP./(1-Wtemp)); %eqivalent to sigma - b*(w+1)/(w-1)

Gval = eval(FLaplace); %FLaplace is an expresion in terms of s

FFTSamples(jidxvec+Ncoeff+1) = (2*bP./(1-Wtemp)).*Gval;
 
%Note the order: FFTSamples(1,1:2*N)
TempCoef = fftshift(fft(fftshift(FFTSamples)))/(2*Ncoeff); 

%Use only part of the TempCoef
LaguerreCoef = TempCoef(Ncoeff+1:2*Ncoeff).*exp(-1i*(0:Ncoeff-1)*pi/(2*Ncoeff));

end %function definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
