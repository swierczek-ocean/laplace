a=3;
b=5;
eps=0.1;

set(0,'DefaultFigureVisible','on')




% [X,Y] = meshgrid(-10:0.1:10,-10:0.1:10);
% Z = imag(master_laplace_fcn(X+1i.*Y,a,b,125,eps));
% Z(:,145:201) = 0;
% figure
% colormap('autumn');
% surf(X,Y,Z)
% xlabel('x')
% ylabel('iy')
% zlabel('imaginary part')


% Years = [2004,1970,2004,1979,1978,2007,1974,1982,1984,1976,1995,1968,1966,1966,1979,1979,1975,1981,1966,1966,1986,1970,1978,1999,1982,1970,2007,1981,2008,1975,1969,1981,2004,2011,1982,1968,1984,1964,2006,2006,1968,1987,2009,1957,1983,1984,2000,2008,2013,1972,1971,1971,1976,1972,1994,1999,1993,1988,2002,1982,1992,1966];
% 
% figure
% histogram(Years,10,'FaceColor','b')
% xlabel('year')
% ylabel('number of papers')
% title('highly cited papers on numerical ILT by decade')
% print('paper_hist','-djpeg')



% n = 100;
% 
% dt = 0:pi/50:(2*pi-pi/50);
% 
% Theta = exp(1i.*dt);
% 
% figure
% plot(Theta,'.','MarkerSize',8)
% axis([-1.5 1.5 -1.5 1.5])
% xlabel('x')
% ylabel('iy')
% title('equidistant partition of unit circle')
% 


% load('error1_100.mat')
% error1 = error;
% load('error111_152.mat')
% error2 = error;
% load('error153_235.mat')
% error3 = error;
% error_total = [error1;error2;error3];
% vector = error_total(:,1:end-1);
% vector = reshape(vector,225*57,1);
% vector = log(vector);
% figure
% histogram(vector,40)
% axis([-40 5 0 2800])
% title('adaptive bromwich distribution of errors')
% xlabel('log error')
% ylabel('count')
% print('NAB_error_dist','-djpeg')




% load('errorW1_100.mat')
% error1 = errorW;
% load('errorW111_153.mat')
% error2 = errorW;
% load('errorW153_235.mat')
% error3 = errorW;
% error_total = [error1;error2;error3];
% vector = error_total(:,1:end-1);
% vector = reshape(vector,225*57,1);
% vector = log(vector);
% figure
% histogram(vector,40)
% axis([-40 5 0 2800])
% title('Weeks 0 distribution of errors')
% xlabel('log error')
% ylabel('count')
% print('Weeks_0_error_dist','-djpeg')








