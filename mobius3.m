function mobius3(n,m,b,sd,radius,R)
colors

X = unifrnd(-10,10,n,1) + 1i.*unifrnd(-10,10,n,1);
sigma1 = max(real(X));
sigma = sigma1 + radius;

Contour1 = sigma + 1i.*normrnd(0,sd,m,1);
Contour2 = sigma + 1i.*unifrnd(-R,R,m,1);

figure
h1=plot(X,'.','MarkerSize',10,'Color',Color(:,9));
hold on
h2=plot([Contour1;Contour2],'.','MarkerSize',7,'Color',Color(:,11));
axis([-13 13 -13 13])
xL = xlim;
yL = ylim;
line([0 0], yL,'Color','k');  %y-axis
line(xL, [0 0],'Color','k');  %x-axis
title('Bromwich Contour plus Singularities')
xlabel('real part s')
ylabel('imaginary part s')
legend([h1(1),h2(1)],'singularities of F(s)','contour')
print('bromwichs','-djpeg')
hold off


MX = (X-sigma-b)./(X-sigma+b);
MC1 = (Contour1-sigma-b)./(Contour1-sigma+b);
MC2 = (Contour2-sigma-b)./(Contour2-sigma+b);


figure
h1=plot(MX,'.','MarkerSize',10,'Color',Color(:,9));
hold on
h2=plot([MC1;MC2],'.','MarkerSize',7,'Color',Color(:,11));
axis([-1.2 1.6 -1.12 1.12])
xL = xlim;
yL = ylim;
line([0 0], yL,'Color','k');  %y-axis
line(xL, [0 0],'Color','k');  %x-axis
title('Mobius Transformation of Bromwich Contour and Singularities')
xlabel('real part s')
ylabel('imaginary part s')
legend([h1(1),h2(1)],'singularities of F(s)','contour')
print('bromwichsm','-djpeg')
hold off


end

