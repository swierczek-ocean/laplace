function mobius(m,n,radius,b,sd,R)
tic();

[s1,s2] = size(b);
s = max(s1,s2);
dist = zeros(s1,s2);
X = unifrnd(-10,10,n,1) + 1i.*unifrnd(-10,10,n,1);
sigma1 = max(real(X));
sigma = sigma1 + radius;

Contour1 = sigma + 1i.*normrnd(0,sd,m,1);
Contour2 = sigma + 1i.*unifrnd(-R,R,m,1);

    figure
    plot(X,'o','MarkerSize',5,'Color','red')
    hold on
    plot(Contour1,'o','MarkerSize',5,'Color','blue')
    plot(Contour2,'o','MarkerSize',5,'Color','green')
    title('Complex Plane with Contour and Singularities')
    xlabel('Real Part')
    ylabel('Imaginary Part')
    print('Complex','-djpeg')

for i=1:s
    
    MX = (X-sigma-b(i))./(X-sigma+b(i));
    MC1 = (Contour1-sigma-b(i))./(Contour1-sigma+b(i));
    MC2 = (Contour2-sigma-b(i))./(Contour2-sigma+b(i));

    figure
    plot(MX,'o','MarkerSize',5,'Color','red')
    hold on
    plot(MC1,'o','MarkerSize',5,'Color','blue')
    plot(MC2,'o','MarkerSize',5,'Color','green')
    axis([-1.1 1.2 -1.1 1.2])
    title(['Mobius Transform with b = ', num2str(b(i))])
    xlabel('Real Part')
    ylabel('Imaginary Part')
    print(['Mobius',num2str(i)],'-djpeg')
    hold off
    
    MXdist = abs(MX);
    dist(i) = min(MXdist);
end

    dist = dist-1;
    figure
    plot(b,dist,'*','MarkerSize',6,'Color','red')
    title('Plot of b versus minimum distance of mobius transformed singularity')
    xlabel('b')
    ylabel('distance')
    print('singdist','-djpeg')


toc()
end

