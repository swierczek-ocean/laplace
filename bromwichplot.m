colors
a=5;
n=3000;
c=3;

y = @(x)sqrt(a^2-(x-c)^2);

D = linspace(-a,a,2*n)';

X = linspace(c-a,c,n)';
Y1 = zeros(n,1);
Y2 = zeros(n,1);

for i=1:n
    Y1(i) = y(X(i));
    Y2(i) = -y(X(i));
end

C = c.*ones(2*n,1);

V = [C,D;X,Y1;X,Y2];

sing = 5;

W1 = normrnd(0.5,1.5,sing,1);
W2 = normrnd(0,1.5,sing,1);
W = [W1,W2];


figure
h1=plot(V(:,1),V(:,2),'o','MarkerSize',3,'Color',Color(:,8));
hold on
h2=plot(W(:,1),W(:,2),'*','MarkerSize',7,'Color',Color(:,9));
axis([-5 7 -7 7])
xL = xlim;
yL = ylim;
line([0 0], yL,'Color','k');  %y-axis
line(xL, [0 0],'Color','k');  %x-axis
title('Bromwich Contour')
xlabel('real part s')
ylabel('imaginary part s')
legend([h1(1),h2(1)],'contour','singularities of F(s)')
print('bromwich','-djpeg')
hold off



















