function Output = contourdist(m,n,dx)
tic();

R = dx:dx:n*dx;
Output = zeros(n,1);
theta = 2*pi/m;

for i=1:n
    sum = 0;
    for j=1:m
        temp1 = cos(j*theta)-cos((j-1)*theta);
        temp2 = 1i*(sin(j*theta)-sin((j-1)*theta));
        deltax = R(i)*(temp1 + temp2);
        z1 = R(i)*cos((j-1)*theta) + 1i*R(i)*sin((j-1)*theta);
        z2 = R(i)*cos(j*theta) + 1i*R(i)*sin(j*theta);
        sum = sum + 0.5*deltax*(1/z1 + 1/z2);
    end
    Output(i) = sum;
end


toc()
end

