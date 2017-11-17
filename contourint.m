function sum = contourint(m,dist,fun)
tic();
theta = 2*pi/m;

sum = 0;
for j=1:m
%     temp1 = cos(j*theta)-cos((j-1)*theta);
%     temp2 = 1i*(sin(j*theta)-sin((j-1)*theta));
%     deltax = dist*(temp1 + temp2);
    z1 = dist*cos((j-1)*theta) + 1i*dist*sin((j-1)*theta);
    z2 = dist*cos(j*theta) + 1i*dist*sin(j*theta);
    deltax = (z2-z1);
    sum = sum + 0.5*deltax*(fun(z1) + fun(z2));
end

%toc()
end

