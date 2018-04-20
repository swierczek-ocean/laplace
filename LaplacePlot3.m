function LaplacePlot3(True,NAB,Weeks,SWeeks,T,t1,t2,t3,char,a,b,eps)
colors

Q = sort(True);
Z = sort(Weeks);
W = sort(SWeeks);
N = sort(NAB);
ps = 14;
ss = 5.5;
Nabc = 34;
Wkc = 16;
Swkc = 12;

if(Q(end)-Q(1)>1e5)&&(Q(1)>0)&&(Z(1)>0)&&(W(1)>0)&&(N(1)>0)
    
    figure
    h1 = semilogy(T,True,'LineWidth',3,'Color',Color(:,35));
    hold on
    h2 = semilogy(t1,real(Weeks),'.','MarkerSize',ps,'Color',Color(:,Wkc));
    semilogy(t1,real(Weeks),'v','MarkerSize',ss,'Color',Color(:,Wkc));
    h3 = semilogy(t2,real(SWeeks),'.','MarkerSize',ps,'Color',Color(:,Swkc));
    semilogy(t2,real(SWeeks),'o','MarkerSize',ss,'Color',Color(:,Swkc));    
    h4 = semilogy(t3,real(NAB),'.','MarkerSize',ps,'Color',Color(:,Nabc));
    semilogy(t3,real(NAB),'d','MarkerSize',ss,'Color',Color(:,Nabc));    
    title(['L',num2str(char),' a=',num2str(a),' b=',num2str(b)])
    xlabel('time')
    ylabel('f(t)')
    legend([h1(1),h2(1),h3(1),h4(1)],'True f(t)','Weeks','Shift Weeks','adaptive','Location','northwest')
    print(['L',num2str(char),' a=',num2str(a),' b=',num2str(b),'_true_vs_ests'],'-djpeg')
    hold off
    
else
    
    figure
    h1 = plot(T,True,'LineWidth',2.75,'Color',Color(:,35));
    hold on
    h2 = plot(t1,real(Weeks),'.','MarkerSize',ps,'Color',Color(:,Wkc));
    plot(t1,real(Weeks),'v','MarkerSize',ss,'Color',Color(:,Wkc));
    h3 = plot(t2,real(SWeeks),'.','MarkerSize',ps,'Color',Color(:,Swkc));
    plot(t2,real(SWeeks),'o','MarkerSize',ss,'Color',Color(:,Swkc));    
    h4 = plot(t3,real(NAB),'.','MarkerSize',ps,'Color',Color(:,Nabc));
    plot(t3,real(NAB),'d','MarkerSize',ss,'Color',Color(:,Nabc));    
    title(['L',num2str(char),' a=',num2str(a),' b=',num2str(b)])
    xlabel('time')
    ylabel('f(t)')
    legend([h1(1),h2(1),h3(1),h4(1)],'True f(t)','Weeks','Shift Weeks','adaptive','Location','northwest')
    print(['L',num2str(char),' a=',num2str(a),' b=',num2str(b),'_true_vs_ests'],'-djpeg')
    hold off
    
end

True1 = master_inverse_laplace_fcn(t1,a,b,char,eps);
ind = find(abs(True1)==0);
WeeksRelError = abs(Weeks-True1)./abs(True1);
WeeksRelError(ind) = abs(Weeks(ind)-True1(ind));

True2 = master_inverse_laplace_fcn(t3,a,b,char,eps);
ind = find(abs(True2)==0);
NABRelError = abs(NAB-True2)./abs(True2);
NABRelError(ind) = abs(NAB(ind)-True2(ind));

True3 = master_inverse_laplace_fcn(t2,a,b,char,eps);
ind = find(abs(True3)==0);
SWeeksRelError = abs(SWeeks-True3)./abs(True3);
SWeeksRelError(ind) = abs(SWeeks(ind)-True3(ind));


figure
h1 = semilogy(t1,real(WeeksRelError),'.','MarkerSize',ps,'Color',Color(:,Wkc));
hold on
semilogy(t1,real(WeeksRelError),'v','MarkerSize',ss,'Color',Color(:,Wkc))
h2 = semilogy(t2,real(SWeeksRelError),'.','MarkerSize',ps,'Color',Color(:,Swkc));
semilogy(t2,real(SWeeksRelError),'o','MarkerSize',ss,'Color',Color(:,Swkc))
h3 = semilogy(t3,real(NABRelError),'.','MarkerSize',ps,'Color',Color(:,Nabc));
semilogy(t3,real(NABRelError),'d','MarkerSize',ss,'Color',Color(:,Nabc))
title(['L ',num2str(char),' a=',num2str(a),' b=',num2str(b),' error'])
xlabel('time')
ylabel('Relative error')
legend([h1(1),h2(1),h3(1)],'Weeks','Shift Weeks','adaptive','Location','northwest')
print(['L',num2str(char),' a=',num2str(a),' b=',num2str(b),'_x_error'],'-djpeg')
hold off

end

