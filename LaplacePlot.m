function error = LaplacePlot(True,NAB,T,char,ub,a,b)
colors

Q = sort(True);
Z = sort(NAB);
ps = 14;
ss = 5.5;
Nabc = 34;


if(Q(end)-Q(1)>1e5)&&(Q(1)>0)&&(Z(1)>0)
    
    figure
    h1 = semilogy(T,True,'LineWidth',3,'Color',Color(:,35));
    hold on
    h2 = semilogy(T,NAB,'.','MarkerSize',ps,'Color',Color(:,Nabc));
    semilogy(T,NAB,'d','MarkerSize',ss,'Color',Color(:,Nabc));
    title(['Laplace function ',num2str(char),' LOC=',num2str(2*ub),' a=',num2str(a),' b=',num2str(b)])
    xlabel('time')
    ylabel('f(t)')
    legend([h1(1),h2(1)],'True f(t)','Bromwich adapt','Location','northwest')
    print(['L',num2str(char),' ',num2str(2*ub),' a=',num2str(a),' b=',num2str(b),'_approx_v_true'],'-djpeg')
    hold off
    
else
    
    figure
    h1 = plot(T,True,'LineWidth',3,'Color',Color(:,35));
    hold on
    h2 = plot(T,NAB,'.','MarkerSize',ps,'Color',Color(:,Nabc));
    plot(T,NAB,'d','MarkerSize',ss,'Color',Color(:,Nabc));
    title(['Laplace function ',num2str(char),' LOC=',num2str(2*ub),' a=',num2str(a),' b=',num2str(b)])
    xlabel('time')
    ylabel('f(t)')
    legend([h1(1),h2(1)],'True f(t)','Bromwich adapt','Location','northwest')
    print(['L',num2str(char),' ',num2str(2*ub),' a=',num2str(a),' b=',num2str(b),'_approx_v_true'],'-djpeg')
    hold off
    
end

ind = find(abs(True)==0);
BromRelError = abs(NAB-True)./abs(True);
BromRelError(ind) = abs(NAB(ind)-True(ind));
error = BromRelError;

figure
semilogy(T,real(BromRelError),'.','MarkerSize',ps,'Color',Color(:,18))
hold on
semilogy(T,real(BromRelError),'d','MarkerSize',ss,'Color',Color(:,Nabc))
title(['Laplace function ',num2str(char),' error,',' LOC=',num2str(2*ub),' a=',num2str(a),' b=',num2str(b)])
xlabel('time')
ylabel('Relative error')
print(['L',num2str(char),' ',num2str(2*ub),' a=',num2str(a),' b=',num2str(b),'_error'],'-djpeg')
hold off

end

