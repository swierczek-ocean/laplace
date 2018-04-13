function error = LaplacePlot(True,NAB,T,char,ub)
colors

Q = sort(True);
Z = sort(NAB);

if(Q(end)-Q(1)>1e5)&&(Q(1)>0)&&(Z(1)>0)
    
    figure
    h1 = semilogy(T,True,'LineWidth',2.75,'Color',Color(:,35));
    hold on
    h2 = semilogy(T,NAB,'.','MarkerSize',18,'Color',Color(:,18));
    semilogy(T,NAB,'*','MarkerSize',8,'Color',Color(:,18));
    title(['Laplace function ',num2str(char),' LOC=',num2str(2*ub)])
    xlabel('time')
    ylabel('f(t)')
    legend([h1(1),h2(1)],'True f(t)','Bromwich adapt')
    print(['L',num2str(char),' ',num2str(2*ub),'_approx_v_true'],'-djpeg')
    hold off
    
else
    
    figure
    h1 = plot(T,True,'LineWidth',2.75,'Color',Color(:,35));
    hold on
    h2 = plot(T,NAB,'.','MarkerSize',18,'Color',Color(:,18));
    plot(T,NAB,'*','MarkerSize',8,'Color',Color(:,18));
    title(['Laplace function ',num2str(char),' LOC=',num2str(2*ub)])
    xlabel('time')
    ylabel('f(t)')
    legend([h1(1),h2(1)],'True f(t)','Bromwich adapt')
    print(['L',num2str(char),' ',num2str(2*ub),'_approx_v_true'],'-djpeg')
    hold off
    
end

ind = find(abs(True)<10e-03);
BromRelError = abs(NAB-True)./abs(True);
BromRelError(ind) = abs(NAB(ind)-True(ind));
error = BromRelError;

figure
semilogy(T,real(BromRelError),'.','MarkerSize',18,'Color',Color(:,18))
hold on
semilogy(T,real(BromRelError),'*','MarkerSize',8,'Color',Color(:,18))
title(['Laplace function ',num2str(char),' error,',' LOC=',num2str(2*ub)])
xlabel('time')
ylabel('Relative error')
print(['L',num2str(char),' ',num2str(2*ub),'_error'],'-djpeg')
hold off

end

