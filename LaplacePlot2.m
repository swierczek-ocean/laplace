function error = LaplacePlot2(True,Weeks,T,char)
colors

Q = sort(True);
Z = sort(Weeks);

if(Q(end)-Q(1)>1e5)&&(Q(1)>0)&&(Z(1)>0)
    
    figure
    h1 = semilogy(T,True,'LineWidth',2.75,'Color',Color(:,35));
    hold on
    h2 = semilogy(T,Weeks,'.','MarkerSize',14,'Color',Color(:,18));
    semilogy(T,Weeks,'*','MarkerSize',7,'Color',Color(:,18));
    title(['Laplace function ',num2str(char)])
    xlabel('time')
    ylabel('f(t)')
    legend([h1(1),h2(1)],'True f(t)','Weeks')
    print(['L',num2str(char),'_approx_v_true_Weeks'],'-djpeg')
    hold off
    
else
    
    figure
    h1 = plot(T,True,'LineWidth',2.75,'Color',Color(:,35));
    hold on
    h2 = plot(T,Weeks,'.','MarkerSize',14,'Color',Color(:,18));
    plot(T,Weeks,'*','MarkerSize',7,'Color',Color(:,18));
    title(['Laplace function ',num2str(char)])
    xlabel('time')
    ylabel('f(t)')
    legend([h1(1),h2(1)],'True f(t)','Weeks')
    print(['L',num2str(char),'_approx_v_true_Weeks'],'-djpeg')
    hold off
    
end

ind = find(abs(True)<10e-03);
BromRelError = abs(Weeks-True)./abs(True);
BromRelError(ind) = abs(Weeks(ind)-True(ind));
error = BromRelError;

figure
semilogy(T,real(BromRelError),'.','MarkerSize',12,'Color',Color(:,18))
hold on
semilogy(T,real(BromRelError),'*','MarkerSize',6,'Color',Color(:,18))
title(['Laplace function ',num2str(char),' error'])
xlabel('time')
ylabel('Relative error')
print(['L',num2str(char),'_error_Weeks'],'-djpeg')
hold off

end

