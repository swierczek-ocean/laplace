function error = LaplacePlot(True,NAB,T,char,ub)
colors

Q = sort(True);

if(Q(end)-Q(1)>1e5)&&(Q(1)>0)
    
    figure
    h1 = semilogy(T,True,'LineWidth',2.5,'Color',Color(:,11));
    hold on
    h2 = semilogy(T,NAB,'*','MarkerSize',7,'Color',Color(:,12));
    title(['L',num2str(char),' loc=',num2str(2*ub)])
    xlabel('time')
    ylabel('f(t)')
    legend([h1(1),h2(1)],'True f(t)','Bromwich adapt')
    print(['L',num2str(char),' loc=',num2str(2*ub)],'-djpeg')
    hold off
    
else
    
    figure
    h1 = plot(T,True,'LineWidth',2.5,'Color',Color(:,11));
    hold on
    h2 = plot(T,NAB,'*','MarkerSize',7,'Color',Color(:,12));
    title(['L',num2str(char),' loc=',num2str(2*ub)])
    xlabel('time')
    ylabel('f(t)')
    legend([h1(1),h2(1)],'True f(t)','Bromwich adapt')
    print(['L',num2str(char),' loc=',num2str(2*ub)],'-djpeg')
    hold off
    
end

ind = find(True==0);
BromRelError = abs(NAB-True)./abs(True);
BromRelError(ind) = abs(NAB(ind)-True(ind));
error = mean(BromRelError);

figure
semilogy(T,real(BromRelError),'*','MarkerSize',7,'Color',Color(:,12))
title(['L',num2str(char),' loc=',num2str(2*ub),' error'])
xlabel('time')
ylabel('Relative error')
print(['L',num2str(char),' error loc=',num2str(2*ub)],'-djpeg')
hold off

end

