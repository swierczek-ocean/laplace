function error = LaplacePlot(True,NAB,T,char,ub)
colors

Q = sort(True);

if(Q(end)-Q(1)>1e5)
    
    figure
    semilogy(T,True,'LineWidth',2.5,'Color',Color(:,11))
    hold on
    semilogy(T,NAB,'*','MarkerSize',7,'Color',Color(:,12))
    title(['L',num2str(char),' cl=',num2str(2*ub)])
    xlabel('time')
    ylabel('f(t)')
    legend('True f(t)','Bromwich adapt')
    print(['L',num2str(char),' cl=',num2str(2*ub)],'-djpeg')
    hold off
    
else
    
    figure
    plot(T,True,'LineWidth',2.5,'Color',Color(:,11))
    hold on
    plot(T,NAB,'*','MarkerSize',7,'Color',Color(:,12))
    title(['L',num2str(char),' cl=',num2str(2*ub)])
    xlabel('time')
    ylabel('f(t)')
    legend('True f(t)','Bromwich adapt')
    print(['L',num2str(char),' cl=',num2str(2*ub)],'-djpeg')
    hold off
    
end

BromRelError = abs(NAB-True)./abs(True);
error = mean(BromRelError);

figure
semilogy(T,real(BromRelError),'*','MarkerSize',7,'Color',Color(:,12))
title(['L',num2str(char),' cl=',num2str(2*ub),' error'])
xlabel('time')
ylabel('Relative error')
print(['L',num2str(char),' cl=',num2str(2*ub),' error'],'-djpeg')
hold off

end

