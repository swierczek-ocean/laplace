function LaplacePlot(True,NAB,T,char)
colors

Q = sort(True);

if(Q(end)-Q(1)>1e5)
    
    figure
    semilogy(T,True,'LineWidth',3,'Color','black')
    hold on
    semilogy(T,NAB,'+','MarkerSize',8,'Color',Color(:,9))
    title(['L',num2str(char)])
    xlabel('time')
    ylabel('f(t)')
    legend('True f(t)','Bromwich adapt')
    print(['L',num2str(char)],'-djpeg')
    hold off
    
else
    
    figure
    plot(T,True,'LineWidth',3,'Color','black')
    hold on
    plot(T,NAB,'+','MarkerSize',8,'Color',Color(:,9))
    title(['L',num2str(char)])
    xlabel('time')
    ylabel('f(t)')
    legend('True f(t)','Bromwich adapt')
    print(['L',num2str(char)],'-djpeg')
    hold off
    
end

BromRelError = abs(NAB-True)./abs(True);

figure
semilogy(T,real(BromRelError),'+','MarkerSize',8,'Color',Color(:,9))
title(['L',num2str(char),' error'])
xlabel('time')
ylabel('Relative error')
print('L1error','-djpeg')
print(['L',num2str(char),'error'],'-djpeg')
hold off

end

