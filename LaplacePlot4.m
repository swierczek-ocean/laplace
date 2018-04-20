function error = LaplacePlot4(True,Weeks,T,char,a,b,sw)
colors

Q = sort(True);
Z = sort(Weeks);
ps = 14;
ss = 5.5;
Swkc = 12;

if(Q(end)-Q(1)>1e5)&&(Q(1)>0)&&(Z(1)>0)
    
    figure
    h1 = semilogy(T,True,'LineWidth',3,'Color',Color(:,35));
    hold on
    h2 = semilogy(T,Weeks,'.','MarkerSize',ps,'Color',Color(:,Swkc));
    semilogy(T,Weeks,'o','MarkerSize',ss,'Color',Color(:,Swkc));
    title(['L',num2str(char),' Shift Weeks method ',num2str(sw),' a=',num2str(a),' b=',num2str(b)])
    xlabel('time')
    ylabel('f(t)')
    legend([h1(1),h2(1)],'True f(t)','Shift Weeks','Location','northwest')
    print(['L',num2str(char),' sw=',num2str(sw),' a=',num2str(a),' b=',num2str(b),'_approx_v_true_SWeeks'],'-djpeg')
    hold off
    
else
    
    figure
    h1 = plot(T,True,'LineWidth',3,'Color',Color(:,35));
    hold on
    h2 = plot(T,Weeks,'.','MarkerSize',ps,'Color',Color(:,Swkc));
    plot(T,Weeks,'o','MarkerSize',ss,'Color',Color(:,Swkc));
    title(['L',num2str(char),' Shift Weeks method ',num2str(sw),' a=',num2str(a),' b=',num2str(b)])
    xlabel('time')
    ylabel('f(t)')
    legend([h1(1),h2(1)],'True f(t)','Shift Weeks','Location','northwest')
    print(['L',num2str(char),' sw=',num2str(sw),' a=',num2str(a),' b=',num2str(b),'_approx_v_true_SWeeks'],'-djpeg')
    hold off
    
end

ind = find(abs(True)==0);
BromRelError = abs(Weeks-True)./abs(True);
BromRelError(ind) = abs(Weeks(ind)-True(ind));
error = BromRelError;

figure
semilogy(T,real(BromRelError),'.','MarkerSize',ps,'Color',Color(:,Swkc))
hold on
semilogy(T,real(BromRelError),'o','MarkerSize',ss,'Color',Color(:,Swkc))
title(['L ',num2str(char),' Shift Weeks method ',num2str(sw),' a=',num2str(a),' b=',num2str(b),' error'])
xlabel('time')
ylabel('Relative error')
print(['L',num2str(char),' sw=',num2str(sw),' a=',num2str(a),' b=',num2str(b),'_error_SWeeks'],'-djpeg')
hold off

end

