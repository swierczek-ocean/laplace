function LaplacePlot(True,NAB,T)
colors

figure
plot(T,True,'LineWidth',3,'Color','black')
hold on
plot(T,NAB,'+','MarkerSize',8,'Color',Color(:,9))
title('L 1')
xlabel('time')
ylabel('f(t)')
legend('True f(t)','Bromwich adapt')
print('L1','-djpeg')
hold off

BromRelError = abs(NAB-True)./abs(True);

figure
semilogy(T,real(BromRelError),'+','MarkerSize',8,'Color',Color(:,9))
title('L 1 Error')
xlabel('time')
ylabel('Relative error')
print('L1error','-djpeg')
hold off

end

