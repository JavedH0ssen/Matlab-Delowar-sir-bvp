x=[331,334,337,340,343]
y=[618.9182,638.9182,659.9182,674.9182,694.9182]
% loglog(x,y);
plot(log10(x),log10(y),'-ok')
xlim([log10(x(1)) log10(x(length(x)))])
ylim([log10(y(1)) log10(y(length(y)))])
xlabel('log(T)');
ylabel('log(Es)')