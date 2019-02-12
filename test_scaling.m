%% test the exact scaling
x = 0.03:0.01:0.6;
x = x'+2;
xi = 1./log(x/2);
ft = fittype('a*x+b');
epsilon = (x-2)/2;
myfit = fit(log(epsilon),log(xi),ft,'StartPoint',[-1 1])
loglog(epsilon,xi,'d')
xlim([0.015 0.3])