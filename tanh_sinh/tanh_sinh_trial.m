clear;
clc;
close all;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% 
% Source:
% D. H. Bailey, X. S. Li and K. Jeyabalan, “A Comparison of Three
% High-Precision Quadrature Schemes,” Experimental Mathematics, 
% vol. 14 (2005), no.3, pg. 317–329
%
% 
% Trial functions:
% 
% * a: \int_{0}^{1} t*log(1 + t) dt = 1/4
% * b: \int_{0}^{pi/2} exp(t).*cos(t) dt = (exp(\pi/2)-1)/2
% * c: \int_{0}^{pi/2} log(cos(t)) dt = -\pi*log(2)/2
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% 
% Case a:
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

tols=(1:1:10);
integr=@(t) t.*log(1 + t);
ares=1/4;
res=zeros(length(tols),1);
a=0;
b=1;
for ij=1:length(tols)
    res(ij)=tanh_sinh_quad(integr,a,b,tols(ij));
end

figure
semilogy(tols,abs(res-ares),'LineWidth',3)
grid minor
xlabel('Tolerance d (10^{-d})')
ylabel('Absolute difference from the analytical solution')
dim = [.6 .7 .3 .1];
str = 'Integral: \int_{0}^{1} tlog(t)dt=0.25';
annotation('textbox',dim,'String',str,'FitBoxToText','on','BackgroundColor','w','FontSize',12);
title('Performance of tanh-sinh quadrature')

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% 
% Case b:
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

tols=(1:1:10);
integr=@(t) exp(t).*cos(t);
ares=(exp(pi/2)-1)/2;
res=zeros(length(tols),1);
a=0;
b=pi/2;
for ij=1:length(tols)
    res(ij)=tanh_sinh_quad(integr,a,b,tols(ij));
end

figure
semilogy(tols,abs(res-ares),'LineWidth',3)
grid minor
xlabel('Tolerance d (10^{-d})')
ylabel('Absolute difference from the analytical solution')
dim = [.5 .7 .3 .1];
str = 'Integral: \int_{0}^{\pi/2} log(cos(t)) dt =1.9052';
annotation('textbox',dim,'String',str,'FitBoxToText','on','BackgroundColor','w','FontSize',12);
title('Performance of tanh-sinh quadrature')

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% 
% Case c:
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

tols=(1:1:10);
integr=@(t) log(cos(t));
ares=-pi*log(2)/2;
res=zeros(length(tols),1);
a=0;
b=pi/2;
for ij=1:length(tols)
    res(ij)=tanh_sinh_quad(integr,a,b,tols(ij));
end

figure
semilogy(tols,abs(res-ares),'LineWidth',3)
grid minor
xlabel('Tolerance d (10^{-d})')
ylabel('Absolute difference from the analytical solution')
dim = [.5 .7 .3 .1];
str = 'Integral: \int_{0}^{\pi/2}log(cos(t))dt =-1.0888';
annotation('textbox',dim,'String',str,'FitBoxToText','on','BackgroundColor','w','FontSize',12);
title('Performance of tanh-sinh quadrature')