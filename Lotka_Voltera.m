function Lotka_Voltera(x1,x2,X_2,B,alpha1,alpha2,beta1,beta2,w1)
%{
x1  - initial values of prey
x2  - initial values of predator
X_2 - predator target
B   - prey restriction
alpha1,alpha2,beta1,beta2,w1 - coefficient
%}
n0 = 1;%left border
nn = 50;%right border
h = 0.01;%step
interval = n0:h:nn;

N = length(interval);

X1 = zeros(1,N);
X2 = zeros(1,N);
X1(1)=x1;
X2(1)=x2;


for n=2:N
psi=X1(n-1)+ B*tanh(X2(n-1) - X_2);
P=(4*exp(2*psi))/((exp(2*psi)+1)^2);
u = -w1^(-1)*psi -(1+B*P)*(alpha1*X1(n-1) - beta1*X1(n-1)*X2(n-1));%+B*P*(X_2); 
    
X1(n)=X1(n-1) + h*(alpha1*X1(n-1) - beta1*X1(n-1)*X2(n-1)+u);
X2(n)=X2(n-1) + h*(-alpha2*X2(n-1)+beta2*X1(n-1)*X2(n-1));


end

Xvals = 1:N;

hold on;
plot(Xvals,X1,'b','LineWidth',1);
plot(Xvals,X2,'--','Color',[.1 .7 .7],'LineWidth',1);
set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',14,'DefaultAxesFontName','Times New Roman');
legend('Prey', 'Predator');
title('Predator - Prey');
grid on;
ylabel('X1,X2');
xlabel('time');
hold off;

end