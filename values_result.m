
%{
x1  - initial values of prey
x2  - initial values of predator
X_2 - predator target
B   - prey restriction
alpha1,alpha2,beta1,beta2,w1 - coefficient
%}
x1 = 6; 
x2 = 3; 
X_2 = 7;
B=9;    

f = fopen('values.txt','wt');
n0 = 1;
nn = 50;
h = 0.01;
interval = n0:h:nn;

N = length(interval);

X1 = zeros(1,N);
X2 = zeros(1,N);
X1(1)=x1;
X2(1)=x2;

for alpha1 = 0.1:0.1:1.5
    
    for alpha2 = 0.1:0.1:1.5
        
        for beta1 = 0.1:0.1:1.5
            
            for beta2 = 0.1:0.1:1.5
                
                for w1 = 0.1:0.1:1 %for w1 = 0:1:10 
                 
                    for n=2:N
                    psi=X1(n-1)+ B*tanh(X2(n-1) - X_2);
                    P=(4*exp(2*psi))/((exp(2*psi)+1)^2);
                    u = -w1^(-1)*psi -(1+B*P)*(alpha1*X1(n-1) - beta1*X1(n-1)*X2(n-1));%+B*P*(X_2); 
    
                    X1(n)=X1(n-1) + h*(alpha1*X1(n-1) - beta1*X1(n-1)*X2(n-1)+u);
                    X2(n)=X2(n-1) + h*(-alpha2*X2(n-1)+beta2*X1(n-1)*X2(n-1));
                    psi_t=X2(n-1)-X_2;
                    
                    proc = N-N/100;%99 percent of the total,number
                    
                    
                    
                    if(n > proc && (abs(psi_t)/X_2)* 100 > 10)%checking deviation at the boundary values 
                    elseif(n==N)% until the counter reaches the end of the interval 
                    fprintf(f, '%4d,%0d,%0d,%0d,%0.2f,%0.2f,%0.2f,%0.2f,%0.2f\n',x1,x2,X_2,B,alpha1,alpha2,beta1,beta2,w1);

                    end
                    
                    end
                    
                end
                
            end
            
        end
           
    end

end

fclose(f);

