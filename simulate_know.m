function [beta_new] =simulate_know(X, beta, lambda, rho, M,S)
[n m]=size(X);
MAX_ITER =50;
f=zeros(size(X,2),1);
Z = zeros(size(X,2),1);
R = zeros(size(X,2),1);
judge=0;
I=eye(m);
for i=1:size(S,2)
      I(S(i),S(i))=0;
end
for i = 1 : n
    judge=M'*X(i,:)'+M*X(i,:)'+judge;
end
        judge=judge/n;
        history(1) = 0;
    for i = 1 : n
        history(1) = (X(i,:)'-f)'*M*(X(i,:)'-f)+history(1);
    end
        history(1) = history(1) / n;            
for k = 2 : MAX_ITER
 %% update f
    f = inv(2 * M + rho*( eye(size(X,2))) ) * ( judge + rho * ( Z-R) );
 %% Update Z
  C = f+R;
  par=lambda/ rho;
 for i=1:size(X,2)
     if abs(C(i))<= par*I(i,i)
         Z(i)=0;
     else
         Z(i)=sign(C(i))*(abs(C(i)) - par*I(i,i));
     end
 end
    %%update U
    Z(S)=0;
     beta_new=Z;
    R=R+f-Z;
    history(k) = 0;
    for i = 1 : n
        history(k) = (X(i,:)'-beta_new)'*M*(X(i,:)'-beta_new)+history(k);
    end
        history(k) = history(k) / n;
end
%plot(1:k, history, 'k', 'MarkerSize', 10, 'LineWidth', 2);
 end





