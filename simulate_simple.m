function [beta_new] =simulate_simple(X, beta,  M)
[n m]=size(X);
MAX_ITER =50;
f=zeros(size(X,2),1);
judge=0;
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
    f =( pinv( 2 * M ) ) * ( judge );
     beta_new=f;  
     history(k) = 0;
    for i = 1 : n
        history(k) = (X(i,:)'-beta_new)'*M*(X(i,:)'-beta_new)+history(k);
    end
        history(k) = history(k) / n;

end
%plot(1:k, history, 'k', 'MarkerSize', 10, 'LineWidth', 2);
end



