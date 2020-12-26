function [beta_new] =simulate_SGL(X, beta, lambda, rho, M,w,Group)
[n m]=size(X);
MAX_ITER =50;
I=zeros(15*4,15);
for i=1:4
   row=find(Group(i,:)~=0);
   [group_n,group_m]=size(row);
   for j=1:group_m
   I(15*(i-1)+row(j),row(j))=1;
   end
end
I_l2sum=zeros(15,15);
for i=1:4
    I_l2sum= I((15*(i-1)+1):15*i,:)'*I((15*(i-1)+1):15*i,:)+ I_l2sum;
end
V = zeros(size(X,2),4);
U = zeros(size(X,2),4);
f = zeros(size(X,2),1);
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
 suml2_V=zeros(size(X,2),1);
 for i=1:4
     suml2_V= I((15*(i-1)+1):15*i,:)'*V(:,i)+suml2_V;
 end
  suml2_U=zeros(size(X,2),1);
 for i=1:4
     suml2_U= I((15*(i-1)+1):15*i,:)'*U(:,i)+suml2_U;
 end
    f = inv(2 * M + rho*( I_l2sum ) )  * ( judge + rho * ( suml2_V-suml2_U) );
    
 %% update V
 for j=1:4
     T=zeros(4,1);
     T=I((15*(j-1)+1):15*j,:)*f+U(:,j);
     par=w(j)/rho;
if norm(T,2)<=par
    V(:,j) = 0;
 else 
   V(:,j) =( 1-par/norm(T,2) ) *T;
 end
   end
  %%update U
     for j=1:4
  U(:,j)=U(:,j)+I((15*(j-1)+1):15*j,:)*f-V(:,j);
     end
  beta_new=sum(V,2);
    history(k) = 0;
    for i = 1 : n
        history(k) = (X(i,:)'-beta_new)'*M*(X(i,:)'-beta_new)+history(k);
    end
        history(k) = history(k) / n;
end
% plot(1:k, history, 'k', 'MarkerSize', 10, 'LineWidth', 2);
 end
 




