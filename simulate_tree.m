function [beta_new] =simulate_tree(X, beta,lambda_tree,rho, M,w_tree,Group_tree_l2,Group_tree_l1,W_tree)
[n m]=size(X);
MAX_ITER =50;
I=zeros(15*9,15);
for i=1:9
   row=find(Group_tree_l2(i,:)~=0);
   [group_n,group_m]=size(row);
   for j=1:group_m
   I(15*(i-1)+row(j),row(j))=1;
   end
end
I_l2sum=zeros(15,15);
for i=1:9
    I_l2sum= I((15*(i-1)+1):15*i,:)'*I((15*(i-1)+1):15*i,:)+ I_l2sum;
end

Group_l1=Group_tree_l1(find(W_tree>0),:);
I_l1=zeros(15*sum(W_tree~=0),15);
   for i=1:sum(W_tree~=0)
   row=find(Group_l1(i,:)~=0);
   [group_l1_n,group_l1_m]=size(row);
   for j=1:group_l1_m
   I_l1(15*(i-1)+row(j),row(j))=1;
   end
   end 
 I_l1sum=zeros(15,15);
for i=1:sum(W_tree~=0)
    I_l1sum= I_l1((15*(i-1)+1):15*i,:)'*I_l1((15*(i-1)+1):15*i,:)+ I_l1sum;
end  
V = zeros(size(X,2),9);
U = zeros(size(X,2),9);
Z = zeros(size(X,2),1);
R = zeros(size(X,2),1);
T = zeros(size(X,2),1);
O = zeros(size(X,2),1);
f = zeros(size(X,2),1);
    judge=0;
for i = 1 : n
    judge=M'*X(i,:)'+M*X(i,:)'+judge;
end
    judge=judge/n;
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
 for i=1:9
     suml2_V= I((15*(i-1)+1):15*i,:)'*V(:,i)+suml2_V;
 end
  suml2_U=zeros(size(X,2),1);
 for i=1:9
     suml2_U= I((15*(i-1)+1):15*i,:)'*U(:,i)+suml2_U;
 end
    f = inv(2 * M + rho*( I_l2sum+2*eye( size(X,2) ) ) ) * ( judge + rho * ( suml2_V+ Z-suml2_U-R+T-O) );
    
 %% update V
 for j=1:9
     P=zeros(15,1);
     P=I((15*(j-1)+1):15*j,:)*f+U(:,j);
     par=lambda_tree*w_tree(j)/rho;
if norm(P,2)<=par
    V(:,j) = 0;
 else 
   V(:,j) =( 1-par/norm(P,2) ) *P;
 end
   end
 

  %% Update Z
 C = f+R;
 W_non=W_tree(W_tree~=0);
 for i=1:size(X,2)
     L(i)=0;
      for j=1:sum(W_tree~=0)
 L(i)=W_non(j)*I_l1((15*(j-1)+i),i)+L(i);
      end
 end
 for i=1:size(X,2)
     if abs(C(i))<=(lambda_tree*L(i)/ rho )
         Z(i)=0;
     else
         Z(i)=sign(C(i))*( abs(C(i)) - ( lambda_tree*L(i)/ rho ) );
     end
 end
     %% Update T
 Y = f+O;
  par=lambda_tree/ rho;
 for i=1:size(X,2)
     if abs(Y(i))<= par
         T(i)=0;
     else
         T(i)=sign(Y(i))*( abs(Y(i)) - par );
     end
 end 
    %%update U
     for j=1:9
  U(:,j)=U(:,j)+I((15*(j-1)+1):15*j,:)*f-V(:,j);
     end
     R=R+f-Z;
%      O=O+f-T;
   beta_new=T;
    history(k) = 0;
    for i = 1 : n
        history(k) = (X(i,:)'-beta_new)'*M*(X(i,:)'-beta_new)+history(k);
    end
        history(k) = history(k) / n;
end
%plot(1:k, history, 'k', 'MarkerSize', 10, 'LineWidth', 2);
 end





