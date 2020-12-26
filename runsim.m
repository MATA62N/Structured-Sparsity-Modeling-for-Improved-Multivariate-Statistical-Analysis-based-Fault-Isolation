load('Rawdata.mat');%The faults in Rawdata are set in variables 2, 3, 15
%% PCA-based process monitoring
    X_train = X_data(1:700,:);     %Training set
    X_test = X_data(701:end,:);      %Test set
    X_mean = mean(X_train);
    X_std = std(X_train);
    [X_row,X_col] = size(X_train);
    X_train=(X_train-repmat(X_mean,X_row,1))./repmat(X_std,X_row,1);        % Normalization
    X_train(isnan(X_train)==1)=0;       %  Replace NAN with 0
    [coeff, score,latent] = pca(X_train);
    var_exp = latent*100/(sum(latent));
    num_pc = find(cumsum(var_exp)>85,1);        % Retain more than 85% variance
    P = coeff;
    T2UCL(:,1)=num_pc*(X_row-1)*(X_row+1)/(X_row*(X_row-num_pc))*finv(0.99, num_pc, X_row-num_pc);     %Significance level is 0.01
    for i = 1:3
        theta(i) = sum((latent(num_pc+1:X_col)).^i);
    end
    h0 = 1 - 2*theta(1)*theta(3)/(3*theta(2)^2);
    Ca = norminv(0.99,0,1);
    QUCL(:,1) = theta(1)*(h0*Ca*sqrt(2*theta(2))/theta(1) + 1 + theta(2)*h0*(h0 - 1)/theta(1)^2)^(1/h0);
    % Online Monitoring
    [n,m] = size(X_test);
    X_test=(X_test-repmat(X_mean,n,1))./repmat(X_std,n,1);
    X_test(isnan(X_test)==1)=0;
    % Calculate T2 statistics, SPE statistics
    [P_r,P_c] = size(P(:,1:num_pc)*P(:,1:num_pc)');
    I = eye(P_r,P_c);
    for i = 1:n
        T2(i,1) = X_test(i,:)*P(:,1:num_pc)*pinv(diag(latent(1:num_pc,:)))*P(:,1:num_pc)'*X_test(i,:)';
        Q(i,1) = X_test(i,:)*(I - P(:,1:num_pc)*P(:,1:num_pc)')*X_test(i,:)';
    end
    % T2 and SPE statistics plotting
    figure('Position',[300 200 810 560])
    subplot(2,1,1);
    plot(1:n,T2(:,1),'k'); 
    box off
    ylabel('T^2','FontSize',20,'Fontname','Times New Roman');
    hold on;
    line([0,n],[T2UCL(:,1),T2UCL(:,1)],'LineStyle','--','Color','r');
    set(gca,'FontSize',26,'Fontname', 'Times New Roman');  
    subplot(2,1,2);
    plot(1:n,Q(:,1),'k');
    box off
    xlabel('Sample number','FontSize',20,'Fontname','Times New Roman')
    ylabel('Q','FontSize',20,'Fontname','Times New Roman');
    hold on;
    line([0,n],[QUCL(:,1),QUCL(:,1)],'LineStyle','--','Color','r');
    set(gca,'FontSize',26,'Fontname', 'Times New Roman');
%%  fault isolation 
M= P(:,1:num_pc)*pinv(diag(latent(1:num_pc,:)))*P(:,1:num_pc)';
X= X_test;
[n,p] = size(X);
density = 100/p;
beta = sprandn(p,1,density);
Group=[];       %Set group information
Group(1,[4 9 13])=1;Group(2,[3 11 15])=1;Group(3,[1 2 6 7 10])=1;Group(4,[5 8 12 14])=1;
group_number = [3 3 5 4];
%Set the weight of the tree structure
%%Set the weight of the l2 norm in the tree structure
sv=zeros(23,1);
sv([1 2 3 4])=1/6;
sv(5)=2/6;
sv(6)=3/6;
sv(7)=4/6;
sv(8)=5/6;
sv(9)=0.9;
gv=1-sv;
w_tree=zeros(9,1);
w_tree(1)=gv(1)*sv(8)*sv(9);
w_tree(2)=gv(2)*sv(7)*sv(8)*sv(9);
w_tree(3)=gv(3)*sv(5)*sv(6)*sv(7)*sv(8)*sv(9);
w_tree(4)=gv(4)*sv(9);
w_tree(5)=gv(5)*sv(6)*sv(7)*sv(8)*sv(9);
w_tree(6)=gv(6)*sv(7)*sv(8)*sv(9);
w_tree(7)=gv(7)*sv(8)*sv(9);
w_tree(8)=gv(8)*sv(9);
w_tree(9)=gv(9);
%%Set the weight of the l1 norm in the tree structure
W_tree=zeros(6,1);
W_tree(1)=sv(1)*sv(8)*sv(9);
W_tree(2)=sv(2)*sv(7)*sv(8)*sv(9);
W_tree(3)=sv(3)*sv(5)*sv(6)*sv(7)*sv(8)*sv(9);
W_tree(4)=sv(4)*sv(9);
W_tree(5)=sv(5)*sv(6)*sv(7)*sv(8)*sv(9);
W_tree(6)=sv(6)*sv(7)*sv(8)*sv(9);
Group_tree_l2=zeros(9,15);
Group_tree_l2(1,[4 9 13])=1;Group_tree_l2(2,[3 11 15])=1;
Group_tree_l2(3,[6 7 10])=1;Group_tree_l2(4,[5 8 12 14])=1;
Group_tree_l2(5,[1 6 7 10])=1;Group_tree_l2(6,[1 2 6 7 10])=1;
Group_tree_l2(7,[1 2 6 7 10 3 11 15])=1;Group_tree_l2(8,[1 2 6 7 10 3 11 15 4 9 13])=1;Group_tree_l2(9,[1:15])=1;
Group_tree_l1=zeros(6,15);
Group_tree_l1(1,[4 9 13])=1;Group_tree_l1(2,[3 11 15])=1;Group_tree_l1(3,[6 7 10])=1;Group_tree_l1(4,[5 8 12 14])=1;
Group_tree_l1(5,1)=1;Group_tree_l1(6,2)=1;
%parameter settings
lambda_l1 =0.1;
lambda_know=0.1;
lambda_GL=0.3;
lambda_SGL=0.7;
alpha=0.95;
rho = 1.2;
lambda_tree=0.32;
%Set the weight of Group lasso
w_GL=zeros(4,1);
for i=1:4
    w_GL(i)=sqrt(group_number(i))*lambda_GL;
end
%Set the weight of Sparse Group lasso
w_SGL=zeros(4,1);
for i=1:4
    w_SGL(i)=sqrt(group_number(i))*lambda_SGL*(1-alpha);
end
%Set support set
S=[11];     %Assuming that variable 11 is normal
%calculate
beta_new_simple=simulate_simple(X, beta, M)
beta_new_l1=simulate_l1(X, beta, lambda_l1, rho, M)
beta_new_know =simulate_know(X, beta, lambda_know, rho, M,S)
beta_new_GL=simulate_GL(X, beta, lambda_GL, rho, M,w_GL,Group)
beta_new_SGL =simulate_SGL(X, beta, lambda_SGL, rho, M,w_SGL,Group,alpha)
beta_new_tree = simulate_tree(X, beta,lambda_tree,rho, M,w_tree,Group_tree_l2,Group_tree_l1,W_tree)
%plot
figure()
subplot(3,2,1)
bar(abs(beta_new_simple))
title('Conventional reconstruction based contribution plot','FontSize',14)
box off
set(gca,'FontSize',18,'Fontname', 'Times New Roman');
axis([0.5 15.5,-inf,inf])
xlabel('Variable Number');
ylabel('Contribution');
subplot(3,2,2)
bar(abs(beta_new_l1))
title('only the l_1 penalty considered','FontSize',14)
box off
set(gca,'FontSize',18,'Fontname', 'Times New Roman');
axis([0.5 15.5,-inf,inf])
xlabel('Variable Number');
ylabel('Contribution');
subplot(3,2,3)
bar(abs(beta_new_know))
title('partially known sparse support','FontSize',14)
box off
set(gca,'FontSize',18,'Fontname', 'Times New Roman');
axis([0.5 15.5,-inf,inf])
xlabel('Variable Number');
ylabel('Contribution');
subplot(3,2,4)
bar(abs(beta_new_GL))
title('Group Lasso','FontSize',14)
box off
set(gca,'FontSize',18,'Fontname', 'Times New Roman');
axis([0.5 15.5,-inf,inf])
xlabel('Variable Number');
ylabel('Contribution');
subplot(3,2,5)
bar(abs(beta_new_SGL ))
title('Sparse Group Lasso','FontSize',14)
box off
set(gca,'FontSize',18,'Fontname', 'Times New Roman');
axis([0.5 15.5,-inf,inf])
xlabel('Variable Number');
ylabel('Contribution');
subplot(3,2,6)
bar(abs(beta_new_tree))
title('Tree-structured sparsity','FontSize',14)
box off
set(gca,'FontSize',18,'Fontname', 'Times New Roman');
axis([0.5 15.5,-inf,inf])
xlabel('Variable Number');
ylabel('Contribution');