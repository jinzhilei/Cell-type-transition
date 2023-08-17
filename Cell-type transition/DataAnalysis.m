

% Number and proportion of cells
clear; clc;
outpath='output/sys/';

N=30;
M=5;

G=zeros(14,M);

for l=1:M
    
    E=zeros(N,4);
    
    F=zeros(N,3);

    for k=1:N
        D=load(strcat(outpath,num2str(k),'_',num2str(l),'.sys'));
        E(k,1)=mean(D(8000:end,2));
        E(k,2)=mean(D(8000:end,3));
        E(k,3)=mean(D(8000:end,4));
        E(k,4)=mean(D(8000:end,5));
      
   
        F(k,1)=100*E(k,2)*1.0/E(k,1);
        F(k,2)=100*E(k,3)*1.0/E(k,1);
        F(k,3)=100*E(k,4)*1.0/E(k,1);
    end
    G(1:4,l)=mean(E); %Mean cell number
    G(5:8,l)=std(E); %Cell number standard deviation
    
    G(9:11,l)=mean(F); %The average ratio of the number of cells
    G(12:14,l)=std(F); %The standard deviation of the ratio of cell numbers
    
end

save('output/figdat/sys.dat','G','-ascii');

% Transition probability

clear; clc;

outpath='output/tran/';

N=30;
M=5;

G=zeros(18,M);

for l=1:M
    
    B=zeros(3,3);
    E=zeros(N,9);
    a=[0 1 3];
    for k=1:N
        A=load(strcat(outpath,num2str(k),'_',num2str(l),'.tran'));
        for i=1:3
            J=find(A(:,2)==a(i));
            Atemp=A(J,1);
            for j=1:3
                Jtemp=find(Atemp==a(j));
                B(i,j)=size(Jtemp,1);
            end
        end
        B1=B';
        C=diag(100./sum(B1))*B;
        D=reshape(C',[9,1]);
        E(k,:)=D;    
    end
    G(1:9,l)=mean(E);   
    G(10:18,l)=std(E);  
    
    %1-9th rows represent the mean transition probability
    %10-18th rows represents the standard deviation of the transition probability

    % 1-th row 1:SC-SC; 5-th row : TA1-TA1;  9-th row: TA2-TA2;
    % 2-th row :SC-TA1; 3-th row: SC-TA2;
    % 4-th row  :TA1-SC;  7-th row: TA2-SC;
    % 6-th row :TA1-TA2;  8-th row: TA2-TA1;
end

save('output/figdat/tran.dat','G','-ascii');

clear; clc;
