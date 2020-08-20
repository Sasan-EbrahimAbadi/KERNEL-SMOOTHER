%% Kernel for Pinf_finite size_random percolation
clear
Y=load('MCS32_giant.mat');
X=load('P32.mat');
sigma=0.01;
three_sigma=3*sigma;
DX=0.001;
Bin=0.1:DX:1;
KERNEL=zeros(1,length(Bin));
 
for B=1:length(Bin)
    MIN=max(0, Bin(B)-three_sigma);
    MAX=min(1, Bin(B)+three_sigma); 
    IND1=find(X >= MIN & X < MAX);
    numerator=0; denominator=0;
    for B1=1:length(IND1)
        numerator=numerator+Y(IND1(B1))*(exp(-(((Bin(B)-X(IND1(B1)))^2)/(2*(sigma^2)))));
        denominator=denominator+exp(-(((Bin(B)-X(IND1(B1)))^2)/(2*(sigma^2))));
        KERNEL(1,B)=numerator/denominator;
        
    end
end
figure
plot(Bin,KERNEL,'o')
grid on
%save('KERNEL.mat','KERNEL')

