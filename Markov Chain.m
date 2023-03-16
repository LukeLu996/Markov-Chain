clear;
%function [chain,state]=markov(T,n,s0,V);
T=[0.45,0.48,0.07;0.05,0.70,0.25;0.01,0.50,0.49]; %  T is transition matrix
[r c]=size(T);   % r is # of rows, c is # of columns of T
n=100000; %  n is number of periods to simulate
V=[1:r]; %  V is the quantity corresponding to each state 
s0=1; %  s0 is initial state (initial probabilities)
X=rand(n-1,1);  % generate (n-1) random numbers drawn from uniform distribution on [0,1], each number to be used in one simulation. 
s=zeros(r,1);   % initiate the state vector "s" to be a rx1 zero vector
s(s0)=1;        % change the "s0"th element of "s" to 1
cum=T*triu(ones(size(T))); % "triu(ones(size(T)))" generates an upper triangular matrix with all elements equal to 1 
% cum is a rxr matrix whose ith column is the cumulative sum from the 1st column to the ith column  
% the ith row of cum is the cumulative distribution for the next period given the current state. 
for k=1:length(X);  % "length(X)" returns the size of the longest dimension of X. "k" indicates the kth simulation.
    state(:,k)=s;     % state is a matrix recording the number of the realized state at time k
  ppi=[0 s'*cum];   % this is the conditional cumulative distribution for the next period given the current state s
   s=((X(k)<=ppi(2:r+1)).*(X(k)>ppi(1:r)))'; 
   % compares each element of ppi(2:r+1) or ppi(1:r) with a scalar X(k), and   
   % returns 1 if the inequality holds and 0 otherwise    
   % this formula assigns 1 when both inequalities hold, and 0 otherwise
end;
   chain=V*state
%for the first column, a1,a2 and a3 are the original numbers of the first
%column of T^L1,but in order to let the first column to be stable within a
%reasonable range, 4 significant digits are kept to get A1,A2 and A3.
L1=1;
a1=T(1,1);
a2=T(2,1);
a3=T(3,1);
A1=roundn(a1,-4);
A2=roundn(a2,-4);
A3=roundn(a3,-4);
while A1~=A2 | A2~=A3 | A1~=A3
    L1=L1+1;
    T1=T^L1;
    a1=T1(1,1);
    a2=T1(2,1);
    a3=T1(3,1);
    A1=roundn(a1,-4);
    A2=roundn(a2,-4);
    A3=roundn(a3,-4);
end
L1=L1 %the number of the power of martix T to make the first column stable

%similar as the first column
L2=1;
b1=T(1,2);
b2=T(2,2);
b3=T(3,2);
B1=roundn(b1,-4);
B2=roundn(b2,-4);
B3=roundn(b3,-4);
while B1~=B2 | B2~=B3 | B1~=B3
    L2=L2+1;
    T2=T^L2;
    b1=T2(1,2);
    b2=T2(2,2);
    b3=T2(3,2);
    B1=roundn(b1,-4);
    B2=roundn(b2,-4);
    B3=roundn(b3,-4);
end
L2=L2 %the number of the power of martix T to make the second column stable
 
L3=1;
c1=T(1,3);
c2=T(2,3);
c3=T(3,3);
C1=roundn(c1,-4);
C2=roundn(c2,-4);
C3=roundn(c3,-4);
while C1~=C2 | C2~=C3 | C1~=C3
    L3=L3+1;
    T3=T^L3;
    c1=T3(1,3);
    c2=T3(2,3);
    c3=T3(3,3);
    C1=roundn(c1,-4);
    C2=roundn(c2,-4);
    C3=roundn(c3,-4);
end
L3=L3 %the number of the power of martix T to make the third column stable
%in order to find the long-run probabilities, the biggest number from
%L1,L2 and L3 are needed as L
if L1>L2 & L1>L3
    L=L1
end
if L2>L1 & L2>L3
    L=L2
end
if L3>L1 & L3>L2
    L=L3
end
T=T^L
P=T(1,:)'