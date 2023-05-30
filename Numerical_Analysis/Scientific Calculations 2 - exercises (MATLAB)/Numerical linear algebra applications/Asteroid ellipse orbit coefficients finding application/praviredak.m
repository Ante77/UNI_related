function [A] = praviredak (A,p,i)
    
    n=max(size(A));
    A(i,1)=p(1)*p(2);
    A(i,2)=p(2)*p(2);
    A(i,3)=p(1);
    A(i,4)=p(2);
    A(i,5)=1;
    
end