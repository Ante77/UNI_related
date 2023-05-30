function DF = df_zad2705(x)
    
    DF(1,1)=2*x(1);
    DF(1,2)=2*x(2);
    DF(2,1)=exp(x(1)-1);
    DF(2,2)=3*(x(2)^2);
end