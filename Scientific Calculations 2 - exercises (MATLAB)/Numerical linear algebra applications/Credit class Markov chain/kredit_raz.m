function [V,D,uk] = kredit_raz (n)

    m=8;
    load('model_kredit_A.mat');

    u0=zeros(m,1);
    u0(1,1)=1;
    uk=u0;
    k=0;
    
    figure()
    axis([0.5 8.5 0 1])
    xlabel('Stanje');
    ylabel('Gustoæa');
    title(sprintf('Iteracija broj %d',k));
    set(gca,'XTick',1:8);
    set(gca,'XTickLabel',{'AAA','AA','A','BBB','BB','B','CCC','D'});
    grid on;
    pause(0.2);
    set(gca,'YTick',0:0.2:1);
    x=1:8;
    
    for i=1:n
        plot(x,uk,'r-');
        uk=A*uk;
        k=k+1;
        axis([0.5 8.5 0 1])
        xlabel('Stanje');
        ylabel('Gustoæa');
        title(sprintf('Iteracija broj %d',k));
        set(gca,'XTick',1:8);
        set(gca,'XTickLabel',{'AAA','AA','A','BBB','BB','B','CCC','D'});
        grid on;
        pause(0.1);
    end
    
    [V,D]=eig(A);

end