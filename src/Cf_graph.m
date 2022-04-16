function []=Cf_graph(AeAt,cf,gamma,fig)
    e=1:0.01:100;
    cf2=zeros(1,length(e));
    cf4=zeros(1,length(e));
    cf10=zeros(1,length(e));
    cf25=zeros(1,length(e));
    cf68=zeros(1,length(e));
    cf150=zeros(1,length(e));
    cf500=zeros(1,length(e));
    cf1000=zeros(1,length(e));
    cfinf=zeros(1,length(e));
    cfmin=zeros(1,length(e));
    cfopt=zeros(1,length(e));
    for i=1:length(e)
        if i==1
            m0=1;
        else
        m0=2; % inital guess
        end
        cf2(i)=Cf_solver(e(i),1/2,m0,gamma);
        cf4(i)=Cf_solver(e(i),1/4,m0,gamma);
        cf10(i)=Cf_solver(e(i),1/10,m0,gamma);
        cf25(i)=Cf_solver(e(i),1/25,m0,gamma);
        cf68(i)=Cf_solver(e(i),1/68,m0,gamma);
        cf150(i)=Cf_solver(e(i),1/150,m0,gamma);
        cf500(i)=Cf_solver(e(i),1/500,m0,gamma);
        cf1000(i)=Cf_solver(e(i),1/1000,m0,gamma);
        cfinf(i)=Cf_solver(e(i),0,m0,gamma);
        cfmin(i)=-0.0445*log(e(i))^2+0.5324*(log(e(i)))+0.1843;
        cfopt(i)=Cf_opt_solver(e(i),0,m0,gamma);
    end
    figure(fig)
    semilogx(e,cf2,'LineWidth',1.5)
    hold on
    semilogx(e,cf4,'LineWidth',1.5)
    semilogx(e,cf10,'LineWidth',1.5)
    semilogx(e,cf25,'LineWidth',1.5)
    semilogx(e,cf68,'LineWidth',1.5)
    semilogx(e,cf150,'LineWidth',1.5)
    semilogx(e,cf500,'LineWidth',1.5)
    semilogx(e,cf1000,'LineWidth',1.5)
    semilogx(e,cfinf,'LineWidth',1.5)
    semilogx(e,cfmin,"--",'LineWidth',1.5)
    semilogx(e,cfopt,"--",'LineWidth',1.5)
    semilogx(AeAt,cf,'LineWidth',1.5)
    hold off
    ylabel("c_f")
    xlabel("A_e/A_t")
    ylim([0.6 2])
    xlim([1 100])
    legend({"P_c/P_a=2","P_c/P_a=4","P_c/P_a=10","P_c/P_a=25","P_c/P_a=68","P_c/P_a=150","P_c/P_a=500","P_c/P_a=1000","P_c/P_a=\infty","Line of Seperation","Optimal c_f","Design"},'Location','northwest','NumColumns',3)
end