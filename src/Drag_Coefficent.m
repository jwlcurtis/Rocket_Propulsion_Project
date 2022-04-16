function [cd1]=Drag_Coefficent(mach)
    if mach<0.6 
        cd1=0.15;
    end
    if mach>=0.6 && mach<1.2
        cd1=-0.12+0.45*mach;
    end
    if mach>=1.2 && mach<=1.8
        cd1=0.76-0.283*mach;
    end
    if mach>1.8 && mach<=4.0
        cd1=0.311-0.034*mach;
    end
    if mach>4.0
        cd1=0.175;
    end
end