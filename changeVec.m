function []= changeVec(PlotHandle,t,yo,yf,color)
 
    S=size(t);
    for i=1:S,
        if(t(i)>0)
           hold on; plot([t(i) t(i)],[yo yf],color);
       end;
    end;
    