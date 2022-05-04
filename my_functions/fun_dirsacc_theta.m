function th = fun_dirsacc_theta(xini,yini,xend,yend)
    th = atan((yend-yini)./(xend-xini));
    for i=1:length(th)
        if      ((yend(i)-yini(i))>0 && (xend(i)-xini(i))<0); th(i)=th(i)+pi;
        elseif  ((yend(i)-yini(i))<0 && (xend(i)-xini(i))<0); th(i)=th(i)+pi;
        elseif  ((yend(i)-yini(i))<0 && (xend(i)-xini(i))>0); th(i)=th(i)+2*pi;
        end
    end
end