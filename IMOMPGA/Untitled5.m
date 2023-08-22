aa=0;
bb=0;
cc = 0;
for i =1:100
    GA = a;
    MPGA = b;
    IMPGA = c;
    
    if GA == 44
        aa= aa+1
    end
    if MPGA == 44
        bb=bb+1
    end
    if IMPGA == 44
        cc=cc+1
    end
end
fprintf('GA:%d, MPGA:%d, IMPGA:%d\n',aa,bb,cc);