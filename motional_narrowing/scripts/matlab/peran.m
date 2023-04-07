clear l td ta id ia bd ba x E ioa iod t ib ks1 ks2 ks3 ks4 io bin

!rm .tot.log .tot.xc tot.log .tot.name 2>/dev/null
!ls *.log > .tot.name
!cat *.log > .tot.log 
!cat *.xc > .tot.xc



tot = load('.tot.log');
c2 = 0;
if exist('.tot.xc','file')
    xtot = load('.tot.xc');
    c2 = 1;
end


if (~rem(length(tot),9))&c2
    if tot(1,2)
        l=length(tot);
        td=tot(1:9:l,1);
        ta=tot(1:9:l,2);
        id=tot(2:9:l,1);
        ia=tot(3:9:l,1);
        bd=tot(4:9:l,1);
        ba=tot(5:9:l,1);
        x =tot(8:9:l,1);
        E =tot(9:9:l,1);
        ioa = ia - ia./ba;
        iod = id - id./bd;
        bin = mean(xtot(:,2)-xtot(:,1));
    else 
        l=length(tot);
        td=tot(1:9:l,1);
        ta=tot(1:9:l,2);
        id=tot(2:9:l,1);
        ia=tot(3:9:l,1);
        bd=tot(4:9:l,1);
        ba=tot(5:9:l,1);
        x =tot(8:9:l,1);
        E =tot(9:9:l,1);
        ioa = ia - ia./ba;
        iod = id - id./bd;
        bin = .005;%mean(xtot(:,2)-xtot(:,1));
    end    
elseif (~rem(length(tot),7))
    if tot(1,2)

        l=length(tot);
        td=tot(1:7:l,1);
        ta=tot(1:7:l,2);
        id=tot(2:7:l,1);
        ia=tot(3:7:l,1);
        bd=tot(4:7:l,1);
        ba=tot(5:7:l,1);
        x =tot(6:7:l,1);
        E =tot(7:7:l,1);
        ioa = ia - ia./ba;
        iod = id - id./bd;
        bin = mean(xtot(:,2)-xtot(:,1));
    else 
        if isnan(tot(5,1))
            n = 6;
        else
            n =7;
        end
        l =length(tot);
        t = tot(1:n:l,1);
        ib= tot(2:n:l,1);
        b = tot(3:n:l,1);
        ks1= tot(4:n:l,1);
        ks2= tot(4:n:l,2);
        io= ib-ib./b;
        if n ==7
            ks3 = tot(5:n:l,1);
            ks4 = tot(5:n:l,2);
        end
    end
end
   

    
