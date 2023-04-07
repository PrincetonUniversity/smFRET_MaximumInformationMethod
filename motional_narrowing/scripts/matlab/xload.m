% plot xc files
function xload(fname,xname,f,nice)

if(~exist('nice'))
	nice=0;
end

cla reset;
gotx=0;
if(nargin==0)
    fname='x.xc';
    xname='x.txt';
    gotx=1;
end

if(nargin>1)
	if(ischar(xname))
		gotx=1;
	end
end

if(nargin<3)
    f=1;
end

if(nargin<4)
	nice=0;
end

if(~exist(fname))
	return;
end
    
ah=gca;
hold on;
set(ah,'Box','on');
x=load(fname);
if(length(x)==0)
    return;
end
tr = (x(end,2)-x(1,1))/length(x);
x=x';

    
r0=51;
if nice
    r0 = 51;
end
x(3:4,:)=x(3:4,:)*r0;

t1=x(1,:);
t2=x(2,:);
tc=t1/2+t2/2;
xc=f*x(3,:);
sc=x(4,:);
xm=f*x(3,:)-f*x(4,:);
xp=f*x(3,:)+f*x(4,:);

xl=[t1;t2;t2;t1;t1];
yl=[xm;xm;xp;xp;xm];
lineh=patch(xl,yl,[.7 .7 .7]);
set(lineh,'LineStyle','none');
lineh=plot(tc,xc,'k-');
if(nice)
    set(lineh,'LineWidth',1);
    set(gca,'FontSize',18,'Linewidth',1.5)
    
    ylabel('Distance (A)','FontSize',18);
    xlabel('time (s)','FontSize',18);
end
ylabel('Distance (A)');
xlabel('time (s)');
Lx = get(gca,'xlim');
Ly = get(gca,'ylim');
str = sprintf('%s %.1f %s','t =',tr*1000,' ms');
if ~nice
    text((Lx(2)-Lx(1))*.45+Lx(1),(Ly(2)-Ly(1))*.9+Ly(1),str);
end

if(gotx)
    xt=load(xname);
    lineh=plot(xt(:,1),xt(:,2),'k-');
    set(lineh,'LineWidth',1);
end


hold off;
