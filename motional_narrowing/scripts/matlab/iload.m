%   IN=iload(FILENAME) loads the time and intensity data into the variable IN.

function [out,coor]=iload(fname,logfile,binsize,nice)

if(nargin<4)
	nice=0;
end

if(nargin==0)
    fname='FRET.fr';
    binsize=0.025;
end
cid=fopen(fname,'r','ieee-le');

if(nargin==2)
    binsize=0;
end

if nargin == 1
    binsize = .025;
end

fseek(cid,0,'bof');
if(fread(cid,1,'*uint32')~=hex2dec('d3edf5f2'))
    error(cat(2,fname,' is not a smurf file.'))
end

cla;
x2 = 0;
y2 = 0;

ver=fread(cid,1,'*uint16');
ft=fread(cid,1,'*uint16');
if(ft==2)
    I=id(cid);
elseif(ft==4)
    I=ii(cid);
elseif(ft==6)
    [I pow coor]=fret(cid);
    %   nbin=sum(I(:,1))./max(I(:,1))
    
    if(ver==1)
        if binsize == 0 
            binsize = .05;
        end
        nbin=sum(I)*12.5e-9/binsize;
        binl=sum(I)./(nbin*80e6);
        [y1 x1]=hist(cumsum(I(find(I(:,1)),1)/80e6),nbin(1));
        y1=y1(1:length(y1)-1)/1000/binl(1);
        x1=x1(1:length(x1)-1);
        
        if(size(I,2)>1)
            [y2 x2]=hist(cumsum(I(find(I(:,2)),2)/80e6),nbin(1));
            y2=y2(1:length(y2)-1)/1000/binl(1);
            x2=x2(1:length(x2)-1);
        end

    elseif(ver==2)
        I=I*12.5e-9;
        ib=I(:,1);
        ib=ib(find(ib));
        dp=diff(ib);
        turnovers=find(dp<-10);
        while(length(turnovers>0))
            ib(turnovers(1)+1:length(ib))=ib(turnovers(1)+1:length(ib))+2^32*12.5e-9;
            turnovers(1)=[];
        end
        if binsize == 0;
            l1 = load(logfile);
            binsize = 50/l1(2,1);
        end
        [y1 x1]=hist(ib,ib(length(ib))/binsize);
        y1=y1/1000/binsize;
        if size(I,2)>1
        	ib=I(:,2);
        	ib=ib(find(ib));
            dp=diff(ib);
            turnovers=find(dp<-10);
            while(length(turnovers>0))
                ib(turnovers(1)+1:length(ib))=ib(turnovers(1)+1:length(ib))+2^32*12.5e-9;
                turnovers(1)=[];
            end
            [y2 x2]=hist(ib,ib(length(ib))/binsize);
            y2=y2/1000/binsize;
        end
    end
    handb=plot(x1,y1,'b-');
	hold on;
	handr=plot(x2,y2,'r-');
	hold off;
	hand=[handb handr];
    if(nice)
        set(hand,'LineWidth',1,'LineStyle','-');
        set(gca,'FontSize',16,'LineWidth',1.5);
        xlabel('time (s)','FontSize',16);
        ylabel('Intensity (kcps)','FontSize',16);
        
    end
    xlabel('time (s)');
    ylabel('Intensity (kcps)');
    title(fname);
    if ~isempty(x1)
    	axis([0 max([max(x1); max(x2)])*1.05 0 max([max(y1); max(y2)])*1.05]);
    end	
    v=axis;
    if((nargin>1)&&(~isempty(logfile)&&(exist(logfile))))
        l1=load(logfile);
        if(~isempty(l1))
            di=[l1(2,1);l1(2,1)]/1000;
            ai=[l1(3,1);l1(3,1)]/1000;
            xl=[v(1);v(2)];
            h=line(xl,di);
            set(h,'LineWidth',1,'LineStyle','--','Color','blue');
            h=line(xl,ai);
            set(h,'LineWidth',1,'LineStyle','--','Color','red');
            dl=[l1(1,1);l1(1,1)];
            al=[l1(1,2);l1(1,2)];
            yl=[v(3);v(4)];
            h=line(al,yl);
            set(h,'Color','red','LineStyle','--');
            if(nice)
                set(h,'LineWidth',1,'LineStyle','--');
            end
            h=line(dl,yl);
            set(h,'Color','blue','LineStyle','--');
            if(nice)
                set(h,'LineWidth',1,'LineStyle','--');
            end
        end
    end
    
    
    if(nargout~=0)
	out=1;
    end
    return;
elseif((ft==6)&&(ver==2))
else
    error(cat(2,fname,' is not a time trajectory file.'))
end

fseek(cid,48,'bof');
p=fread(cid,1,'double');

fclose(cid);

if(nargout==0)
    plot(I(:,1),I(:,2),'b');
    ylabel('Intensity (cps)')
    xlabel('Time (s)')
    title(strcat('Power: ',num2str(p),'\muW'))
else
    out=I;
end

%----------------------------
function out=id(cid)

fseek(cid,24,'bof');

fseek(cid,0,'eof');
N=(ftell(cid)-1024)/4;

fseek(cid,1024,'bof');
c=fread(cid,N,'uint32');

nbin=sum(c)/max(c);

[y,x]=hist(cumsum(c/80e6),nbin);

I=y*80e6/max(c);
t=x;

out=[t' I'];

%------------------------------
function out=ii(cid)

fseek(cid,24,'bof');
dt=fread(cid,1,'*double');

fseek(cid,0,'eof');
N=(ftell(cid)-1024)/8;

fseek(cid,1024,'bof');
c=fread(cid,N,'double');

t=[dt:dt:N*dt];

out=[t' c];
%-------------------------------
function [out,pow,coor]=fret(cid)

N=0;
L=zeros(4,1);
for i=1:4
    fseek(cid,4+4*i,'bof');
    L(i)=fread(cid,1,'uint32');
    if L(i)
        N=N+1;
    else 
        break;
    end
end


fseek(cid,128,'bof');
act=fread(cid,4,'double');
meas=fread(cid,4,'double');

maxl=max(L);
t=[1:1:maxl];
d=zeros(maxl,N);

fseek(cid,48,'bof');
pow=fread(cid,1,'double');

fseek(cid,56,'bof');
coor=fread(cid,2,'double');

fseek(cid,1024,'bof');
for i=1:N
    dt=fread(cid,L(i),'uint32');
    d(1:L(i),i)=dt;
end

out=d;
