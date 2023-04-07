%function m=hload(varargin);
%  load a Probability distribution from boots

function  m=hload(varargin);

ro =51;

if (length(varargin) == 1)&iscellstr(varargin)
    bfile = (char(varargin));
    hfile = strcat(bfile,'.hist');
    lfile = strcat(bfile,'.log');

    if fopen(hfile)+1 
        h = load(hfile);
    else 
        h = load(bfile);
    end
    n = 2;
    if sum(sum(h(:,3))>0)
        n = 3;
    end
    h(:,1) = h(:,1)*ro;
    h(:,2:end) = h(:,2:end)/ro;
    plot(h(:,1),h(:,n))
    axis([.2*ro 1.8*ro max(h(:,n))*-.2 max(h(:,n))*1.3])
    ylabel(bfile,'fontsize',8)
    xlabel('Distance','Fontsize',5)
    if n == 3
	ep = h(:,n)+h(:,4);
	em = h(:,n)-h(:,4);
	line(h(:,1),[ep em],'linestyle','--','color','blue')
    end
    m = sum(h(:,1).*h(:,n).*(h(2,1)-h(1,1)));
    xL = get(gca,'xlim');
    tp = xL(2) -(xL(2)-xL(1))*.3;
    text(tp,max(h(:,n))*1.05,sprintf('%s %.2f %s','x = ',m,' A'),'fontsize',8);
    var = sum(h(:,1).^2.*h(:,n)*(h(2,1)-h(1,1)))-m.^2;
    text(tp,max(h(:,n))*.90,sprintf('%s %.2f','var = ',var,' A^2'),'fontsize',8);
    if fopen(lfile)+1
        l = load(lfile)
        text(1.3,max(h(:,n))*.9,sprintf('%.1f %s',l(1,2),' s'),'fontsize',8);
    end
    fclose('all');
    
else
    if iscellstr(varargin)
        h = load(char(varargin(1)));
    else
        varargin = varargin{1};
        h = load(char(varargin(1)));
    end
    n = 2;
    ne = 2;
    if isnan(h(1,3))
        n = 3;
	ne =4;
    end
    xL = [h(1,1) h(end,1)];
    dx = h(2,1)-h(1,1);
    h = h(:,n);
    xall = xL;
    for i = 2:length(varargin)
	t = load(char(varargin(i)));
        xall(i,:)=[t(1,1) t(end,1)];
        xL = [min([xL t(1,1)]) max([xL t(end,1)])];
    end
    x = xL(1):dx:xL(2);
    h = [];
    for i = 1:length(varargin);
	t = load(char(varargin(i)));
        n = 2;
	ne = 2;
	if ~isnan(t(1,3))
		n = 3;
		ne=4;
	end	
	dt = t(2,1)-t(1,1);
        e = spline(t(:,1),t(:,ne),x);
        t = spline(t(:,1),t(:,n),x);
        rm = find(x<xall(i,1)|x>xall(i,2));
        t(rm) = 0;
        h(i,:) = t;%*dx/sum(t*dx);
        he(:,i) = t-e;
        he2(:,i) = t+e;
    end
    check = sum(h,2);
    x = x*ro;
    plot(x,h/ro,'linewidth',2)
    line(x,he/ro,'linestyle','--')
    line(x,he2/ro,'linestyle','--')
    legend(strrep(varargin,'\','\\'))
    set(gca,'FontSize',16,'LineWidth',1);
    xlim([10 90])
    xlabel('Distance (A)','fontsize',16)
    ylabel('Probability','fontsize',16)
    m=0;
end