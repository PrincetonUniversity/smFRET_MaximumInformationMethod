function alpha_plot_results(results)

final_dir=sload('.final_dir');
work_dir=sload('.work_dir');
alphas=load('.alphas');
sample=sload('.sample');
boot_num=str2num(sload('.boot_num'));
boot_dir=sload('.boot_dir');
output_file=sload('.output_file');
r0=load('.r0');

load([sample '.mat']);
tau=out(:,2);

if nargin == 0
	results=output_file;
end
load([results '.mat']);
k0=parm(1);
k1=parm(2);
a=parm(3);
b=parm(4);
sigma = parm(5:length(tau)+4);
s2 = sigma.*sigma*2;

%==============================================================

for i = 1:length(tau)
	figure(1)
	if alphas(i) < 10 
		alpha_str = sprintf('0%1i',alphas(i));
	else
		alpha_str = sprintf('%2i',alphas(i));
	end
	
	hist_file = load([sample alpha_str '.out']);
	x0(i,:) = hist_file(:,1)*r0;
	y0(i,:) = hist_file(:,3)./(sum(hist_file(:,3))*(x0(i,2)-x0(i,1)));
	dist = GevaSkinner(k0,k1,tau(i)); 

	dx = dist(2,1)-dist(1,1);
	nl = floor(abs(-2.0/dx)); 
	nr = floor(abs(2.0/dx));
	xl = dist(1,1) - (nl:-1:1)'.*dx;
	xr = dist(length(dist(:,1)),1) + (1:nr)'.*dx;
	x = [xl; dist(:,1); xr];
	
	g = exp(-x.^2./s2(i))./sqrt(2*pi)./sigma(i);
	gg = conv(dist(:,2),g);
	model = gg(1:length(x));
	
	xm = a.*x + b;
	ym = interp1(xm, model, x0(1,:)); % use interpolation
	ym = ym / (sum(ym)*(x0(i,2)-x0(i,1))); % renormalize
	
	
	subplot(3,3,i)
	plot(x0(i,:),y0(i,:),'o',x0(i,:),ym,'r')
	title([sample alpha_str '.out' ' - ' num2str(ceil(tau(i)*10000)/10) ' ms']) 
	if i ==1;
	xL = get(gca,'xlim');
	yL = get(gca,'ylim');
	else
	set(gca,'xlim',xL,'ylim',yL);
	end
	
	if exist(['data/' alpha_str '/ctb.mat'])
		figure(4)
		load(['data/' alpha_str '/ctb.mat'])
		
		FRET = r0*(1./FRET-1).^(1/6);
		
		xx = linspace(0,2*r0,41);
		yy = hist(FRET,xx);
		yy = yy ./ sum(yy) ./ (xx(2)-xx(1));
	
		subplot(3,3,i)
		bar(xx,yy,'w')
		line(x0(i,:),hist_file(:,2)/r0,'linestyle','--','color','b')
		line(x0(i,:),y0(i,:),'linestyle','-','marker','o','color','b')
		axis tight
		title([num2str(ceil(tau(i)*10000)/10) ' ms'])
	end
	i = i+1;
end
print('-djpeg','Figure_1.jpg')

%=============================================================

figure(3)
plot(x0(1,:),y0,'b')
print('-djpeg','Figure_3.jpg')

zL = get(gca,'yLim');

%==============================================================


figure(2)
j=15;
colormap('hsv');
h=waterfall(x0(1,:),tau*1000,y0);
set(gca,'xlim',[.2 1.8]*r0,'ylim',[tau(1)*1000-.5 tau(end)*1000+1],'zLim',zL)
xlabel('Distance (A)')
ylabel('Time Resolution (ms)')
zlabel('Probability Denisity')
set(gca,'LineWidth',1.5,'FontSize',18)
set(h,'LineWidth',2)
view(56,16)
print('-djpeg','Figure_2.jpg')


%==============================================================s

if exist([boot_dir '/run1/' output_file '.mat'])
	boot_count = 1;
	for bi = 1:boot_num
		if exist([boot_dir '/run' num2str(bi) '/' output_file '.mat'])
			load([boot_dir '/run' num2str(bi) '/' output_file '.mat']);
			bparm(:,boot_count) = parm;
			boot_count=boot_count+1;
		else %exit if the next run is not found
			break
		end	
	end
	sparm = std(bparm');
	fprintf(1,'  kopen    = %8.3f (%.3f)\n', k0,sparm(1));
	fprintf(1,'  kclose   = %8.3f (%.3f)\n', k1,sparm(2));
	fprintf(1,'  xopen    = %8.3f (%.3f)\n', b+a,sqrt(sparm(3).^2+sparm(4).^2));
	fprintf(1,'  xclose   = %8.3f (%.3f)\n\n', b,sparm(4));
	for k=1:length(tau),
		fprintf(1,'  sigma(%i) = %8.3f (%.3f)\n', k, sigma(k)*a,sqrt(sparm(k+4).^2+sparm(3).^2));
	end
	
	fprintf('\n  fopen   = %8.3f (%.3f)\n',k0/(k0+k1),std(bparm(1,:)./(bparm(1,:)+bparm(2,:))))
	fprintf('  fclose  = %8.3f (%.3f)\n',k1/(k0+k1),std(bparm(2,:)./(bparm(1,:)+bparm(2,:))))
	
	fprintf(1,'\n  %i bootstraping runs used\n\n',boot_count-1)
else
	fprintf(1,'  kopen    = %8.3f\n', k0);
	fprintf(1,'  kclose   = %8.3f\n', k1);
	fprintf(1,'  xopen    = %8.3f\n', b+a);
	fprintf(1,'  xclose   = %8.3f\n\n', b);
	for k=1:length(tau),
		fprintf(1,'  sigma(%i) = %8.3f\n', k, sigma(k)*a);
	end



	fprintf('\n  fopen   = %8.3f\n',k0/(k0+k1))
	fprintf('  fclose  = %8.3f\n',k1/(k0+k1))
	end

function str=sload(fname)
fid = fopen(fname);
str = fgetl(fid);
fclose(fid);
