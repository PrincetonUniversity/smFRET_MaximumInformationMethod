function select_traj

dates = read_dates;


log_dir = sload('.log_dir',1);


nt = sload('.nt');
sigma = sload('.sigma');
sigmaI = 1.95;	%region I
corbinIII = sload('.corbinIII');
corbinII = sload('.corbinII');
corbinI = .002;	%region I


for date_cnt = 1:length(dates)
	
	good_dir  = [dates(date_cnt).dir 'good/'];
	
	good_savefile = [log_dir '/' dates(date_cnt).date '_pass_manual.txt'];
	bad_savefile = [log_dir '/' dates(date_cnt).date '_fail_manual.txt'];
	
	list = [];
	
	
	fprintf('Working on Directory %s\n============================================\n\n',good_dir);
	
	
	traj = read_traj(good_dir);

	if isempty(traj)
		fprintf('No good trajectories were identified on this day\n\n')
		%create an empty file anyway
		if 0%exist(good_savefile)
			[good_save,bad_save] = readf(good_savefile,bad_savefile);
			ii = 1;
			for i = 1:length(good_save)
				list(ii).good = good_save(i).good;
				ii = ii+1;
			end	
				
			for i = 1:length(good_save)
				list(ii).bad = bad_save(i).bad;
				ii = ii+1;
			end	
		end		
		write(list,good_savefile,bad_savefile);
			 
		continue;
	end
	

	keepold=0;
	empty = [];
	
	if exist(good_savefile,'file') 
		keepold = [];
		while isempty(keepold)
			keepold = input('Log from Previous File Selection Found\n   2 - view good trajectories only\n   1 - keep selections  \n   0 - delete selections\n:');
			if isempty(keepold)
				keepold = 100;
			end	
			switch keepold
				case 1
					[good_save,bad_save] = readf(good_savefile,bad_savefile);	
				case 2
					[good_save,bad_save] = readf(good_savefile,bad_savefile);
				case 0
					
				otherwise
					fprintf('Input was not understood!\n\n')
					keepold = [];
			end
		end		
	end
		
	t_cnt = 1;	


		
	while t_cnt <= length(traj)
		
		if keepold == 1
			found = test_old(traj(t_cnt).base,good_save,bad_save);
			if found == 1
				list(t_cnt).good = traj(t_cnt).base;
				list(t_cnt).bad = [];
				write(list,good_savefile,bad_savefile);
				t_cnt = t_cnt+1;
				continue;
			end
			if found == 2
				list(t_cnt).good = [];
				list(t_cnt).bad = traj(t_cnt).base;
				write(list,good_savefile,bad_savefile);
				t_cnt = t_cnt+1;
				continue;
			end
		end	
		
		if keepold == 2
			found = test_old(traj(t_cnt).base,good_save,bad_save);
			if found == 2
				list(t_cnt).good = [];
				list(t_cnt).bad = traj(t_cnt).base;
				write(list,good_savefile,bad_savefile);
				t_cnt = t_cnt+1;
				continue;
			end
		end	
		
		cont = 0;
		plot_type = 1;
		traj_zoom=[];
		while cont == 0
			cla;
			clf;
						
			%================
			mean_time_res = load(traj(t_cnt).xc);
			if length(mean_time_res) == 0
				list(t_cnt).good = [];
				list(t_cnt).bad = traj(t_cnt).base;
				write(list,good_savefile,bad_savefile);	
				cont = 1;05
				continue;
			elseif size(mean_time_res,1) == 1
				mean_time_res = mean_time_res(2) - mean_time_res(1);
			else	
				mean_time_res = mean(diff(mean_time_res(:,1)))	;
			end	
			subplot(3,2,1)
			iload(traj(t_cnt).fr,traj(t_cnt).log,mean_time_res*10)
			title(traj(t_cnt).name)
			
	
			%================
			sp = 3;
			if plot_type == 1
				sp = 2;
			end
			subplot(3,2,sp)
			log = load(traj(t_cnt).log);
			if isempty(traj_zoom)
				traj_zoom=log(1,1)*1.1;
			end	
			xL = [0 traj_zoom];
			iload(traj(t_cnt).fr,traj(t_cnt).log,mean_time_res*4)
			set(gca,'xlim',xL);
			title([])
			if plot_type == 1
				title(sprintf('%i of %i files',t_cnt,length(traj)))
			end	
			
			%================e
			sp = 5;
			if plot_type == 1
				sp = 6;
			end	
			subplot(3,2,sp)
			xload(traj(t_cnt).xc)
			set(gca,'xlim',[0 log(1,2)*1.1]);
			title(good_dir)
			
			switch plot_type
				case 1
					subplot(3,2,4)
					reg_I_plot(traj(t_cnt),corbinI,sigmaI,nt)
					
					subplot(3,2,3)
					cor_plot(traj(t_cnt),corbinII,2)
					title('Region II Auto-correlation')
					
					subplot(3,2,5)
					cor_plot(traj(t_cnt),corbinIII,3)
					title('Region III Auto-correlation')
									
				case 2
					coplot(traj(t_cnt).sam);
					
					subplot('position',[.55 .7 .4 .2])
        				ipplot(traj(t_cnt).fr);
        				title(sprintf('%i of %i files',t_cnt,length(traj)));
        				
        			case 3
					subplot('position',[.55 .05 .4 .4])
					if exist(traj(t_cnt).hist)
						hload(traj(t_cnt).hist,[good_dir 'hist10.out']);
					end	
					title([])
					legend off
					axis square
					
					subplot('position',[.55 .55 .4 .4])
					imloads(traj(t_cnt).base)
					title(sprintf('%i of %i files',t_cnt,length(traj)))
			end
			
			fprintf('# %i of %i files: %s\n',t_cnt,length(traj),traj(t_cnt).base)
			fprintf('Choose one of the following\n')
			fprintf('===========================\n')
			fprintf('  enter - keep trajectory\n')
			fprintf('  0     - discard trajectory\n')
			fprintf('  1     - Display correlation results\n')
			fprintf('  2     - Display cox-oakes test\n')
			fprintf('  3     - Display distribution\n')
			fprintf('  4     - Change trajectory zoom\n')
			fprintf('  7     - exit\n')
			fprintf('  9     - go back\n')
			test = input(':');

			if isempty(test)
				test = 100;
			end
			
			switch test
				case 1
					plot_type=1;
				case 2
					plot_type=2;
				case 3
					plot_type=3;
				case 4
					traj_zoom=abs(input('Input new zoom in seconds\n:'))
						
				case 100
					list(t_cnt).good = traj(t_cnt).base;
					list(t_cnt).bad = [];
					write(list,good_savefile,bad_savefile);
					cont = 1;
				case 0
					list(t_cnt).good = [];
					list(t_cnt).bad = traj(t_cnt).base;
					write(list,good_savefile,bad_savefile);					
					cont = 1;
				case 9
					num = [];
					while isempty(num)
						num = abs(input('Go Back how many files?\n:'));
						if num > 0;
							t_cnt = t_cnt-num;
							if t_cnt < 1
								t_cnt = 1;
							end
						else
							fprintf('Input not understood!\n\n')
							num = [];
						end
					end
					traj_zoom=[];
				case 7
					return;		
				otherwise	
					fprintf('Input not understood!\n\n')
			end				
					
		end
		t_cnt = t_cnt+1;
		if ~isempty(empty)
			t_cnt = empty(1);
			empty(1) = [];
		end
	end
	!/home/$USER/Projects/alpha/scripts/sort_cleanup.sh
end




function dates = read_dates
fid = fopen('data_dirs.txt');
tmp = fgetl(fid);
i = 1;
while tmp ~= -1
	dates(i).dir = tmp;
	dates(i).date = tmp(18:25);
	tmp = fgetl(fid);
	i = i+1;
end
fclose(fid);


function traj = read_traj(good_dir)
traj = [];

pass_file = [good_dir 'pass_corr.txt'];

if exist(pass_file)
	fid = fopen(pass_file);
	tmp = fgetl(fid);
	i = 1;
	while tmp ~= -1
		traj(i).name  = tmp;
		traj(i).fr= [good_dir tmp];
		traj(i).log = strrep(traj(i).fr,'.fr','.log');
		traj(i).xc = strrep(traj(i).fr,'.fr','.xc');
		traj(i).sam = strrep(traj(i).fr,'.fr','.sam');
		traj(i).hist = strrep(traj(i).fr,'.fr','.hist');
		traj(i).base = strrep(traj(i).fr,'.fr','');
		tmp = fgetl(fid);
		i = i+1;
	end
	fclose(fid);
else
	traj = [];
end	



function reg_I_plot(traj,corbinI,sigma,nt)

d_auto = frcorr(traj.fr,corbinI,traj.log,5);
a_auto = frcorr(traj.fr,corbinI,traj.log,6);
			
d_test = test_auto(d_auto,sigma,nt);
a_test = test_auto(a_auto,sigma,nt);
opts = optimset('Display','off');
tau_guess = [length(d_auto)*3/4 length(d_auto)/3 100 5];
mt = [];
if d_test

	for j = 1:length(tau_guess)

		[parms(:,j),fval(j),exitflag(j)] = fminsearch(@exp_fit,[tau_guess(j) d_auto(1,2)],opts,d_auto);
			
	end
	j = find(fval == min(fval));
	parms = parms(:,j);
	parms(1) = floor(abs(parms(1)));
	mt(1) = parms(1);
end
tau_guess = [length(a_auto)*3/4 length(a_auto)/3 100 5];
if a_test
	
	for j = 1:length(tau_guess)

		[parms(:,j),fval(j),exitflag(j)] = fminsearch(@exp_fit,[tau_guess(j) a_auto(1,2)],opts,a_auto);

	end	
		j = find(fval == min(fval));
	parms = parms(:,j);
	parms(1) = floor(abs(parms(1)));
	mt(2) = parms(1);
end

cross = frcorr(traj.fr,corbinI,traj.log,7);
if sum(mt)
	cbin = max(mt);
	cross_e = frcorr(traj.fr,corbinI,traj.log,8,cbin);
	cross(:,3) = cross_e*ones(size(cross(:,3)));
else
	cbin = nt;
end

c_test = test_cross(cross,sigma,cbin);

xL = [cross(1,1) cross(end,1)];
if xL(2) == xL(1)
	xL(2) = xL(2)*2;
end
plot(cross(:,1),cross(:,2)./cross(:,3),'b',xL,[2 2],'k--',xL,[-2 -2],'k--',xL,[0 0],'k:')
xlabel('tau')
ylabel('C_x_y/\sigma')
title('Region I cross-correlation')
yL = get(gca,'ylim');
yL = get(gca,'ylim');
if yL(1) > -3
	yL(1) = -3;
end
if yL(2) < 3
	yL(2) = 3;
end

set(gca,'xlim',xL,'ylim',yL)

if c_test == -1
	text(xL(2)*.1,2.15,'Significant negative-correlation to 95% confidence')
elseif c_test == 1
	text(xL(2)*.1,-2.25,'Significant positive-correlation to 95% confidence')
	line(xL,yL,'color','r')
	line(xL,[yL(2) yL(1)],'color','r')
end		

function fail = test_auto(data,sigma,nt)
if size(data,1) < nt
	nt = size(data,1);
end	
fail = 0;

Z = abs(sum(data(1:nt,2))/nt);
C = sigma*data(1,3)/sqrt(nt);
if Z>C
	fail=1;
end
	

function out = exp_fit(parms,data)
parms(1) = floor(abs(parms(1)));
efit = exp_gen(parms,data(:,1));
out = sum((efit' - data(:,2)).^2);
	
function efit = exp_gen(parms,x)
if parms(1) == 0
	efit = zeros(size(x));
	efit = efit';
else
	x = 1:length(x);
	efit = parms(2)	* exp(-x./parms(1));
end
	

function fail = test_cross(data,sigma,nt)
if size(data,1) < nt
	nt = size(data,1);
end	

Z = sum(data(1:nt,2))/nt;
C = sigma*data(1,3)/sqrt(nt);
	
fail = 0;	
if Z > C
	fail = 1;
elseif Z < -C
	fail = -1;
end	
	

function cor_plot(traj,corbin,region)

if region == 2
	Dreg = 1;
	Areg = 3;
else
	Dreg = 2;
	Areg = 4;
end	

donor = frcorr(traj.fr,corbin,traj.log,Dreg);
acceptor = frcorr(traj.fr,corbin,traj.log,Areg);

if size(donor,1) < 2
	return
end	

xL = [donor(1,1) donor(end,1)];

plot(donor(:,1),donor(:,2)/donor(1,3),'b',acceptor(:,1),acceptor(:,2)/acceptor(1,3),'r',...
	xL,[2 2],'k--',xL,[-2 -2],'k--',xL,[0 0],'k:')
yL = get(gca,'ylim');
if yL(1) > -3
	yL(1) = -3;
end
if yL(2) < 3
	yL(2) = 3;
end		


set(gca,'xlim',xL,'ylim',yL)


xlabel('tau')
ylabel('C_x_x/\sigma')	


function [good,bad] = readf(gname,bname)
good = [];
bad = [];

empty = [];

gid = fopen(gname);
i = 1;
to = fgetl(gid);
while to~=-1
    if ischar(to)
    	good(i).good = to;
    	i = i+1;
    end
    to = fgetl(gid);
end
fclose(gid);

if exist(bname)
	i = 1;
	gid = fopen(bname);
	to = fgetl(gid);
	while to~=-1
		if ischar(to)
			bad(i).bad = to;
			i = i+1;
		end
		to = fgetl(gid);
	end
	fclose(gid);
end	


function found = test_old(fr,good,bad)
found = 0;
for i = 1:length(good)
	if strcmp(good(i).good,fr)
		found = 1;
		return
	end
end
for i = 1:length(bad)
	if strcmp(bad(i).bad,fr)
		found = 2;
		return
	end
end		


function skip = testan(an,fname);
skip = 0;
for i = 1:length(an);

    if(length(double(char(an(i))))==length(double(fname)))
        if (double(char(an(i)))==double(fname))
        skip=1;
        end
    end
end

function write(fr,gname,bname)
gid = fopen(gname,'w');
bid = fopen(bname,'w');
for i = 1:length(fr)
	if ~isempty(fr(i).good)
		fprintf(gid,'%s\n',fr(i).good);
	elseif ~isempty(fr(i).bad)
		fprintf(bid,'%s\n',fr(i).bad);		
	end
end
fclose('all');


function coplot(fname)
if exist(fname)
	fid = fopen(fname);
	for i = 1:6
	cotmp = [];
	in =1;
	tmp = fgetl(fid);
	while(~isnumeric(tmp))
		tmp = str2num(tmp);
		if(isempty(tmp))
		break;
		end
		cotmp(in) = tmp;
		tmp = fgetl(fid);
		in = in+1;
	end
	switch i
		case 1
		cofiles.b10 = cotmp;
		case 2
		cofiles.b50 = cotmp;
		case 3
		cofiles.b100 = cotmp;
		case 4
		cofiles.b500 = cotmp;
		case 5
		cofiles.b1000=cotmp;
		case 6
		cofiles.b5000=cotmp;
	end
	end
	fclose(fid);
	
	subplot('position',[.55 .025 .4 .1])
	draw_coplot(cofiles.b5000,'5000')
	
	subplot('position',[.55 .125 .4 .1])
	draw_coplot(cofiles.b1000,'1000')
	set(gca,'xticklabel',[])
	
	subplot('position',[.55 .225 .4 .1])
	draw_coplot(cofiles.b500,'500')
	set(gca,'xticklabel',[])
	
	subplot('position',[.55 .325 .4 .1])
	draw_coplot(cofiles.b100,'100')
	set(gca,'xticklabel',[])
	
	subplot('position',[.55 .425 .4 .1])
	draw_coplot(cofiles.b50,'50')
	set(gca,'xticklabel',[])
	
	subplot('position',[.55 .525 .4 .1])
	draw_coplot(cofiles.b10,'10')
	set(gca,'xticklabel',[])
end

function draw_coplot(data,pnum)
cla;

if ~isempty(data)
    bw = 3.49*std(data)*length(data)^(-1/3);
    hx=min(data):bw:max(data);
    if length(hx)<4
        hx = 4;
    end

    [y,x]=hist(data,hx);
    bar(x,y)
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','w','EdgeColor','k')
end
set(gca,'xlim',[-17 17])
yL = get(gca,'ylim');
set(gca,'yLim',[0 yL(2)])
text(-16,yL(2)*.9,sprintf('%s %s %s','binsize =',pnum,'photons'))


function ipplot(fname)

logf = load(strrep(fname,'.fr','.log'));
ts = logf(1,2);
tb = logf(1,1);

binsize = 100/logf(2,1);
if (((tb-ts)/binsize < 500)&(logf(2)*40 < 100))
	binsize = (tb - ts)/40;
end

I = pload(fname);

ib=I(:,1);
ib=ib(find(ib));
dp=diff(ib);
turnovers=find(dp<-10);
while(length(turnovers>0))
    ib(turnovers(1)+1:length(ib))=ib(turnovers(1)+1:length(ib))+2^32*12.5e-9;
    turnovers(1)=[];
end

ib = ib(find((ib>ts)&(ib<=tb)));

[y,x] = hist(ib,ib(end)/binsize);

Ib = mean(y);

bw = 0:max(y);
[y,x] = hist(y,bw);


xp = 0:round(max(x));
for i = 1:length(xp)
    fact = factorial(xp(i))*(x(2)-x(1));
    pd(i) = exp(-Ib)*exp(log(Ib)*xp(i))./fact;
end

pd = pd*sum(y); 
bar(x,y)
h = findobj(gca,'Type','patch');
set(h,'FaceColor','w','EdgeColor','k')
line(xp,pd,'linestyle','--','color','red')
xL = get(gca,'xlim');
yL = get(gca,'ylim');
text([abs(.1*xL(1))+xL(1)],yL(2)*.9,sprintf('%s %.1f %s','binsize =',binsize*1000,'ms'))


function num=sload(fname,str_test)
fid = fopen(fname);
str = fgetl(fid);
if nargin == 1
	num = str2num(str);
else
	num = str;
end		
fclose(fid);
