function sort_corsort
global nt 

nt = str2num(sload('.nt'));
sigma = str2num(sload('.sigma'));
corbinIII = str2num(sload('.corbinIII'));
corbinII = str2num(sload('.corbinII'));


if ~exist('pass_sn.txt')
	exit
end

fr = read_files('pass_sn.txt');

prerun = 0;
if exist('pass_corr.txt')
	good = read_files('pass_corr.txt')
	prerun=1;
else
	good = [];
end		
if exist('fail_corr.txt')
	fail = read_files('fail_corr.txt')
	prerun=1;
else
	fail = [];
end

cor = dir('*.cor');

for i = 1:length(fr)
	if prerun
		exists = check_run(good,fail,fr(i).name);
		if exists
			continue
		end	
	end

	fr(i).name
	auto = frcorr(fr(i).name,corbinII,strrep(fr(i).name,'.fr','.log'),1);
	if test_auto(auto,sigma,nt)
		fail(end+1).name = fr(i).name;
		continue
	end	

	auto = frcorr(fr(i).name,corbinII,strrep(fr(i).name,'.fr','.log'),3);
	if test_auto(auto,sigma,nt)
		fail(end+1).name = fr(i).name;
		continue
	end
	
	auto = frcorr(fr(i).name,corbinIII,strrep(fr(i).name,'.fr','.log'),2);
	if test_auto(auto,sigma,nt)
		fail(end+1).name = fr(i).name;
		continue
	end	

	auto = frcorr(fr(i).name,corbinIII,strrep(fr(i).name,'.fr','.log'),4);
	if test_auto(auto,sigma,nt)
		fail(end+1).name = fr(i).name;
		continue
	end	

	good(end+1).name = fr(i).name;

end



write_file_list(good,'pass_corr.txt');
write_file_list(fail,'fail_corr.txt');


function fr = read_files(fname)
fr =[];
fid=fopen(fname);

tmp = fgetl(fid);
cnt = 1;
while tmp ~= -1
	fr(cnt).name = tmp;
	tmp = fgetl(fid);
	cnt = cnt+1;
end
fclose(fid);


function write_file_list(fnames,save_name)
fid=fopen(save_name,'w');

for i = 1:length(fnames)
	fprintf(fid,'%s\n',fnames(i).name);
end
fclose(fid);


function str=sload(fname)
%reads the first line of a .config file
fid = fopen(fname);
str = fgetl(fid);
fclose(fid);


function exists=check_run(good,fail,fname)
exists = 0;
for i = 1:length(good)
	if strcmp(good(i).name,fname)
		exists = 1;
		return
	end
end
for i = 1:length(fail)
	if strcmp(fail(i).name,fname)
		exists = 1;
		return
	end
end	
	

function fail = test_auto(data,sigma,nt)
	
	if size(data,1) < 4
		fail =0;
		return
	end
	
	if size(data,1) < nt
		nt = size(data,1);
	end	
	fail = 0;
	
	Z = abs(sum(data(1:nt,2))/nt);
	C = sigma*data(1,3)/sqrt(nt);

	if Z>C
		fail=1;
	end
