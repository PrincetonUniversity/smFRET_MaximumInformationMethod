function bootstrap_check_runs


sample =sload('.sample');
alpha = load('.alphas');
boot_dir = sload('.boot_dir');
final_dir = sload('.final_dir');
boot_num = load('.boot_num');
work_dir = sload('.work_dir');
output_file = sload('.output_file');


cd(final_dir)
for ai = 1:length(alpha)
	if alpha(ai) < 10
		alpha_num = ['0' num2str(alpha(ai))];
	else
		alpha_num = num2str(alpha(ai));
	end

	tmp = load([sample alpha_num '.out']);
	full(ai).h = tmp(:,[1 3]);
end
load([output_file '.mat']);
fullr = parm;
cd(boot_dir)
ri=1;
run_num = ['run' num2str(ri)];
while exist(run_num)
	cd(run_num)
	
	for ai = 1:length(alpha)
	if alpha(ai) < 10
		alpha_num = ['0' num2str(alpha(ai))];
	else
		alpha_num = num2str(alpha(ai));
	end
		tmp = load([sample alpha_num '.out']);
		boot(ai).h = tmp(:,[1 3]);
	end
	
	for pi = 1:length(alpha)
		subplot(3,3,pi)
		plot(full(pi).h(:,1),full(pi).h(:,2),boot(pi).h(:,1),boot(pi).h(:,2))
	end
	fprintf('Run %i\n',ri)
	if exist([output_file '.mat'])
		load([output_file '.mat']);
		bootr = parm;
		print_results(fullr,bootr)
	end
	pause
	ri=ri+1;
	run_num = ['run' num2str(ri)];
	cd(boot_dir)
end
cd ../../

function print_results(fullr,bootr)
d = bootr-fullr;
fprintf('\n%15s %15s %15s %15s\n','parameter','full dataset','bootstrap set','difference')
fprintf('%15s %15s %15s %15s\n','----------','---------------','---------------','---------------')
fprintf('%15s %15.0f %15.0f %15.0f\n','kopen',fullr(1),bootr(1),d(1))
fprintf('%15s %15.0f %15.0f %15.0f\n','kclose',fullr(2),bootr(2),d(2))
fprintf('%15s %15.2f %15.2f %15.2f\n','xopen',fullr(3)+fullr(4),bootr(3)+bootr(4),d(3)+d(4))
fprintf('%15s %15.2f %15.2f %15.2f\n\n','xclose',fullr(4),bootr(4),d(4))

for i = 5:length(fullr)
	sig = ['sigma(' num2str(i-4) ')'];
	fprintf('%15s %15.2f %15.2f %15.2f\n',sig,fullr(i)*51,bootr(i)*51,d(i)*51)
end


function str=sload(fname)
fid = fopen(fname);
str = fgetl(fid);
fclose(fid);
