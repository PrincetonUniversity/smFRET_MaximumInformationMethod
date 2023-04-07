function hist_error_gen


sample =sload('.sample');
alpha = load('.alphas');
boot_dir = sload('.boot_dir');
final_dir = sload('.final_dir');
boot_num = load('.boot_num');
r0=load('.r0'); 


for ai = 1:length(alpha)
	count=1;
	if alpha(ai) < 10
		alpha_num = ['0' num2str(alpha(ai))];
	else
		alpha_num = num2str(alpha(ai));
	end
	for ri = 1:boot_num
		boot_file=[boot_dir '/run' num2str(ri) '/' sample alpha_num '.out'];
		if exist(boot_file)
			tmp = load(boot_file);
            avedist(count) = sum(tmp(:,1).*tmp(:,3).*(tmp(2,1)-tmp(1,1))); 
			bdata(:,count) = tmp(:,3);
			count=count+1;
		end
	end
	tmp = load([sample alpha_num '.out']);
	data = [tmp(:,1) tmp(:,2) tmp(:,3) std(bdata')' tmp(:,5)];
	fn = [sample alpha_num '.e.out'];
    fdistmean = [sample alpha_num '.meanX'];
    distmean = [ro*mean(avedist) ro*std(avedist')];
	fprintf('%s\n',fn)
    fprintf('Average distance = %3.1f (%1.1f) A @ alpha%2.0f.\n',ro*mean(avedist),ro*std(avedist'),alpha(ai))
	save(fn,'data','-ASCII')
    save(fdistmean,'distmean','-ascii')
end

fprintf('\n%i bootstrapping runs were used in calculating error bars\n',count-1)

function str=sload(fname)

fid = fopen(fname);
str = fgetl(fid);
fclose(fid);