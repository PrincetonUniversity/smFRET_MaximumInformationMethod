function process_global_fit

final_dir=sload('.final_dir');
alphas=load('.alphas');
sample=sload('.sample');
output_file=sload('.output_file');
r0 = load('.r0');

load([sample '.mat'])
tau=out(:,2);

if exist(strcat(output_file,'.mat'))
	load([output_file '.mat'])
	k0=parm(1);
	k1=parm(2);
	a=parm(3);
	b=parm(4);
	sigma=parm(5:end);

elseif exist(strcat(final_dir,'/',output_file,'.mat')) 
	load(strcat(final_dir,'/',output_file,'.mat'))
	fprintf('yes\n')
	k0=parm(1);
	k1=parm(2);
	a=parm(3);
	b=parm(4);
	sigma=parm(5:end);

else 
	k0=200
	k1=400
	a=0.4*r0;
	b=0.8*r0;
	sigma = 0.1639*ones(length(tau),1); 
end

while length(sigma) < length(alphas)
	sigma(end+1) = sigma(end);
end
if length(sigma) > length(alphas)
	sigma = sigma(1:length(alphas))
end

x0 = zeros(100,length(tau));
y0 = zeros(100,length(tau));
for i = 1:length(alphas)
	if alphas(i) < 10 
		alpha_str = sprintf('0%1i',alphas(i));
	else
		alpha_str = sprintf('%2i',alphas(i));
	end
	d1 = load(strcat(sample,alpha_str,'.out'));
	x0(:,i) = d1(:,1).*r0; % scale to angstrom
	y0(:,i) = d1(:,3)./(sum(d1(:,3))*(x0(2,i)-x0(1,i))); % renormalize to 1
end


[parm,exitf] = fit_GevaSkinner_model_global(k0,k1,tau,sigma,a,b,x0,y0);
if exitf==1
	save([output_file '.mat'],'parm');
end
  
wd=cd;
if strcmp(wd,final_dir)
	alpha_plot_results
end

exit



function str=sload(fname)
fid = fopen(fname);
str = fgetl(fid);
fclose(fid);
