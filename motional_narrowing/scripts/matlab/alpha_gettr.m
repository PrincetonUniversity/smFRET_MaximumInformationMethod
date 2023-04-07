function alpha_gettr
%function alpha_gettr
%   this function will return the time resolution for different alpha values
%   output format:
%	[alpha value, mean time resolution, std time resolution, total number of bins]



final_dir=sload('.final_dir');
work_dir=sload('.work_dir');
alphas=load('.alphas');
sample=sload('.sample');


for i = 1:length(alphas)

	if alphas(i) < 10 
		alpha_str = sprintf('0%1i',alphas(i));
	else
		alpha_str = sprintf('%2i',alphas(i));
	end
	

	cd(strcat(work_dir,'/',alpha_str,'/'))
	

	!cat *.xc > .xc.tot
	xc = load('.xc.tot');
	!rm .xc.tot
	xc = xc(:,2)-xc(:,1); %this is the length of each bin
	out(i,1:4)= [alphas(i)/100 mean(xc) std(xc)/length(xc) length(xc)];
end

cd(final_dir)

save(strcat(sample,'.mat'),'out')

function str=sload(fname)
fid = fopen(fname);
str = fgetl(fid);
fclose(fid);