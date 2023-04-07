function plot_alpha_series(bname)


load([bname '.mat']);
tr = out(:,2)*1000;

dists = dir([bname '*.out']);
dists = dists(end).name;
jj = strrep(dists,bname,'');
jj = strrep(jj,'.out','');
jj = str2num(jj);

for i = 1:length(tr)
	tmp = load(sprintf('%s%02i%s',bname,jj,'.out'));
	x = tmp(:,1);
	y(:,i) = tmp(:,3);
	jj=jj-1;
end

colormap('hsv');
h=waterfall(x',tr',y')
colormap('hsv');
figspec([2 9])
set(h,'LineWidth',2)
view(56,16)

function figspec(yL)
caxis(gca,[0 3]);
