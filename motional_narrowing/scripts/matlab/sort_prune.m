function sort_prune


for i = 1:4
	if exist(['m' num2str(i)]) == 7
		cd(['m' num2str(i)])
		sort_fn
		cd ../
	end	
end	



function sort_fn

ratiosigma = 5;

names = {};


peran;

names = readnames;

[x,o] = sort(iod);
y = ioa(o);
xl = [x(1) x(end)];

total = length(x);

ta = ta(o);
t = ta;
names = names(o);

bd = bd(o);
ba = ba(o);


for i = 1:length(names)
    xc = load(char(strrep(names(i),'.log','.xc')));
    binn(i) = length(xc);
    bint(i) = mean(xc(:,2)-xc(:,1));
end
    

if ratiosigma

	del=length(x)*sum(x.^2)-sum(x)^2;
	b = (sum(x.^2)*sum(y)-sum(x)*sum(y.*x))/del;
	m = (length(x)*sum(x.*y)-sum(x)*sum(y))/del;
	
	dy = sqrt((length(x)-2)^-1*sum((y-b-m*x).^2));
	
	db=dy*sqrt(sum(x.^2)/del);
	dm=dy*sqrt(length(x)/del);
	
	
	yr = m*x+b;
	yre= dm*x+db;
	yrb=[yr-ratiosigma*yre,yr+ratiosigma*yre];
	
	out = [];
	keep = [];
	for i = 1:length(x)
	if (y(i)<yrb(i,1))|(y(i)>yrb(i,2))
		out(end+1) = i;
	else 
		keep(end+1) = i;
	end
	end
	
	
	
	
	ox = x(out);
	oy = y(out);
	outnames = names(out);
	writenames(outnames,'outliers.txt')
	
	x = x(keep);
	y = y(keep);
	names = names(keep);
	
	
	yr = m*x+b;
	yre= dm*x+db;
	yrb=[yr-ratiosigma*yre,yr+ratiosigma*yre];
	
	if length(ox) ==0
		ox = 0;
		oy = 0;

	end
end


function names = readnames
names = {};
fid = fopen('.tot.name');
i =1;
to = fgetl(fid);
while 1
    if ~ischar(to)
        break
    end
    names(i) = {to};
    to = fgetl(fid);
    i = i+1;
end
fclose(fid);

function writenames(names,fname)
fid = fopen(fname,'w');
for i = 1:length(names)
    fprintf(fid,'%s\n',char(strrep(names(i),'log','')));

end
fclose(fid);

