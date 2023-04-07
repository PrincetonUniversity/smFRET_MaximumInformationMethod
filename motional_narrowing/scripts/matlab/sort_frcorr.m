function results = ftcorr(fname,t_bin,logf,mode);
% t- binsize for correlation in s
% logf - to find bleach time
%
% Written per J. Chem. Phys. 2008, 128, 214101
%
% modes:
%   1 - donor auto, region 2
%   2 - donor auto, region 3
%   3 - acceptor auto, region 2
%   4 - acceptor auto, region 3


[d_raw,a_raw] = load_data(fname,t_bin);


logf=load(logf);

switch mode
	case 1
		%   1 - donor auto, region 2
		i_start = find(d_raw > logf(1,2));
		i_start = i_start(1);
		
		i_stop = find(d_raw < logf(1,1));
		i_stop = i_stop(end);	
		
		I = d_raw(i_start:i_stop);
		
		
		[I,x] = hist(I,(I(end)-I(1))/t_bin);
	case 2
		%   2- donor auto, region 3
		i_start = find(d_raw > logf(1,1));
		i_start = i_start(1);	
	
		i_stop = length(d_raw);
		
		I = d_raw(i_start:i_stop);
		
		[I,x] = hist(I,(I(end)-I(1))/t_bin);
	case 3
		%   3- acceptor auto, region 2
		i_start = find(a_raw > logf(1,2));
		i_start = i_start(1);	
	
		i_stop = find(a_raw < logf(1,1));
		i_stop = i_stop(end);	
		
		I = a_raw(i_start:i_stop);
		

		[I,x] = hist(I,(I(end)-I(1))/t_bin);
	case 4
		%   4 - acceptor auto, region 3
		i_start = find(a_raw > logf(1,1));
		i_start = i_start(1);	
	
		i_stop = length(a_raw);
		
		I = a_raw(i_start:i_stop);
		 		
		[I,x] = hist(I,(I(end)-I(1))/t_bin);
end



N = length(I);
dI = fft(I);
dI = dI .* conj(dI);
dI(1) = 0;                 
g = ifft(dI)/length(dI);  
g(1) = [];                
g = g(1:floor(N/2));    
tau = [1:length(g)] .* t_bin;


c1 = 1/N^3;                
c2 = -4/N^3;               
c3 = (N-3)*(N+1)/N^3;      
c4 = -2*(N^2-2*N-6)/N^3;  
c5 = (N^2-2*N-6)/N^3;      
nt = (1:25);               

EX1 = mean(I);
EX2 = mean(I.^2);
EX3 = mean(I.^3);
EX4 = mean(I.^4);

A = c1.*EX4 + c2*EX3*EX1 + c3*EX2^2 + c4*EX2*EX1^2 + c5*EX1^4;

sCxx = sqrt(A);

ge = ones(size(tau))*sCxx;

results = [tau',g',ge'];

%==================================================
function [d_raw,a_raw] = load_data(fname,t_bin)


I = pload(fname);
dI = I(find(I(:,1)),1);
aI = I(find(I(:,2)),2);

dp=diff(dI);
turnovers=find(dp<-10);
while(length(turnovers>0))
	dI(turnovers(1)+1:length(dI))=dI(turnovers(1)+1:length(dI))+2^32*12.5e-9;
	turnovers(1)=[];
end

d_raw = dI;

ap=diff(aI);
turnovers=find(ap<-10);
while(length(turnovers>0))
	aI(turnovers(1)+1:length(aI))=a_bin(turnovers(1)+1:length(aI))+2^32*12.5e-9;
	turnovers(1)=[];
end

a_raw = aI;