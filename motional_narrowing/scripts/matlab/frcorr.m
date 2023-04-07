function results = frcorr(fname,t_bin,logf,mode,mt)
% photon-photon auto-correlation
% t- binsize for correlation in s
% logf - to find bleach time
% Written per J. Chem. Phys. 2008, 128, 214101, for auto correlation statistical test
% also per J. Phys. Chem. B 2008, 112, 13962, for cross correlation statistical test
%
% modes:
%   1 - donor auto, region 2
%   2 - donor auto, region 3
%   3 - acceptor auto, region 2
%   4 - acceptor auto, region 3
%   5 - donor auto, region 1 
%   6 - acceptor auto, region 1
%   7 - cross, region 1
%   8 - recalculate error bars for new N_effective


[d_raw,a_raw] = load_data(fname,t_bin);


logf=load(logf);

switch mode
	case 1
		i_start = find(d_raw > logf(1,2));
		i_start = i_start(1);
		
		i_stop = find(d_raw < logf(1,1));
		i_stop = i_stop(end);	
		
        if i_start < i_stop
            I = d_raw(i_start:i_stop);
            [I,x] = hist(I,(I(end)-I(1))/t_bin);
        else
            I=[];
        end
        
	case 2
		i_start = find(d_raw > logf(1,1));
		i_start = i_start(1);	
	
		i_stop = length(d_raw);
        
        if i_start < i_stop
            I = d_raw(i_start:i_stop);
            [I,x] = hist(I,(I(end)-I(1))/t_bin);
        else
            I=[];
        end
            
	case 3
		i_start = find(a_raw > logf(1,2));
		i_start = i_start(1);	
	
		i_stop = find(a_raw < logf(1,1));
		i_stop = i_stop(end);	
		
        if i_start < i_stop
            I = a_raw(i_start:i_stop);
            [I,x] = hist(I,(I(end)-I(1))/t_bin);
        else
            I=[];
        end
        
	case 4
		i_start = find(a_raw > logf(1,1));
		i_start = i_start(1);	
	
		i_stop = length(a_raw);
		
        if i_start < i_stop
            I = a_raw(i_start:i_stop);
            [I,x] = hist(I,(I(end)-I(1))/t_bin);
        else
            I=[];
        end
	
	case 5
		i_start = 1;

		i_stop = find(d_raw <= logf(1,2));
		
		i_stop = i_stop(end);
        
        if i_start < i_stop
            I = d_raw(i_start:i_stop);
            [I,x] = hist(I,(I(end)-I(1))/t_bin);
        else
            I=[];
        end
        
	case 6
		i_start = 1;
		i_stop = find(a_raw <= logf(1,2));
		
		i_stop = i_stop(end);
        
		if i_start < i_stop
            I = a_raw(i_start:i_stop);
            [I,x] = hist(I,(I(end)-I(1))/t_bin);
        else
            I=[];
        end
		
	case 7
		i_start = 1;
		i_stop = find(d_raw <= logf(1,2));
		i_stop = i_stop(end);
		
        if i_start < i_stop
            Id = d_raw(i_start:i_stop);
            [Id,xd] = hist(Id,(Id(end)-Id(1))/t_bin);
        else
            Id=[];
        end
		
		i_stop = find(a_raw <= logf(1,2));
		i_stop = i_stop(end);
		
        if i_start < i_stop
            Ia = a_raw(i_start:i_stop);
            [Ia,xa] = hist(Ia,xd);
        else
            Ia=[];
        end
        
	case 8
		i_start = 1;
		i_stop = find(d_raw <= logf(1,2));
		i_stop = i_stop(end);
		
        if i_start < i_stop
            Id = d_raw(i_start:i_stop);
            [Id,xd] = hist(Id,(Id(end)-Id(1))/t_bin);
        else
            Id=[];
        end
		
		i_stop = find(a_raw <= logf(1,2));
		i_stop = i_stop(end);
		
        if i_start < i_stop
            Ia = a_raw(i_start:i_stop);
            [Ia,xa] = hist(Ia,xd);
        else
            Ia=[];
        end
        
end

if mode < 7

	N = length(I);
	if N == 0
		results = 0;
		return
	end	
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
elseif mode == 7 

	N = length(Id);
    
    if N == 0 || length(Ia)==0
		results = 0;
		return
    end	
    
	dId = fft(Id);
	dIa = fft(Ia);
	dI = dId .* conj(dIa);
	dI = dI/length(dI);
	dI(1) = 0;                
	g = ifft(dI);  
	g(1) = [];                
	g = g(1:floor(N/2));     
	tau = [1:length(g)] .* t_bin;

	ge = (N-1)./N.^2*(mean(Id.^2) - mean(Id).^2).*(mean(Ia.^2) - mean(Ia).^2);
	
	ge = ones(size(g))*sqrt(ge);

	results = [tau',g',ge'];
else
	
	N = length(Id);
    if N == 0 || length(Ia)==0
		results = 0;
		return
    end	 
	Neff = N/mt;
	results = sqrt((Neff-1)/Neff.^2*(mean(Id.^2) - mean(Id).^2).*(mean(Ia.^2) - mean(Ia).^2));	
end

function [d_raw,a_raw] = load_data(fname,t_bin)



I = photonload(fname);
dI = I(find(I(:,1)),1);
if size(I,2) == 1
	aI = dI; %hack for single channel data
else
	aI = I(find(I(:,2)),2);
end

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
	aI(turnovers(1)+1:length(aI))=aI(turnovers(1)+1:length(aI))+2^32*12.5e-9;
	turnovers(1)=[];
end

a_raw = aI;