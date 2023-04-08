#include "frag.h"
#include <gsl/gsl_sf_gamma.h>

int main(int argc, char* argv[]){
        double confidence1 = 0.95;
        double confidence2 = 0.95;
	int v=1;
	int o=0;
	int n=1;
	double alpha=.1;
	bool check=false,constt=false,efficiency=false;
	string infile="x.fr";
	string outfile="x.xc";
	string logfile="x.log";
	string tfile="x.ph";

	while(1){
		o=getopt(argc, argv, "w:u:v:prthgi:o:a:b:l:ed:cn:x:t:");
		if(o==-1)
			break;

		switch(o){
	        case 'w':
		        confidence1 = atof(optarg);
			break;
		case 'u':
		        confidence2 = atof(optarg);
			break;
		case 'i':
			//specify input file (unused)
			infile=optarg;
			break;
		case 'e':
			//calculate efficiency instead of distance
			efficiency=true;
			break;
		case 'o':
			//specify output file
			outfile=optarg;
			break;
		case 'l':
			//specify log file
			logfile=optarg;
			break;
		case 'a':
			//specify error
			alpha=atof(optarg);
			break;
		case 'v':
			//specify verbosity
			v=atoi(optarg);
			break;
		case 'c':
			//check validity and exit without calculating distances
			check=true;
			break;
		case 'n':
			//specify dimensionality of calculation
			n=atoi(optarg);
			break;
		case 't':
			//constant time binning: alpha is in seconds
			constt=true;
			break;
		case 'x':
			//specify crosstalk
			SMURF_XT=smurf_fret_xt_parse(optarg);
			break;
		case 'h': //fall through
		default:
			usage();
			return 1;
			break;
		}
	}

	int nin=argc-optind;
	if(nin>1){
		cout << "Only one input file can be accepted." << endl;
		exit(1);
	}

	if(nin==1) infile=argv[optind];
	string bfile=infile;
	bfile.erase(infile.length()-2,2);

	if(nin==1 && outfile=="x.xc"){
		if(efficiency) outfile=bfile+"xce";
		else outfile=bfile+"xc";
	}

	if(nin==1 && logfile=="x.log")
		logfile=bfile+"log";

	if(nin==1 && tfile=="x.ph")
		tfile=bfile+"ph";

	if(v>0) cout << "Loading " << infile << "..." << endl;
	s_fretdata fragbait(infile,tfile);
	if(v>0) fragbait.print_stats();
	//set the confidence level for the CoxOakes test for the distribution of interphoton distances
	//all of the coxOakes tests are done inside .isvalid() so we have to set it here.
	fragbait.CONFIDENCE1 = confidence1;
	fragbait.CONFIDENCE2 = confidence2;
	
	bool isValid = fragbait.isvalid();
	fragbait.print_log(logfile);
	if(!isValid){
		if(v>0) cout << "Not FRET!-in frag" << endl;
		return 1;
	}

	ofstream lout(logfile.c_str(),ios::app);
	//uint32_t lt=fragbait.b.nt;
	s_x tx(n);

	tx=fragbait.box(0,fragbait.b.ta_bleach);
	if(v>0) cout << "Average x: " << tx.x << ' ' << char(177) << ' ' << tx.std << endl;
	lout << tx.x << ' ' << tx.std << endl;

	tx=fragbait.boxe(0,fragbait.b.ta_bleach);
	if(v>0) cout << "Average E: " << tx.x << ' ' << char(177) << ' ' << tx.std << endl;
	lout << tx.x << ' ' << tx.std << endl;
	if(check) return 0;


	if(constt){
		ofstream xout(outfile.c_str());
		for(double t=0;t<fragbait.b.ta_bleach;t+=alpha){
			tx=fragbait.box(t,t+alpha);
			if(!tx.error) tx.print(xout);
		}
		return 0;
	}

	if(efficiency) outfile=outfile.replace(outfile.length()-2,2,"ec");
	ofstream xout(outfile.c_str());

	if(alpha>0){
		if(v>0) cout << "Using constant-information binning." << endl;
		if(efficiency) tx=fragbait.miestep(0,alpha,n);
		else tx=fragbait.mipstep(0,alpha,n);
		for(double t1=tx.tmax;tx.tmax<fragbait.b.ta_bleach;t1=tx.tmax){
			if(!tx.error) tx.printp(xout);
			if(efficiency) tx=fragbait.mipstep(t1,alpha,n);
			else tx=fragbait.mipstep(t1,alpha,n);
		}
	} else {
        if(v>0) cout << "Using constant-time binning." << endl;
        alpha*=-1;
		for(double t1=0, t2=alpha;t2 < fragbait.b.ta_bleach;t1+=alpha,t2+=alpha){
			if(efficiency) tx=fragbait.boxe(t1,t2);
			else tx=fragbait.box(t1,t2);
			if(!tx.error) tx.printp(xout);
		}
	}

//	if(efficiency){
//		outfile=outfile.replace(outfile.length()-2,2,"ec");
//		ofstream xout(outfile.c_str());
//		tx=fragbait.miestep(0,alpha,n);
//		
	//for(double t1=tx.tmax;tx.tmax<fragbait.b.ta_bleach;t1=tx.tmax){
///			if(!tx.error) tx.printp(xout);
//			tx=fragbait.miestep(t1,alpha,n);
//		}
//	}

	return 0;

	//Find average over total file
/*	s_fr np=read_num(tin,lt-1); //BUG: should use last not-bleached photon

	s_x tot;
	tot.tmin=0;
	tot.tmax=np.time;
	tot.x=smurf_MLE2B(np.cpcd, np.cpca, b);
	tot.std=smurf_STD2B(tot.x, tot.tmax-tot.tmin, b);

	ofstream lout(logfile);
	lout << b.td_bleach << ' ' << b.ta_bleach << endl;
	lout << b.IdB << ' ' << b.IaB << endl;
	lout << b.bd << ' ' << b.ba << endl;

	cout << "Average x: " << tot.x << ' ' << char(177) << ' ' << tot.std << endl;
	if((b.IdB<0)||(b.IaB<0)||(b.bd<0)||(b.ba<0)||(isnan(b.IdB))||(isnan(b.IaB))||(isnan(b.ba))||(isnan(b.bd))){
		cout << "Bad..." << endl;
		exit(1);
	}
*/
	//Setup recursion
//	s_x test_x;
//	s_x enc_x;

//	int nnan=1;
//	double test_int[2];
//	double enc_int[2];
//	double times[lt];
//	times[0]=0;
//	uint32_t ltimes=1;

//	double ltest=tot.tmax;
//	uint32_t ntest;
//	double dt=.1;
//	double split=0;
//	double tlim=0.005;//5.0/b.IdB;
//	uint32_t nntimes=0;
//	dt=0.001;

//	while(times[ltimes]<tot.tmax){
//		times[ltimes]+=dt;
//		find_split(tin,0,tot.tmax,times,ltimes,b);
//		if(split>0){
//			times[ltimes]=split;
//			ltimes++;
//			cout << split << endl;
//		}
//	}

/*	while(dt>tlim){
		ltest*=0.8;
		ntest=(uint32_t)(tot.tmax/ltest+1);
		dt=(tot.tmax-ltest)/(2*ntest-1);

		//fragment
		for(int i=0;i<ntest;i++){
			test_int[0]=0+i*dt;
			test_int[1]=test_int[0]+ltest;
			find_enc(test_int,enc_int,times,ltimes);
			test_x=read_dist(tin,test_int[0],test_int[1],b);
			enc_x=read_dist(tin,enc_int[0],enc_int[1],b);

			if(isnan(test_x.x)) continue;

			if(test_x!=enc_x){
				times[ltimes]=test_int[0];
				ltimes++;
				times[ltimes]=test_int[1];
				ltimes++;
			}
		}

		if(ltimes<6) continue;

		sort_times(times,ltimes);

		//consolidate
//		for(int i=0;i<ltimes-3;i++){
//			enc_int[0]= times[i+0];
//			test_int[0]=times[i+1];
//			test_int[1]=times[i+2];
//			enc_int[1]= times[i+3];
//
//			test_x=read_dist(tin,test_int[0],test_int[1],b);
//			enc_x=read_dist(tin,enc_int[0],enc_int[1],b);
//
//			if(isnan(enc_x.x)){
//				cout << "Unexpected NaN!" << endl;
//			}
//
//			if((test_x!=enc_x)&&(isfinite(test_x.x))){
//				continue;
//			} else {
//				times[i+2]=times[i+1]/2+times[i+2]/2;
//				times[i+1]=0;
//				i++;
//			}
//		}

//		nntimes=ltimes;
//		for(int i=1;i<ltimes;i++){
//			if(times[i]==0){
//				i++;
//				space++;
//			}
//			times[i-space]=times[i];
//		}
	}
*/
//	sort_times(times,ltimes);
//	print_times(tin, xout, b, times, ltimes);

//	cout << "Done: " << ltimes-1 << " points." << endl;
}

void find_split(fstream& tin, double tmin, double tmax, double times[],uint32_t& ltimes, s_cal b){
	double dt=(tmax-tmin)/50;
	s_x test1;
	s_x test2;
	double split=0;
	double ttest=0;
	double ttmax=0;

	for(double t=tmin;t<tmax;t+=dt){
		test1=read_dist(tin,tmin,t,b);
		test2=read_dist(tin,t,tmax,b);

		if((test1.isinvalid())||(test2.isinvalid()))
			continue;

//		ttest=abs(test1.x-test2.x)/sqrt(test1.std*test1.std+test2.std*test2.std);
		if(ttest>ttmax){
			split=t;
			ttmax=ttest;
		}
	}

	test1=read_dist(tin,tmin,split,b);
	test2=read_dist(tin,split,tmax,b);

	if((test1.std<.03)||(test1!=test2))
		find_split(tin,tmin,split,times,ltimes,b);

	if((test1.std<.03)||(test2.std<.03)||(test1!=test2)){
		times[ltimes]=split;
		ltimes++;
	}

	if((test2.std<.03)||(test1!=test2))
		find_split(tin,split,tmax,times,ltimes,b);
}

//Given times t1 and t2, data file f, and calibration info b, find the position, std, etc.
s_x read_dist(fstream& f, double t1, double t2, s_cal b){
	f.seekg(0,ios::end);
	uint32_t N=f.tellg()/sizeof(s_fr);

	s_fr p1a, p1b, p2a, p2b, pn;

	uint32_t n1a=0;
	uint32_t n1b=N-1;
	uint32_t n2a=0;
	uint32_t n2b=N-1;
	uint32_t nn=(uint32_t)(PHI*n1b);

	//Find p1a and p1b
	p1a=read_num(f, n1a);
	p1b=read_num(f, n1b);
	pn=read_num(f, nn);

	while(n1b-n1a>1){
		if(pn.time>t1){
			n1b=nn;
			p1b=read_num(f,n1b);
		} else if(pn.time<t1){
			n1a=nn;
			p1a=read_num(f,n1a);
		} else{
			n1a=nn;
			p1a=pn;
			n1b=nn;
			p1b=pn;
			break;
		}

		nn=(uint32_t)((1-PHI)*n1a+PHI*n1b);
		pn=read_num(f,nn);
	}

	//Find p1a and p1b
	nn=(uint32_t)(PHI*n2b);
	p2a=read_num(f, n2a);
	p2b=read_num(f, n2b);
	pn=read_num(f, nn);

	while(n2b-n2a>1){
		if(pn.time>t2){
			n2b=nn;
			p2b=read_num(f,n2b);
		} else if(pn.time<=t2){
			n2a=nn;
			p2a=read_num(f,n2a);
		} else{
			n2a=nn;
			p2a=pn;
			n2b=nn;
			p2b=pn;
			break;
		}

		nn=(uint32_t)((1-PHI)*n2a+PHI*n2b);
		pn=read_num(f,nn);
	}

	//double T=t2-t1;
	double cd=p2a.cpcd-p1b.cpcd+1;
	double ca=p2a.cpca-p1b.cpca+1;
	s_x nx;
	nx.tmin=t1;
	nx.tmax=t2;
	nx.x=smurf_MLE2B(cd, ca, b);
	nx.std=smurf_STD2B(nx.x, nx.tmax-nx.tmin, b);

	return nx;
}

//??
void find_enc(double test[], double enc[], double times[], int n){
	enc[0]=0;
	enc[1]=0;
	for(int i=0;i<n;i++)
		if(times[i]>enc[1]) enc[1]=times[i];

	for(int i=0;i<n;i++){
		if((times[i]<=test[0])&&(times[i]>enc[0])){
			enc[0]=times[i];
		}
		if((times[i]>=test[1])&&(times[i]<enc[1])){
			enc[1]=times[i];
		}
	}
}

void sort_times(double t[], int n){
	double tt[n];

	double max=0;
	for(int i=0;i<n;i++)
		if(t[i]>max) max=t[i];

	for(int i=0;i<n;i++)
		tt[i]=max;
	tt[0]=0;

	for(int i=1;i<n-1;i++)
		for(int j=0;j<n;j++){
//cout << tt[i] << ' ' << tt[i-1] << ' ' << t[j] << endl;
			if((t[j]<tt[i])&&(t[j]>tt[i-1]))
				tt[i]=t[j];
		}

	for(int i=0;i<n;i++)
		t[i]=tt[i];
}

//Print the x/std data given by the data in tf, calibration data b, and consecutive interval timings
//times to the ostream of.
void print_times(fstream& tf, ostream& of, s_cal b, double times[], int n){
	s_x tx;
	for(int i=0;i<n-1;i++){
		tx=read_dist(tf,times[i],times[i+1],b);
		tx.print(of);
	}
}

//Apply the Maximum Information Method to a trajectory contained in the file
//open in f, with calibration information b and desired standard deviation a.
void mim(fstream &f, double times[], uint32_t& ltimes, double a, s_cal b){
	times[0]=0;

	for(int ltimes=1;times[ltimes-1]<b.ta_bleach;ltimes++){
		times[ltimes]=mistep(f,times[ltimes-1],a,b);
	}
}

//Given a start time t1, find the end time t2 given std a and data f and b.
double mistep(fstream &f, double t1, double a, s_cal b){
	uint32_t cnum=findtime(f,b.nt-1,t1);
	s_fr fp=read_num(f,cnum);
	cnum++;
	s_fr lp=read_num(f,cnum);
	s_fr lpm;
	double x,std,tt;

	while(lp.time<=b.ta_bleach){
		x=smurf_MLE2B(lp.cpcd-fp.cpcd,lp.cpca-fp.cpca,b);
		std=smurf_STD2B(x,lp.time-t1,b);

		if(std<a){
			cnum--;
			lpm=read_num(f,cnum);
			x=smurf_MLE2B(lpm.cpcd-fp.cpcd,lpm.cpca-fp.cpca,b);
			std=smurf_STD2B(x,lp.time-t1,b);

			if((std>a)||(isnan(x))||(isnan(std))){
				return lp.time;
			}

			tt=t1+pow(smurf_STD2B(x,1,b)/a,2); //Information Rate
			if(tt>lp.time) cout << "Error in frag.cpp:mistep()" << endl;
			if(tt>b.ta_bleach) return b.ta_bleach;
			return tt;
		}

		cnum++;
		lp=read_num(f,cnum);
	}
	return b.ta_bleach;
}

void usage(){
	cerr << "frag: fret analyzer" << endl;
	cerr << endl;
	cerr << " -a SIGMA:\t\tSpecify desired accuracy." << endl;
	cerr << " -b FILE: \t\tSpecify base file name, without extension." << endl;
	cerr << "          \t\tDefault extensions: .fr, .ph, .log, .xc" << endl;
	cerr << " -c:      \t\tCheck validity and exit immediately." << endl;
	cerr << " -h:      \t\tDisplay this message." << endl;
	cerr << " -v NUM:  \t\tSet verbosity." << endl;
	cerr << " -n NUM:  \t\tSet dimensionality" << endl << endl;
}

/*
uint32_t process_data(fstream& fin, fstream &tout){
	//Find total photon numbers
	uint32_t ld;
	uint32_t la;

	fin.seekg(8);
	fin.read((char*)&(ld),sizeof(uint32_t));
	fin.read((char*)&(la),sizeof(uint32_t));
	uint32_t lt=la+ld;

	uint32_t ds=1024;
	uint32_t as=1024+ld*sizeof(uint32_t);

	uint32_t dta=0,dtd=0;
	uint64_t tta=0,ttd=0;

	//Clock rate
	uint32_t cr;
	fin.seekg(24);
	fin.read((char*)&cr,sizeof(cr));
	if(cr==0) cr=80000000;

	fin.seekg(ds);
	fin.read((char*)&dtd,sizeof(uint32_t));
	fin.seekg(as);
	fin.read((char*)&dta,sizeof(uint32_t));

	s_fr np;
	uint32_t cpcd=0;
	uint32_t cpca=0;

	while((cpca<la)&&(cpcd<ld)){
		if(tta+dta<ttd+dtd){
			tta+=dta;
			np.type=ACCEPTOR;
			np.time=tta/double(cr);
			np.cpca++;
			cpca++;
			fin.seekg(as+cpca*sizeof(uint32_t));
			fin.read((char*)&dta,sizeof(uint32_t));
			if(np.time==0){
				np.cpca--;
				continue;
				cout << "Acceptor dt=0" << endl;
			}
		} else if(ttd+dtd<=tta+dta){
			ttd+=dtd;
			np.type=DONOR;
			np.time=ttd/double(cr);
			np.cpcd++;
			cpcd++;
			fin.seekg(ds+cpcd*sizeof(uint32_t));
			fin.read((char*)&dtd,sizeof(uint32_t));
			if(np.time==0){
				np.cpcd--;
				continue;
			}
		}
		tout.write((char*)&np,sizeof(np));
	}

	return np.cpca+np.cpcd;
}
*/
