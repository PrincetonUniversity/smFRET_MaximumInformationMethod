#include <xcorrlib.h>
void usage();
int main(int argc, char* argv[]){
	
	int o;
	double T=1;
	double dt=.1;
	int v=1;
	int dim=1;
	int N=50; // (HY) number of time correlation points in log10 space
	string ofile="corr.out";
	bool constant_time=true;
	bool gaussian=false;
	double alpha=.1;
	double X1=-1;
	double X2=-1;

	while(1){
		o=getopt(argc,argv,"hgi:a:t:T:x:o:n:N:Y:Z:cg");
		if(o==-1) {
			break;
		}

		switch(o){
		case 'c':
			constant_time=false;
			break;
		case 'g':
			gaussian=true;
			break;
		case 'a':
			alpha=atof(optarg);
			break;
		case 't':
			dt=atof(optarg);
			break;
		case 'T':
			T=atof(optarg);
			break;
		case 'o':
			ofile=optarg;
			break;
		case 'N': // (HY-20060712)
			N=atoi(optarg);
			break;
		case 'n':
			dim=atoi(optarg);
			break;
		case 'x':
			SMURF_XT=smurf_fret_xt_parse(optarg);
			break;
		case 'Y':
			X1=atof(optarg);
			break;
		case 'Z':
			X2=atof(optarg);
			break;
		case 'h': // (HY-20060712) fall through
		default:
			usage();
			return 1;
			break;
		}
	}
	
	int nin=argc-optind;
	s_fretdata* farr[nin];

	if(nin==1 && ofile=="corr.out"){
		string prof=argv[optind];
		ofile=prof.replace(prof.length()-2,2,"xcor");
	}

	double totT=0;
	for(int i=0;i<nin;i++){
		if(v>0) cout << "Loading " << argv[i+optind] << endl;
		farr[i]=new s_fretdata(argv[i+optind]);

		if(!farr[i]->isvalid()){
			if(v>0) cout << "Not FRET!" << endl;
			continue;
		}
		totT+=farr[i]->b.ta_bleach;
	}

//	int N=(int)(T/dt)+1;
//	int N=5*log10(T/dt)+1;
        cout << "Correlation lentgh = " << N << endl;
	xcorr finland(farr, nin, 1);
	finland.set_times(T,dt,N);
	  // "set_times" assigns set global variables Tcorr = T, tcorr = t, and
	  //             Ncorr = 5*log10(T/t)+1 [xcorrlib.cpp]
	finland.set_alpha(alpha); // [xcorrlib.h] substitution macro, alpha = a;
	double tc[N];
	double xct[N];
	double xce[N];

	if((X1<0)||(X2<0)){
		if(constant_time){
			if(v>0) cout << "Using constant-time binning." << endl;
			finland.xcorr_calc_ct(tc,xct,xce); // (HY) this is the function that computes
		} else {
			if(v>0) cout << "Using constant-information binning." << endl;
			finland.xcorr_calc(tc, xct);
		}
	} else {
		if(gaussian){
			if(v>0) cout << "Using constant-time binning." << endl;
			finland.xcorr_calc_ct(tc, xct, X1, X2);
		} else {
			if(v>0) cout << "Using constant-time binning." << endl;
			finland.xcorr_calc_ct_pos(tc, xct, X1, X2); // (HY) position constrained correlation
		}
	}

	if(v>0) cout << endl << "1/e time: " << finland.xcorr_e(tc, xct, N) << endl;
	if(v>0) cout << endl << "integral: " << finland.xcorr_int(tc, xct, N) << endl;
	double te=2*finland.xcorr_int2(tc,xct,N)/totT;

	ofstream corout(ofile.c_str());
	for(int i=0;i<N;i++) corout << tc[i] << '\t' << xct[i] << '\t' << xce[i] << endl;
}

void usage() {
	cerr << "Usage: xcorr [OPTION]... FILES" << endl << endl
	     << "-h\tThis message" << endl
	     << "-x\tSpecify cross talk" << endl
	     << "-t\tSpecity constant-time increment" << endl
	     << "-T\tSpecify time duration to be correated" << endl
	     << "-N\tSpecify number of time point in the log-time averaged correlation function" << endl
	     << endl;
}
