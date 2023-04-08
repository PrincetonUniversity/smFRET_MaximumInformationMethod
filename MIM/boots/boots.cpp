#include "boots.h"

int main(int argc, char* argv[]) {    
    gsl_rng *rnd;
    int nboot=0;
    int ndim=1;
    double rtime=0;
    int nind;
    int v=1;
    int o=0;
    double alpha=0.1;
    double lambda=0;
    double tol=0.1;
    double ss=0.001;
    double nstd=1.0;
    double mean=0;
    int nbt=0;
    int NH=100;
    //This is the flag that determines whether only 'X2' minimization or 'X2+H' 
    //minimization is used in the method MaxEnt_interate(int flag) in hist.cpp...
    int lo=2;
    bool noopt=false,iqrc=false,stopafteriqr=false,mem=false,efficiency=false; 
    bool bmem=false,gotsfile=false,adjust_mean=false,bootstrapfiles=false;
    int zwanzig=0;
    string infile="x.fr",outfile="hist.out",logfile,tfile,iqrfile="iqr.out";
    string sfile;
    int niter=10;
    int Ndiag=25;
    int fileBootnum = 1;

    if (argc == 1) {
      // no arguments were supplied. print out instructions.
      usage();
      exit(1);
    }
    //parse options
    while(1){
        o=getopt(argc, argv, "hgi:v:prt:o:a:l:ep:s:n:E:qx:z:b:B:m::ZN:d:S:M");
        if(o==-1) {
            break;
	}
        switch(o){
        case 'q':
            iqrc=true;
            break;
        case 'N':
            NH=atoi(optarg);
            break;
        case 'M':
            adjust_mean=true;
            mean=atof(optarg);
            break;
        case 'z':
            zwanzig=1;
            if(optarg) nstd=atof(optarg);
            break;
        case 'Z':
            zwanzig=2;
            break;
        case 'd':
            zwanzig=3;
            if(optarg) Ndiag=atoi(optarg);
            break;
        case 'S':
            gotsfile=true;
            sfile=optarg;
            break;
        case 'i':
            infile=optarg;
            break;
        case 'e':
            efficiency=true;
            break;
        case 'E':
            bootstrapfiles=true;
            if(optarg) fileBootnum=atoi(optarg);
            break;
        case 'o':
            outfile=optarg;
            break;
        case 'p':
            tfile=optarg;
            break;
        case 'l':
            tol=atof(optarg);
            break;
        case 'a':
            alpha=atof(optarg);
            break;
        case 't':
            nbt=atoi(optarg);
            break;
        case 's':
            ss=atof(optarg);
            break;
        case 'v':
            v=atoi(optarg);
            break;
        case 'n':
            ndim=atoi(optarg);
            break;
        case 'b':
            nboot=atoi(optarg);
            break;
        case 'B':
            bmem=true;
            nboot=atoi(optarg);
            break;
        case 'x':
            SMURF_XT=smurf_fret_xt_parse(optarg);
            break;
        case 'm':
            mem=true;
//            if(nboot==0) nboot=10;
            if(optarg) lambda=atof(optarg);
            break;
         case 'h': //fall through
        default:
            usage();
            return 1;
            break;
        }
    }

    //Initialization of Random Number generator for file bootstrapping....
    rnd = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rnd,time(NULL));
    //Vectors that will hold data between different iterations of file bootstrapping...
    vector<double> hxv;
    vector< vector<double> > hyv, ey2v, hdcv, hytv;
    
    int nin=argc-optind;
    s_fretdata* farr[nin];
    s_fretdata* sorter[nin];

    //set output files
    if(nin==1 && outfile=="hist.out"){
        string prof=argv[optind];
        if(efficiency)
            outfile=prof.replace(prof.length()-2,2,"histe");
        else
            outfile=prof.replace(prof.length()-2,2,"hist");
    }

    if(outfile=="hist.out" && efficiency)
        outfile="histe.out";

    //load data files
    //nin: number of total data files
    for(int i=0;i<nin;i++) {
        if(v>0) cout << "Loading " << argv[i+optind] << endl;
        sorter[i]=new s_fretdata(argv[i+optind]);
    }

    if (nin == 1) {
      //If there is only one file we can't bootstrap the files...
      fileBootnum = 1;
      bootstrapfiles = false;
    }

    for(int q=1;q<=fileBootnum;q++) {  

    hist abby(ndim,NH,HIST_SCALE_LOG);
    s_x nd;

    if(mem){
        //lambda doesn't need to be set or separately optimized but I will leave this here for now...
        abby.set_lambda(lambda);
    }

    string names[nin];

    if (bootstrapfiles) {
        cout << "Bootstrapping " << nin << " files " << fileBootnum << " times.  Iteration Num: " << q << endl;
            //I stored all the files in sorter so that I could bootstrap them
        //into farr down here and farr is the vector that is actually 
        //fed into abbey which is a hist.cpp....
        for (int i=0;i<nin;i++) {
            int randN = gsl_rng_uniform_int(rnd,nin);
            farr[i] = sorter[randN];
            names[i] = argv[randN+optind];
        }
    } else {
        for (int i=0;i<nin;i++) {
            farr[i] = sorter[i];
        }
    }

    gsl_vector* btsd=gsl_vector_calloc((int)pow((double)abby.NumHist(abby.getlength()),(double)ndim));
    double dm;

//    ofstream xxpOut("xxp.out");

    //calculate histogram
    if(alpha>0){
        if(v>0) cout << "Filling histogram with constant-information binning: \n";
        for(int j=0;j<nin;j++){
            // if(v>0) cout << names[j] << endl;

            //find mean adjustment
            if(adjust_mean){
                dm=0;
                if(efficiency) nd=farr[j]->miestep(0,alpha,ndim);
                else nd=farr[j]->mipstep(0,alpha,ndim);
                for(double t1=nd.tmax;t1 < farr[j]->b.ta_bleach;t1=nd.tmax){
                    if(nd.tmax==0||isnan(gsl_vector_get(nd.c,0))) continue;
                    dm+=gsl_vector_get(nd.c,0)*(nd.tmax-nd.tmin);
                    if(efficiency) nd=farr[j]->miestep(t1,alpha,ndim);
                    else nd=farr[j]->mipstep(t1,alpha,ndim);
                }
                dm/=farr[j]->b.ta_bleach;
                dm=mean-dm;
            }

            if(efficiency) nd=farr[j]->miestep(0,alpha,ndim);
            else nd=farr[j]->mipstep(0,alpha,ndim);
//            xxpOut << nd.x << ' ' << gsl_vector_get(nd.c,0) << ' ' << gsl_vector_get(nd.c,1) << ' ' << sqrt(gsl_matrix_get(nd.cS,1,1)) <<  endl;
            for(double t1=nd.tmax;t1 < farr[j]->b.ta_bleach;t1=nd.tmax){
	    //while (nd.tmax < farr[j]->b.ta_bleach) {
                if(nd.tmax==0||isnan(gsl_vector_get(nd.c,0))) continue;
		  // (nd.c, 0) stores the MLE(x) starting with t1
                //if ( nd.tmax>0 ) {
                    if (adjust_mean) nd.adjust_mean(dm);
                    abby.adddatum(nd);
                    //if(efficiency) nd=farr[j]->miestep(t1,alpha,ndim);
                    //else nd=farr[j]->mipstep(t1,alpha,ndim);
		//}
                if(efficiency) nd=farr[j]->miestep(nd.tmax,alpha,ndim);
                else nd=farr[j]->mipstep(nd.tmax,alpha,ndim);
//                xxpOut << nd.x << ' ' << gsl_vector_get(nd.c,0) << ' ' << gsl_vector_get(nd.c,1) << ' ' << sqrt(gsl_matrix_get(nd.cS,1,1)) << endl;
            }
        }
    } else {
        if(v>0) cout << "Filling histogram with constant-time binning: ";
        //abbey.use_ind_var()  I think this is current not compatible so I will need to update this later for constant time binning case....
        alpha*=-1;
        for(int j=0;j<nin;j++){
            //cout << argv[j+optind] << endl;

            //find mean adjustment
            if(adjust_mean){
                dm=0;
                for(double t1=0, t2=alpha;t2 < farr[j]->b.ta_bleach;t1+=alpha,t2+=alpha){
                    if(efficiency) nd=farr[j]->boxe(t1,t2);
                    else nd=farr[j]->box(t1,t2);
                    if(nd.tmax==0||isnan(gsl_vector_get(nd.c,0))) continue;
                    dm+=gsl_vector_get(nd.c,0)*(nd.tmax-nd.tmin);
                }
                dm/=farr[j]->b.ta_bleach;
                dm=mean-dm;
            }
            
            for(double t1=0, t2=alpha;t2 < farr[j]->b.ta_bleach;t1+=alpha,t2+=alpha){
                if(efficiency) nd=farr[j]->boxe(t1,t2);
                else nd=farr[j]->box(t1,t2);
                if(nd.tmax==0||isnan(gsl_vector_get(nd.c,0))) continue;
                abby.adddatum(nd);
            }
        }
    }

    if(v>0) cout << "Constructing Histogram. " << endl;
    //build a histogram in abby using the data and variances....
    //MakeHist input does nothing at the moment...
    abby.MakeHist(-1.0);

    //bootstrapping is the only implementation of error estimation in hist.cpp currently...
    //I remove a lot of code that was here for correlation function error estimation and for loading histogram errors
    //I think it would be easy to add these back but it would take a little bit of rewriting and testing so I
    //will leave that for later...
        if(v>0) cout << "Bootstrapping." << endl;
        cout.flush();
    
    if(nboot) abby.Bootstrap(nboot);
    //MakeHist input does nothing at the moment...
    //abby.MakeHist(-1.0);
    //DoMEM is the equivalent of test_MEM.m it does the repeated Skilling construction and deconvolution stuff...
    if(mem && !nboot) abby.DoMEM(lo);

    if (bootstrapfiles) {
         abby.BootPrint(&hxv, &hyv, &hdcv, &hytv);
         outfile.erase(4);
         char buf[30];
         sprintf(buf, "%d", q);
         outfile.append(buf);
         outfile.append(".out");
         abby.Print(outfile.c_str(),bmem);
    } else {
         abby.Print(outfile.c_str(),bmem);
    }

    }

    if (bootstrapfiles) {
    int posSize=hxv.size()+1;
    //Now we will go through results of bootstrapping the files and calculate the mean and std dev 
    //for all of the various output distributions and then output a masterHist.out file.
    double hyM[posSize], hdcM[posSize], hytM[posSize];
    double hyMS[posSize], hdcMS[posSize], hytMS[posSize];
    double hdcSS[posSize];
    //initialize arrays to zero...
    for (uint32_t r=0;r<=hxv.size();r++) {
        hyM[r]=0;
        hdcM[r]=0;
        hytM[r]=0;
        hyMS[r]=0;
        hdcMS[r]=0;
        hytMS[r]=0;
        hdcSS[r]=0;
    }

    for (uint32_t c=0;c<hxv.size();c++) {
      for (uint32_t s=0;s<fileBootnum;s++) {
        hyM[c] += hyv[s][c];
        hdcM[c] += hdcv[s][c];
        hytM[c] += hytv[s][c];
      }
    }  
    
    //Now we get the mean using the sums from above...
        for (uint32_t n=0;n<hxv.size();n++) {
            hyM[n] = hyM[n]/fileBootnum;
            hytM[n] = hytM[n]/fileBootnum;
            hdcM[n] = hdcM[n]/fileBootnum;
        }

        for (uint32_t c=0;c<hxv.size();c++) {
          for (uint32_t s=0;s<fileBootnum;s++) {
              hdcSS[c] += pow((double)hdcv[s][c]-hdcM[c],2.0);
          }
        }

    //Now to calculate the standard deviation...
    for (uint32_t c=0;c<hxv.size();c++) {
        hdcMS[c] = (double)sqrt(hdcSS[c]/(fileBootnum-1));
    }

    //Now we will print out the masterHist.out file with 9 columns where 
    //every other column is the std of the previous column....
    outfile = "masterhist.out";
    ofstream mout(outfile.c_str());
    mout << "%Histogram generated by boots using bootstrapping of files." << endl;
    mout << scientific;
    for(int i=0;i<hxv.size();i++) {
      mout << hxv[i] << '\t' << hyM[i] << '\t' << hdcM[i] << '\t' << hdcMS[i] << '\t' << hytM[i] << '\t' << endl;
    }
    mout.close();
    }
    
    //clean up - was causing error so I got rid of it..
    for(int i=0;i<nin;i++){
      farr[i]->processed.close(); // close temporary files
      unlink(farr[i]->procfilename.c_str()); // delete temporary files
	  //
      //delete farr[i];
      //delete sorter[i];
    }

    return 0;
}

void usage(){
    cerr << endl << "Usage: boots [OPTION]... FILES" << endl
//    cerr << "\nboots\tCalculate distance histograms from FRET data" << endl
        << "  -o\tSpecify output file, default is base.hist or hist.out" << endl
        << "  -a\tSpecify smoothing parameter" << endl
        << "  -b[N]\tSpecify boot strap error estimation, where N is number of boot straps" << endl
	<< "  -v[n]\tVerbose mode when n > 0" << endl
        << "  -x\tSpecify crosstalk (-x0 or -x1)" << endl
        << "  -m\tDeconvolute histogram using entropy-regularized deconvolution" << endl
        << "  -E\tBootstrap on individual files" << endl
        << "  -h\tThis message" << endl
        << "Example: boots -x1 -a.1 -m -b50 Good_Data.fr" << endl << endl;
}

