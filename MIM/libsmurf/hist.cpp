#include "hist.h"

hist::hist(int ndim, int nh, int ns) {
    length=0; //This is the number of data points starting
    lavail=50; //The amount that all vectors will be incremented to add space for new data
    dim=ndim; //The dimensionality either 1 for x or 2 for x and x'
    scale=ns;  //Log scale or Linear scale..

    numhist=nh;    //The number of bins that will be used
    n=NumHist(lavail);
    //We should remove this later when adding 2D maps...

    hist_alloc_space();
    fill_hx();
    clear();
    //For testing purposes...
    //test();
}

void hist::test() {
    //This method is for debugging the code...
    //First thing to test is the matrix multiplication...
    gsl_matrix* one = gsl_matrix_calloc(3,2);
    gsl_matrix_set(one,0,0,23);
    gsl_matrix_set(one,1,0,2);
    gsl_matrix_set(one,2,0,-7);
    gsl_matrix_set(one,0,1,7);
    gsl_matrix_set(one,1,1,-9);
    gsl_matrix_set(one,2,1,34);
    gsl_matrix* two = gsl_matrix_calloc(2,3);
    gsl_matrix_set(two,0,0,6);
    gsl_matrix_set(two,0,1,-3);
    gsl_matrix_set(two,0,2,8);
    gsl_matrix_set(two,1,0,5);
    gsl_matrix_set(two,1,1,-2);
    gsl_matrix_set(two,1,2,7);
    //gsl_matrix* result;
    //    MatrixMul(one,CblasNoTrans,two,CblasNoTrans,&result);
    //    PrintMatrix(result);    
    /* 
       expected Result from above is:
       173  -83  233
       -33   12   -47
       128   -47  182
       and the function gives this result...
    */
    //MatrixMul(one,CblasTrans,two,CblasTrans,&result);
    //PrintMatrix(result);
    //The above lines should work but don't but perhaps we would never want to do something like take 
    //the transpose of both matrices...
  //Now we will test the function eig to make sure it is still intact....
  gsl_matrix* g = gsl_matrix_calloc(3,3);
  gsl_matrix* M = gsl_matrix_calloc(3,3);
  //Took some of the output from test_MEM to test this eig function...
  gsl_matrix_set(g,0,0,1.0);
  gsl_matrix_set(g,0,1,0.6372);
  gsl_matrix_set(g,0,2,0.4984);
  gsl_matrix_set(g,1,0,0.6372);
  gsl_matrix_set(g,1,1,1.0);
  gsl_matrix_set(g,1,2,-0.2311);
  gsl_matrix_set(g,2,0,0.4984);
  gsl_matrix_set(g,2,1,-0.2311);
  gsl_matrix_set(g,2,2,1.0);
  //Now for M...
  gsl_matrix_set(M,0,0,0.6271);
  gsl_matrix_set(M,0,1,0.3992);
  gsl_matrix_set(M,0,2,0.5263);
  gsl_matrix_set(M,1,0,0.3992);
  gsl_matrix_set(M,1,1,0.5208);
  gsl_matrix_set(M,1,2,0.0440);
  gsl_matrix_set(M,2,0,0.5263);
  gsl_matrix_set(M,2,1,0.0440);
  gsl_matrix_set(M,2,2,0.8041);
  /*
     W = {-1.000  -0.4515  -0.2781
           0.6981  0.7900   1.0000
           0.6094  1.0000  -0.1993}
   
     d = {9.0832       0       0
               0  0.9854       0
               0       0  2.1665}
  */
  //Here we have to use a function from the LAPACK++ Library to simultaneously diagonalize g and M.
  gsl_matrix* W;
  gsl_matrix* d;
  //[W,d] = eig(g,M)
  eig(&W,&d,g,M);
  //End of LAPACK++ Stuff...
}

hist::~hist() {
    hist_free_space();
}

void hist::hist_alloc_space(){

    rnd = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rnd,time(NULL));

    value=new gsl_vector*[dim];
    dev=new gsl_vector*[dim];
    hx=new gsl_vector*[dim];

    for(int i=0;i<dim;i++){
        value[i]=gsl_vector_alloc(lavail);
        dev[i]=gsl_vector_alloc(lavail);
        hx[i]=gsl_vector_calloc((int)gsl_pow_int(n,dim));
    }

    //This is false until the method bootstrap is called...
    //Needed so that if bootstrap is never called but DoMEM is called
    //sigma will get set properly...
    bootInit = false;

    weight=gsl_vector_alloc(lavail);
    cov=gsl_vector_alloc(lavail);
    hy=gsl_vector_calloc((int)gsl_pow_int(n,dim));
    thy=gsl_vector_calloc((int)gsl_pow_int(n,dim));    
    hdc=gsl_vector_calloc((int)gsl_pow_int(n,dim));
    hyt=gsl_vector_calloc((int)gsl_pow_int(n,dim));
    sigma=gsl_vector_calloc((int)gsl_pow_int(n,dim));
    //Need this so that when DoMEM is not called there is still a matrix
    //That will be freed in hist_free_space..Other wise you will get a seg fault....
    Response = gsl_matrix_calloc(1,1);

    minx=gsl_vector_calloc(dim);
    maxx=gsl_vector_calloc(dim);
    binsize=gsl_vector_calloc(dim);
    meanAlpha=gsl_vector_calloc(dim);
}

void hist::hist_free_space(){

    gsl_rng_free(rnd);
    for(int i=0;i<dim;i++){
        gsl_vector_free(dev[i]);
        gsl_vector_free(value[i]);
        gsl_vector_free(hx[i]);
    }
    delete[] hx;
    delete[] value;
    delete[] dev;

    gsl_vector_free(weight);
    gsl_vector_free(cov);
    gsl_vector_free(hy);
    gsl_vector_free(thy);
    gsl_vector_free(hdc);
    gsl_vector_free(hyt);

    //Free vectors used for Skillings mem deconvolution method...
    gsl_matrix_free(Response);

    gsl_vector_free(minx);
    gsl_vector_free(maxx);
    gsl_vector_free(binsize);
    gsl_vector_free(meanAlpha);
    gsl_vector_free(sigma);
}

void hist::hist_alloc_space_aug(){
    hx=new gsl_vector*[dim];

    for(int i=0;i<dim;i++){
        hx[i]=gsl_vector_calloc((int)gsl_pow_int(n,dim));
    }

    hy=gsl_vector_calloc((int)gsl_pow_int(n,dim));
    thy=gsl_vector_calloc((int)gsl_pow_int(n,dim));
    hdc=gsl_vector_calloc((int)gsl_pow_int(n,dim));
    hyt=gsl_vector_calloc((int)gsl_pow_int(n,dim));
    sigma=gsl_vector_calloc((int)gsl_pow_int(n,dim));

    minx=gsl_vector_calloc(dim);
    maxx=gsl_vector_calloc(dim);
    binsize=gsl_vector_calloc(dim);
    meanAlpha=gsl_vector_calloc(dim);
}

void hist::hist_free_space_aug(){

    for(int i=0;i<dim;i++){
        gsl_vector_free(hx[i]);
    }
    delete[] hx;
    gsl_vector_free(hy);
    gsl_vector_free(thy);
    gsl_vector_free(hdc);
    gsl_vector_free(hyt);

    gsl_vector_free(sigma);
    gsl_vector_free(meanAlpha);
    gsl_vector_free(minx);
    gsl_vector_free(maxx);
    gsl_vector_free(binsize);
}

void hist::fill_hx(){
    if(length==0) return;
    //We can remove the end point determination here because MaxEnt_iter takes care of that..
    double x,x2,sd,w;

    for(int i=0;i<dim;i++){
        x=0;
        x2=0;
        w=0;
        for(uint j=0;j<length;j++){
            x+=gsl_vector_get(weight,j)*gsl_vector_get(value[i],j);
            x2+=gsl_vector_get(weight,j)*gsl_pow_int(gsl_vector_get(value[i],j),2);
            w+=gsl_vector_get(weight,j);
        }
        x/=w;
        x2/=w;
        sd=pow(x2-gsl_pow_int(x,2),0.5);
        gsl_vector_set(minx,i,x-3*sd);
        gsl_vector_set(maxx,i,x+3*sd);
    }

    if(getenv("MINX0")) gsl_vector_set(minx,0,atof(getenv("MINX0")));
    if(getenv("MAXX0")) gsl_vector_set(maxx,0,atof(getenv("MAXX0")));
    if(dim>1){
		//cout << "setting MINX1 and MAXX1 from environment variables... " << gsl_vector_get(minx,1) << " "<< gsl_vector_get(maxx,1) << endl;
        if(getenv("MINX1")) gsl_vector_set(minx,1,atof(getenv("MINX1")));
        if(getenv("MAXX1")) gsl_vector_set(maxx,1,atof(getenv("MAXX1")));
    }

    gsl_vector_memcpy(binsize,maxx);
    gsl_vector_sub(binsize,minx);
    gsl_vector_scale(binsize,1.0/n);

	uint32_t ps[dim+1],ei;
    double tb,tl,ts;

	for(int i=0;i<dim+1;i++) {
		ps[i]=(uint32_t)gsl_pow_int(n,i);
	}

    for(int i=0;i<pow(n,dim);i++){
		for(int j=0;j<dim;j++) {
			ei=i%ps[j+1]/ps[j];
            tb=gsl_vector_get(binsize,j);
            tl=gsl_vector_get(maxx,j);
            ts=gsl_vector_get(minx,j);
            gsl_vector_set(hx[j],i,ts+(ei)*tb);
		}
    }
}

void hist::adddatum(s_x d){
    gsl_vector* sd=gsl_vector_alloc(dim);
    for(int i=0;i<dim;i++){
        gsl_vector_set(sd,i,pow(gsl_matrix_get(d.cS,i,i),0.5));
    }

    adddatum(d.tmax-d.tmin, d.c, sd, 0);

    gsl_vector_free(sd);
}

void hist::clear(){
    length=0;
}

void hist::adddatum(double w, gsl_vector* v, gsl_vector* d){
    adddatum(w,v,d,0);
}
 
void hist::adddatum(double w, gsl_vector* v, gsl_vector* d, double c){
    if(length+1>=lavail){
        augmentlavail();
    }

    gsl_vector_set(weight,length,w);
    gsl_vector_set(cov,length,c);

    for(int i=0;i<dim;i++){
        gsl_vector_set(value[i],length,gsl_vector_get(v,i));
        gsl_vector_set(dev[i],length,gsl_vector_get(d,i));
    }

    length++;
}

void hist::augmentlavail(){
    const int inc=50;

    gsl_vector* ph=gsl_vector_alloc(lavail);
    gsl_vector_view tv;

    gsl_vector_memcpy(ph,weight);
    gsl_vector_free(weight);
    weight=gsl_vector_calloc(lavail+inc);
    tv=gsl_vector_subvector(weight,0,lavail);
    gsl_vector_memcpy(&tv.vector,ph);

    gsl_vector_memcpy(ph,cov);
    gsl_vector_free(cov);
    cov=gsl_vector_calloc(lavail+inc);
    tv=gsl_vector_subvector(cov,0,lavail);
    gsl_vector_memcpy(&tv.vector,ph);

    for(int i=0;i<dim;i++){
        gsl_vector_memcpy(ph,value[i]);
        gsl_vector_free(value[i]);
        value[i]=gsl_vector_calloc(lavail+inc);
        tv=gsl_vector_subvector(value[i],0,lavail);
        gsl_vector_memcpy(&tv.vector,ph);

        gsl_vector_memcpy(ph,dev[i]);
        gsl_vector_free(dev[i]);
        dev[i]=gsl_vector_calloc(lavail+inc);
        tv=gsl_vector_subvector(dev[i],0,lavail);
        gsl_vector_memcpy(&tv.vector,ph);
    }

    lavail+=inc;
    n=NumHist(lavail);

    hist_free_space_aug();
    hist_alloc_space_aug();
    fill_hx();
    gsl_vector_free(ph);
}

//Could be useful for just printing data but isn't used by any functions right now...
//This will contain the same info as the .xc file from frag...
void hist::Print_data(ostream& dout){
    for(uint i=0;i<length;i++){
        dout << i << ' ';
        for(int j=0;j<dim;j++){
            dout << gsl_vector_get(value[j],i) << ' ' << gsl_vector_get(dev[j],i) << ' ';
        }
        dout << endl;
    }
}

void hist::MakeHist(double d){
    if(length==0){
        cerr << "MakeHist(): No data!" << endl;
        exit(1);
    }

    fill_hx();

    gausshist(false);

    gsl_vector_memcpy(hy,thy);
}

double hist::std(){
    double sx=0,sx2=0,n=0,tm;

    for(uint i=0;i<length;i++){
        tm=gsl_vector_get(weight,i)*gsl_vector_get(value[0],i);
        sx+=tm;
        sx2+=tm*gsl_vector_get(value[0],i);
        n+=gsl_vector_get(weight,i);
    }

    return pow(sx2/n-gsl_pow_int(sx/n,2),0.5);
}

double hist::mean(){
    double m=0,n=0;

    for(uint i=0;i<length;i++){
        m+=gsl_vector_get(value[0],i)*gsl_vector_get(weight,i);
        n+=gsl_vector_get(weight,i);
    }

    return m/n;
}

void hist::toProb(gsl_vector* ne){
    double total=0;
    for(int i=0;i<n;i++) {
        total+=gsl_vector_get(ne,i);
    }
    for(int i=0;i<n;i++) {
        gsl_vector_set(ne,i,gsl_vector_get(ne,i)/total);
    }
}

void hist::normalize(gsl_vector* ne){
    double dxdy=1;
    for(int k=0;k<dim;k++){
        dxdy*=gsl_vector_get(binsize,k);
    }

    double total=0;
    for(int i=0;i<gsl_pow_int(n,dim);i++){
        total+=gsl_vector_get(ne,i);
    }
    total*=dxdy;

    gsl_vector_scale(ne,1.0/total);
}

void hist::gausshist(bool bootstrapping){
    //bootstrapping=false by default as defined in hist.h...
    double inc;
    double tinc;
    double tot=0;
    for(uint i=0;i<length;i++) tot+=gsl_vector_get(weight,i);

    gsl_vector* tx=gsl_vector_alloc(dim);
    gsl_vector* tv=gsl_vector_alloc(dim);
    gsl_vector_set_all(thy,0.0);

    int bin;
    double min=gsl_vector_get(minx,0);
    double bs=gsl_vector_get(binsize,0);

    gsl_vector_set_all(meanAlpha,0);
    for(int i=0;i<dim;i++){
        for(uint j=0;j<length;j++){
            //this is setting td so that it is the sum of all the deviations squared...
            gsl_vector_set(meanAlpha,i,gsl_vector_get(meanAlpha,i)+gsl_pow_int(gsl_vector_get(dev[i],j),2));
        }
        //This is dividing the sum of the deviations squared by the number of data points
        //and taking the square root...
        gsl_vector_set(meanAlpha,i,pow(gsl_vector_get(meanAlpha,i)/length,0.5));
        //I want this meanAlpha vector so I can use it when constructing the response matrix...
    }

    if(dim==2) initHKernel2d();

    //This if statement is where the raw histogram is contructed considering all data to be gaussian distributed
    //with width alpha...
    if(!bootstrapping){
        for(int i=0;i<gsl_pow_int(n,dim);i++){
            //This is showing the percentage complete and has a spinner that goes along with it...
            spinner(i,(int)gsl_pow_int(n,dim),25);
            //inc is going to be the sum of the contributions as a result of adding alpha gaussian error 
            //to all of the values for all of the bins of the deconvoluted histrogram...
            inc=0;
            //Here we are simply setting tx so that it is equal to the position of the ith bin...
            for(int k=0;k<dim;k++){
                gsl_vector_set(tx,k,gsl_vector_get(hx[k],i));
            }
            //We will go through all the data and determine if each data point is a gaussian of width alpha
            //how much does it contribute to the current bin tx...
            for(uint j=0;j<length;j++){
                //Here we are simply setting tv equal to the jth data point extracted from value...
                for(int k=0;k<dim;k++){
                    gsl_vector_set(tv,k,gsl_vector_get(value[k],j));
                }
                //This is for distance and velocity where td is going to have both position and velocity information
                //We have to consider the std of velocity...
                if(dim==2){
                    bin=int((gsl_vector_get(tx,0)-min)/bs);
                    gsl_vector_set(meanAlpha,1,gsl_vector_get(kvp,bin));
                }

                //The HKernel method constructs a gaussian using varience value td and determines how much it fills
                //the tx bin...
                tinc=HKernel(tx,tv,meanAlpha,0);
                if(isnan(tinc)){
                    cout << __FILE__ << ' ' << __LINE__ << ": nan from HKernel." << endl;
                    cout << gsl_vector_get(tx,0) << '\t' << gsl_vector_get(tv,0) << '\t' << gsl_vector_get(meanAlpha,0) << endl;
                    continue;
                }
                //tinc is the amount the current data point fills the current bin (ith bin) and weight
                //is the size of the bin in time.  All the data is weighted by its size in time.
                inc+=tinc*gsl_vector_get(weight,j);
            }
            //now we have finished filling one bin will the effect of making all the data points guassians
            //we will add that to the final probability distribution thy before going to the next bin and filling it...
            gsl_vector_set(thy,i,gsl_vector_get(thy,i)+inc);
        }
    }
    // **********************
    //bootstrapping in one dimension
    
    //if guasshist is called by the method bootstrap then bootstrapping is positive and this if loop is entered instead
    //of the one above.  It is called reapeatedly the number of times you want to bootstrap...
    if(bootstrapping){
        for(int i=0;i<gsl_pow_int(n,dim);i++){
            //The amount of stuff that will go into one bin...
            inc=0;
            //set the current bin position tx to the ith position from the hx vector...
            for(int k=0;k<dim;k++){
                gsl_vector_set(tx,k,gsl_vector_get(hx[k],i));
            }
            
            //Here will select data points at random until we have reached length (with is the total number of data points)
            //These randomly selected data points will be used the same as above to construct gaussians and fill the
            //current bin by determining how much each gaussian contributes...
            //We will do this for every bin tx...
            for(uint j=0;j<length;j++){
                //Here we will be selecting random values...
                int randomNum=gsl_rng_uniform_int(rnd,length);
                for(int k=0;k<dim;k++){
                    gsl_vector_set(tv,k,gsl_vector_get(value[k],randomNum));
                }
                //for chaning variances in 2D...
                if(dim==2){
                    bin=int((gsl_vector_get(tx,0)-min)/bs);
                    gsl_vector_set(meanAlpha,1,gsl_vector_get(kvp,bin));
                }
                //Here once again HKernel is constructing the gaussian for the given piece of data and
                //then determining exactly how much it will contirbute to the current bin...
                tinc=HKernel(tx,tv,meanAlpha,0);
                if(isnan(tinc)){
                    cout << __FILE__ << ' ' << __LINE__ << ": nan from HKernel." << endl;
                    continue;
                }
                //Here we are adding the contribution of the current data point to the sum of all contributions...
                inc+=tinc*gsl_vector_get(weight,randomNum);
            }
            gsl_vector_set(thy,i,gsl_vector_get(thy,i)+inc);
        }
    }

// **********************

//    normalize(thy);
//Remember this line was causing artifact in 2D maps because it would changed one row of the matrix
//but it is needed for deconvolution so remember to rewrite toProb for 2D deconvolution...
    if (dim == 1) {
		toProb(thy);
	}
    gsl_vector_free(tx);
    gsl_vector_free(tv);
}

//This function was added back on 112606 for 2D map construction.....
void hist::initHKernel2d(){
	if(dim!=2) return;
	static int bin=0;
	if(bin) return;
	kvp=gsl_vector_calloc(n);
	gsl_vector* kvpn=gsl_vector_calloc(n);

	double tv=0;
	double min=gsl_vector_get(minx,0);
	double bs=gsl_vector_get(binsize,0);

	for(uint i=0;i<length;i++){
		tv=gsl_vector_get(value[0],i);
		bin=(int)floor((tv-min)/bs);
		if(gsl_vector_get(dev[1],i)>100) continue;
		gsl_vector_set(kvp,bin,gsl_vector_get(kvp,bin)+gsl_vector_get(dev[1],i));
		gsl_vector_set(kvpn,bin,gsl_vector_get(kvpn,bin)+1);
	}

	for(int i=0;i<n;i++){
		gsl_vector_set(kvp,i,gsl_vector_get(kvp,i)/gsl_vector_get(kvpn,i));
	}

	bool change=true;
	while(change){
		change=false;
		for(int i=1;i<n-1;i++){
			if(isnan(gsl_vector_get(kvp,i))){
				change=true;
				if(isnan(gsl_vector_get(kvp,i-1)) && isnan(gsl_vector_get(kvp,i+1)))
					continue;
				else if(isnan(gsl_vector_get(kvp,i-1)))
					gsl_vector_set(kvp,i,gsl_vector_get(kvp,i+1));
				else if(isnan(gsl_vector_get(kvp,i+1)))
					gsl_vector_set(kvp,i,gsl_vector_get(kvp,i-1));
				else
					gsl_vector_set(kvp,i,gsl_vector_get(kvp,i+1)/2+gsl_vector_get(kvp,i-1)/2);
			}
		}
	}		

	gsl_vector_set(kvp,0,gsl_vector_get(kvp,1));
	gsl_vector_set(kvp,n-1,gsl_vector_get(kvp,n-2));

	gsl_vector_free(kvpn);
}

//This function constructs a gaussian of varience d centered at v and sees how much it contributes to 
//the bin x...
double hist::HKernel(gsl_vector* x,gsl_vector* v, gsl_vector* d,double c){
    double z=0;
    double norm;
    static int cachemade=0;
    int cachesize=4000;
    double cachemax=10;
    double cached=cachemax/cachesize;

    if(!cachemade){
cout << "Generating exp cache... ";
        cache=new double[cachesize];
        for(int i=0;i<cachesize;i++){
            cache[i]=gsl_sf_exp(-i*cached);
        }
        cachemade=true;
cout << "done." << endl << flush;
    }

    if(dim == 1) {
        norm=gsl_vector_get(d,0)*M_SQRT2*M_SQRTPI;
        z=gsl_pow_2((gsl_vector_get(x,0)-gsl_vector_get(v,0))/gsl_vector_get(d,0))/2;
        if(z>cachemax) return 0;
        return cache[(int)(z/cached)]/norm;
    } else if(dim == 2) {
        norm=2*M_PI*gsl_vector_get(d,0)*gsl_vector_get(d,1);
        z=gsl_pow_2((gsl_vector_get(x,0)-gsl_vector_get(v,0))/gsl_vector_get(d,0));
        z+=gsl_pow_2((gsl_vector_get(x,1)-gsl_vector_get(v,1))/gsl_vector_get(d,1));
        if(z>cachemax/2) return 0;
        return cache[(int)(z/2/cached)]/norm;
    } else{
        cerr << "Dimensionality too high!" << endl;
        exit(1);
    }

    return 0;
}


void hist::Bootstrap(int nb) {
    //This method uses bootstrapping inside repeated calls to gausshist in order to 
    //determine the error in the histogram for each bin which is then put into the vector sigma...
    //If the number of bootstrap cycles is less than 2 exit...
    if(nb<2) return;

    bootInit = true;

    gsl_matrix* hdcM = gsl_matrix_alloc(nb,n);
    gsl_matrix* hyM = gsl_matrix_alloc(nb,n);
    gsl_matrix* hytM = gsl_matrix_alloc(nb,n);
    gsl_vector* Histmean = gsl_vector_calloc(n);
    gsl_vector* pmeanSq = gsl_vector_calloc(n);

    //I guess I am going to put all the histograms into a matrix and do the std stuff after the for loop...
    for(int i=0;i<nb;i++){
        //spinner(i,nb,25);
        //Every time this gets called thy is reset to the current histogram determined by bootstrapping...
        gausshist(true);
        gsl_vector_memcpy(hy,thy);
        DoMEM(2);
        for (int p=0;p<n;p++) {
                gsl_matrix_set(hdcM,i,p,gsl_vector_get(hdc,p));
                gsl_matrix_set(hyM,i,p,gsl_vector_get(hy,p));
                gsl_matrix_set(hytM,i,p,gsl_vector_get(hyt,p));
        }
    }
    PrintMatrix(hdcM,"hdcAll.out");
    PrintMatrix(hyM,"hyAll.out");
    PrintMatrix(hytM,"hytAll.out");
    for (int i=0;i<n;i++) {
        for (int p=0;p<nb;p++) {
            gsl_vector_set(Histmean,i,gsl_vector_get(Histmean,i)+gsl_matrix_get(hdcM,p,i));
        }
    }
    //now lets calculate sigma for every bin...
    //first step find the mean value for every bin...
    for (int p=0;p<n;p++) {
        gsl_vector_set(Histmean,p,gsl_vector_get(Histmean,p)/nb);
    }
    //lets set the deconvoluted output to be the mean of all of the bootstrap cycles...
    for (int p=0;p<n;p++) {
        gsl_vector_set(hdc,p,gsl_vector_get(Histmean,p));
    }
    //Now for the sum of the individual points minus the mean squared...
    double summertime; // = stupid programming problems....
    for (int i=0;i<n;i++) {
        for (int p=0;p<nb;p++) {
            summertime = gsl_vector_get(pmeanSq,i)+gsl_pow_2(gsl_matrix_get(hdcM,p,i)-gsl_vector_get(Histmean,i));
            gsl_vector_set(pmeanSq,i,summertime);
        }
    }
    //Now scale substracting one so that sigma is not underestimated...
    //Once again stupid gsl function wasn't doing what it was supposed to do....
    double scaleFactor = nb-1;
    for (int p=0;p<n;p++) {
        gsl_vector_set(pmeanSq,p,gsl_vector_get(pmeanSq,p)/scaleFactor);
    }
    //Now set the global variable sigma...
    for (int p=0;p<n;p++) {
        gsl_vector_set(sigma,p,sqrt(gsl_vector_get(pmeanSq,p)));
    }
    
    gsl_vector_free(pmeanSq);
    gsl_vector_free(Histmean);
    gsl_matrix_free(hdcM);
    gsl_matrix_free(hyM);
    gsl_matrix_free(hytM);           
    Print_MEM();
}
//The following function until DoMEM are left over from Lucas's code...
//I am going to leave them because they could be useful utility functions for debugging...
void hist::setNumHist(int num){
    numhist=num;
}

void hist::setHistX(gsl_vector* x,int ddd){
    gsl_vector_memcpy(hx[ddd],x);

    gsl_vector_set(minx,ddd,gsl_vector_get(hx[ddd],0));
    gsl_vector_set(maxx,ddd,gsl_vector_get(hx[ddd],(int)gsl_pow_int(n,dim)-1));
    gsl_vector_set(binsize,ddd,(gsl_vector_get(maxx,ddd)-gsl_vector_get(minx,ddd))/n);
}

void hist::setHistXAveDev(double ddd, int dim){
    gsl_vector_set_all(dev[dim],ddd);
    if(length==0) length++;
}

void hist::setHistRaw(gsl_vector* y){
    gsl_vector_memcpy(hy,y);
}

void hist::getHistRaw(gsl_vector* y){
    gsl_vector_memcpy(y,hy);
}

void hist::getHistDC(gsl_vector* y){
    gsl_vector_memcpy(y,hdc);
}

void hist::getHistX(gsl_vector* x, int dim){
    gsl_vector_memcpy(x,hx[dim]);
}

double hist::getHistX(int ind, int dim){
    return gsl_vector_get(hx[dim],ind);
}

void hist::DoMEM(int flag) {
    if ((flag > 2) || (flag < 1)) {
        cout << "flag must be set to 1 or 2" << endl;
        exit(1);
    }
    Print_MEM();
    //This should go before converting to probability because it depends on probability density...
    buildResponse();
    //change hy to probability instead of probability density...
    double hySum=0;
    for (int p=0;p<n;p++) {
        hySum += gsl_vector_get(hy,p);
    }
    gsl_vector_scale(hy,1/hySum);
    gsl_vector_scale(sigma,1/hySum);
    for (int w=0;w<n;w++) {
        if (gsl_vector_get(hy,w)<2.2*gsl_pow_int(10,-16)) {
            gsl_vector_set(hy,w,2.2*gsl_pow_int(10,-16));
        }
        if (gsl_vector_get(sigma,w)<2.2*gsl_pow_int(10,-16)) {
            gsl_vector_set(sigma,w,2.2*gsl_pow_int(10,-16));
        }
    }
    //need to be after hy is changed to probability because hdc is set equal to hy and should be in probability...
    init_hdc();

    for (int d=0;d<n;d++) {
        gsl_vector_set(sigma,d,1);
    }

    double success;
    //Initialize MaxEnt_iterate...
    success = MaxEnt_iterate(0);
    if (success) {
        cout << "MaxEnt Initialization succeeded..." << endl;
    } else {
        cout << "MaxEnt Initialization failed..." << endl;
        exit(1);
    }
    cout << "******************** Starting MaxEnt Iteration ********************" << endl;
    
    cout << scientific << endl;

    //iteration...
    int iter=1;
    double oldTestParam = 0;
    double old_chisq = 0;
    int converged = 0;
    double chisq = 0;
    while(!converged) {
        chisq = MaxEnt_iterate(flag);
        cout << iter << '\t' << "TEST = " << testparam << "\t" << "Chisq = " << chisq << endl;
        //test convergence conditions...
        if (flag==1) {
            if(fabs(testparam-oldTestParam)/testparam < 0.0001) {
                cout << "converged because TEST < 0.0001...\n" << endl;
                converged = 1;
            }
        } else if (flag==2) {
            if (chisq < 0.00001) {
                cout << "converged because chisq < 0.00001...\n" << endl;
                converged = 1;
            }
            if (fabs(chisq-old_chisq)/chisq < 0.0001) {
                   cout << "converged because Delta chisq < 0.0001..." << endl;
                   converged = 1;
            }
        }
	if ( iter > MAX_MAXENT_ITER ) {
	  cout << "teration aborted, exceeding maximum interation of " << MAX_MAXENT_ITER << "...\n" << endl;
	  converged = 1;
	}
        //house keeping...
        old_chisq = chisq;
        oldTestParam = testparam;
        iter += 1;
        Print_MEM();
    }
    toProb(hy);
    toProb(hyt);
    toProb(hdc);
}
void hist::init_hdc(){
    //For now this sets hdc to start as the raw histogram generated simply with the data...
    //We could start with it flat...
    for(int i=0;i<gsl_pow_int(n,dim);i++)
        gsl_vector_set(hdc,i,gsl_vector_get(hy,i));
}
/*
double hist::X2(){
    double x2=0;
    for(int i=0;i<numBins;i++){
        double valers=gsl_pow_2((gsl_vector_get(hy,v1+i)-gsl_vector_get(hyt,v1+i))/gsl_pow_2(gsl_vector_get(sigma,v1+i)));
        x2 += valers;
    }
    
    return x2/gsl_pow_int(n,dim);
}

double hist::H(){
    double h=0;

    for(int i=0;i<gsl_pow_int(n,dim);i++){
        if(gsl_vector_get(hdc,i)==0) continue;
        h+=gsl_vector_get(hdc,i)*gsl_sf_log(fabs(gsl_vector_get(hdc,i)));
    }

    return lambdaM*h/gsl_pow_int(n,dim);
}
*/
double hist::MaxEnt_iterate(int flag) {
    //f0 is going to be the first possible deconvoluted test funciton...
    //This will be hdc...
    //R is the transformation matrix which will stay in this method for now...
    //D is the experimental data which is stored inside hy...
    //sigma is the errors in the histogram determined by bootstrapping...
    //flag is the only input from the matlab script that I brought over...
    //flag==0 means initialization.  This is default as defined in hist.h..
    //flag==1 means MaxEnt method...
    //flag==2 means chi-squared minimization...
    
    //v1 and v2 are global variables so that chisq will work...
    static int v1;
    static int v2;
    static gsl_matrix* R;
    static gsl_matrix* Rt;
    static gsl_vector* sigma2;
    static gsl_matrix* ddC;
    static double C_exp;
    static int numBins;   //Once we determine the number of bins used to deconvolve we must differentiate it from the global variable n...

    double chiSq = 0;

    //initialization which could be separated out in the future...
    if (flag == 0) {
        //we need to make sure that the ends of the area being deconvoluted have sigma
        //larger than the machine acurracy...
        v1=-1;
        v2=-1;
        double Sigmin,Sigmax,hymin,hymax; 
        int vv1 = -1; 
        int vv2 = -1;
        double epsilon = sqrt(2.2*gsl_pow_int(10,-16));
        for (int s=0;s<n;s++) {
            Sigmin = gsl_vector_get(sigma,s);
            Sigmax = gsl_vector_get(sigma,n-s-1);
            hymin = gsl_vector_get(hy,s);
            hymax = gsl_vector_get(hy,n-s-1);
            //    printf("left side: %f  right side: %f\n",(test1*pow(10,12)),(test2*pow(10,12)));
            //In matlab eps was 2.2 * 10e-16 but for now lets make it large        
             if ((Sigmin>epsilon) && v1==-1) {
                v1 = s;   
            }
            if ((hymin>epsilon) && vv1==-1) {
                vv1 = s;
            }
            if ((Sigmax>epsilon) && v2==-1) {
                v2 = n-s-1;
            }
            if ((hymax>epsilon) && vv2==-1) {
                vv2 = n-s-1;
            }
        }
        if (vv1 > v1) {
            v1 = vv1;
        }
        if (vv2 > v2) {
            v2 = vv2;
        }
        numBins = v2-v1+1;
        //printf("xmin: %d,xmax: %d\n",v1,v2);
        //printf("numBins: %d\n",numBins);
        sigma2=gsl_vector_calloc(numBins);
        for (int w=0;w<numBins;w++) {
            gsl_vector_set(sigma2,w,gsl_pow_2(gsl_vector_get(sigma,w+v1))); 
        }
        gsl_matrix_view R_view = gsl_matrix_submatrix(Response,v1,v1,numBins,numBins);
        R = gsl_matrix_calloc(numBins,numBins);
        gsl_matrix_memcpy(R,&R_view.matrix);
        Rt = gsl_matrix_calloc(numBins,numBins);
        gsl_matrix_transpose_memcpy(Rt,R);
        gsl_matrix* tmp_R = gsl_matrix_calloc(numBins,numBins);
        gsl_matrix_memcpy(tmp_R,R);
        gsl_matrix* tmp_Rt = gsl_matrix_calloc(numBins,numBins); 
        gsl_matrix_memcpy(tmp_Rt,Rt);
        for (int g=0;g<numBins;g++) {
            for (int e=0;e<numBins;e++) {
                gsl_matrix_set(tmp_Rt,g,e,gsl_matrix_get(tmp_Rt,g,e)*gsl_vector_get(sigma,v1+e)); 
                gsl_matrix_set(tmp_R,e,g,gsl_matrix_get(tmp_R,e,g)*gsl_vector_get(sigma,v1+e));
            }
        }
        gsl_matrix* tmp_RTR;
        MatrixMul(tmp_Rt,CblasNoTrans,tmp_R,CblasNoTrans,&tmp_RTR);
        gsl_matrix_scale(tmp_RTR,2.0);
        //Here ddC hasn't yet been initialized so we are simple setting the pointer
        //to point to the tmp_RTR matrix...
        ddC = tmp_RTR;
        C_exp= numBins + 3.29*sqrt((double)numBins);
        gsl_matrix_free(tmp_Rt);
        gsl_matrix_free(tmp_R);
        //We don't want to free tmp_RTR because ddC points to the same matrix and will be used below...
        return 1;
   } else {
       //Set up some constants...
       double l0 = sqrt(0.1); //distance constraint so that the maximum x^2<=10^2...
       double a_min = 0;
       double a_max = 1;
       double EigenVectorTol = 0.00001;
       //Set up microiteration subspace...
       gsl_matrix* f=gsl_matrix_calloc(numBins,1);
       gsl_matrix* D=gsl_matrix_calloc(numBins,1);
       //hdc is the deconvoluted test histogram
       //hy is the actual data constructed using gaussians...
       for (int q=0;q<numBins;q++) {
           gsl_matrix_set(f,q,0,gsl_vector_get(hdc,v1+q));
           gsl_matrix_set(D,q,0,gsl_vector_get(hy,v1+q));
       }
       gsl_matrix* F;
         MatrixMul(R,CblasNoTrans,f,CblasNoTrans,&F);
       double entsum = 0;
       double ctimes = 0;
       double ctimes2 = 0;
       for (int y=0;y<numBins;y++) {
           ctimes = gsl_matrix_get(f,y,0);
           ctimes2 = gsl_sf_log(ctimes);
           if (isnan(ctimes2)||(gsl_isinf(ctimes2) == -1)) {
               ctimes2 = -20;
           }
           entsum += ctimes*ctimes2;
       }
       double A = gsl_sf_exp(entsum);
       double C0 = 0;
       gsl_matrix* tmp_pdC = gsl_matrix_alloc(numBins,1);
       for (int p=0;p<numBins;p++) {
           C0+=gsl_pow_2(gsl_matrix_get(F,p,0)-gsl_matrix_get(D,p,0))/gsl_vector_get(sigma2,p);
           gsl_matrix_set(tmp_pdC,p,0,(gsl_matrix_get(F,p,0)-gsl_matrix_get(D,p,0))/gsl_vector_get(sigma2,p));
       }
       gsl_matrix* dS = gsl_matrix_alloc(numBins,1);
       gsl_matrix* dC;
       MatrixMul(Rt,CblasNoTrans,tmp_pdC,CblasNoTrans,&dC);
       gsl_matrix_scale(dC,2.0);
       double logA = gsl_sf_log(A);
       double dS_abs = 0;
       double dC_abs = 0;
       for(int p=0;p<numBins;p++) {
           double carr = gsl_sf_log(gsl_matrix_get(f,p,0));
           gsl_matrix_set(dS,p,0,logA-carr);
           dS_abs += gsl_matrix_get(f,p,0)*gsl_pow_2(gsl_matrix_get(dS,p,0));
           dC_abs += gsl_matrix_get(f,p,0)*gsl_pow_2(gsl_matrix_get(dC,p,0));
       }
       dS_abs = sqrt(dS_abs);
       dC_abs = sqrt(dC_abs);

       //Now we can construct the search vectors stored in a new matrix e...
       
       gsl_matrix* e = gsl_matrix_calloc(numBins,3);
       double cf = 0;
       for (int p=0;p<numBins;p++) {
           cf = gsl_matrix_get(f,p,0);
           gsl_matrix_set(e,p,0,cf*gsl_matrix_get(dS,p,0));
           gsl_matrix_set(e,p,1,cf*gsl_matrix_get(dC,p,0));
       }
       gsl_matrix* ddCe1;
       gsl_matrix_view eSub = gsl_matrix_submatrix(e,0,0,numBins,1); 
       MatrixMul(ddC,CblasNoTrans,&eSub.matrix,CblasNoTrans,&ddCe1);
       gsl_matrix* ddCe2;
       gsl_matrix_view eSub2 = gsl_matrix_submatrix(e,0,1,numBins,1);
       MatrixMul(ddC,CblasNoTrans,&eSub2.matrix,CblasNoTrans,&ddCe2);
       //Could do this with matrix operations but it seems slower...
       for (int p=0;p<numBins;p++) {
           gsl_matrix_set(e,p,2,((gsl_matrix_get(f,p,0)*gsl_matrix_get(ddCe1,p,0))/dS_abs)-((gsl_matrix_get(f,p,0)*gsl_matrix_get(ddCe2,p,0))/dC_abs));
       }
       //Rescale e matrix...
       double sume1 = 0;
       double sume2 = 0;
       double sume3 = 0;
       for (int p=0;p<numBins;p++) {
           sume1 += gsl_pow_2(gsl_matrix_get(e,p,0));
           sume2 += gsl_pow_2(gsl_matrix_get(e,p,1));
           sume3 += gsl_pow_2(gsl_matrix_get(e,p,2));
       }
       sume1 = sqrt(sume1);
       sume2 = sqrt(sume2);
       sume3 = sqrt(sume3);
       //Could use a scale operation but I would have to define more vectors/matrixs or use row views might be better...
       for (int p=0;p<numBins;p++) {
           gsl_matrix_set(e,p,0,gsl_matrix_get(e,p,0)/sume1);
           gsl_matrix_set(e,p,1,gsl_matrix_get(e,p,1)/sume2);
           gsl_matrix_set(e,p,2,gsl_matrix_get(e,p,2)/sume3);
       }
       gsl_matrix* g;
       MatrixMul(e,CblasTrans,e,CblasNoTrans,&g);
       //All of the stuff below seems to be required for the line...
       //M = transpose(e) * [ddC*e(:,1) ddC*e(:,2) ddC*e(:,3)]
       gsl_matrix *ddCe1s,*ddCe2s,*ddCe3s;
       gsl_matrix_view eS1 = gsl_matrix_submatrix(e,0,0,numBins,1);
       gsl_matrix_view eS2 = gsl_matrix_submatrix(e,0,1,numBins,1);
       gsl_matrix_view eS3 = gsl_matrix_submatrix(e,0,2,numBins,1);
       MatrixMul(ddC,CblasNoTrans,&eS1.matrix,CblasNoTrans,&ddCe1s);
       MatrixMul(ddC,CblasNoTrans,&eS2.matrix,CblasNoTrans,&ddCe2s);
       MatrixMul(ddC,CblasNoTrans,&eS3.matrix,CblasNoTrans,&ddCe3s);
       gsl_matrix* ddCe = gsl_matrix_alloc(numBins,3);
       for (int p=0;p<numBins;p++) {
           gsl_matrix_set(ddCe,p,0,gsl_matrix_get(ddCe1s,p,0));
           gsl_matrix_set(ddCe,p,1,gsl_matrix_get(ddCe2s,p,0));
           gsl_matrix_set(ddCe,p,2,gsl_matrix_get(ddCe3s,p,0));
       }
       gsl_matrix* M;
       MatrixMul(e,CblasTrans,ddCe,CblasNoTrans,&M);
       /*
       printf("g = \n");
       PrintMatrix(g);
       printf("M = \n");
       PrintMatrix(M);
       */
       //Here we have to use a function from the CLAPACK Library to solve the general eigen value problem.
       //A*x=lambda*B*x for the eigenvalues lambda and eigen vectors x.
       gsl_matrix* W;
       gsl_matrix* d;
       //matlab function is [W,d] = eig(g,M).  We wrote a c wrapper for the function dggev_ so that we can simply
       //feed in gsl_matrix* and get out the desired result...
       eig(&W,&d,g,M);
       //End of CLAPACK Stuff....
       /*
       gsl_matrix* gm = gsl_matrix_calloc(2,2);
       gsl_matrix* Mm = gsl_matrix_calloc(2,2);
       gsl_matrix_set(gm,0,0,-0.25);
       gsl_matrix_set(gm,1,0,-0.7);
       gsl_matrix_set(gm,0,1,0.3);
       gsl_matrix_set(gm,1,1,2);
       printf("gm = \n");
       PrintMatrix(gm);
       PrintMatrix(gm,"gm.txt");
       gsl_matrix_set(Mm,0,0,-1.5);
       gsl_matrix_set(Mm,1,0,0.35);
       gsl_matrix_set(Mm,0,1,2);
       gsl_matrix_set(Mm,1,1,0.4);
       //gsl_matrix_set(Mm,0,0,1);
       //gsl_matrix_set(Mm,1,0,1);
       //gsl_matrix_set(Mm,0,1,1);
       //gsl_matrix_set(Mm,1,1,1);
       printf("Mm = \n");
       PrintMatrix(Mm);
       PrintMatrix(Mm,"Mm.txt");
       eig(&W,&d,gm,Mm);
       printf("W = \n");
       PrintMatrix(W);
       printf("d = \n");    
       PrintMatrix(d);
       //Does W diagonlize g?
       gsl_matrix* wpg;
       MatrixMul(W,CblasTrans,g,CblasNoTrans,&wpg);
       gsl_matrix* wpgw;
       MatrixMul(wpg,CblasNoTrans,W,CblasNoTrans,&wpgw);
       printf("wpgw = \n");
       PrintMatrix(wpgw);
       //Does W diagonalize M?
       gsl_matrix* wpM;
       MatrixMul(W,CblasTrans,M,CblasNoTrans,&wpM);
       gsl_matrix* wpMw;
       MatrixMul(wpM,CblasNoTrans,W,CblasNoTrans,&wpMw);
       printf("wpMw = \n");
       PrintMatrix(wpMw);

 
       printf("resetting W and d for debugging purposes\n");
       gsl_matrix_set(W,0,0,0.2954);
       gsl_matrix_set(W,1,0,-1);
       gsl_matrix_set(W,2,0,-0.7646);
       gsl_matrix_set(W,0,1,-0.7273);
       gsl_matrix_set(W,1,1,1);
       gsl_matrix_set(W,2,1,0.6624);
       gsl_matrix_set(W,0,2,-1);
       gsl_matrix_set(W,1,2,0.9447);
       gsl_matrix_set(W,2,2,-0.8289);
       printf("W = \n");
       PrintMatrix(W);

       gsl_matrix_set(d,0,0,6.3693);
       gsl_matrix_set(d,1,1,1.3942);
       gsl_matrix_set(d,2,2,0.5699);
       printf("d = \n");
       PrintMatrix(d);
       */

       gsl_matrix* tmp_g;
       MatrixMul(g,CblasNoTrans,W,CblasNoTrans,&tmp_g);
       gsl_matrix* tmp_W;
       MatrixMul(W,CblasTrans,tmp_g,CblasNoTrans,&tmp_W);
       gsl_vector_view gde = gsl_matrix_diagonal(tmp_W);
       //The search direction will then be along those that have
       //significant eigen values in gde...
       //produe new orthogonal vectors...
       gsl_matrix* e_orthonormal;
       MatrixMul(e,CblasNoTrans,W,CblasNoTrans,&e_orthonormal);
       //at this point, e_orthonormal(:,i) and e_orthonormal(:,j) are
       //orthogonal, that is, sum(e_orthonormal(:,i).*e_orthonormal(:,j))=0
       //re-normalize and store the orthonormal vectors in e...
       double e_orthSq1 = 0;
       double e_orthSq2 = 0;
       double e_orthSq3 = 0;
       for (int p=0;p<numBins;p++) {
           e_orthSq1 += gsl_pow_2(gsl_matrix_get(e_orthonormal,p,0));
           e_orthSq2 += gsl_pow_2(gsl_matrix_get(e_orthonormal,p,1));
           e_orthSq3 += gsl_pow_2(gsl_matrix_get(e_orthonormal,p,2));
       }
       e_orthSq1 = sqrt(e_orthSq1);
       e_orthSq2 = sqrt(e_orthSq2);
       e_orthSq3 = sqrt(e_orthSq3);
       for (int p=0;p<numBins;p++) {
           gsl_matrix_set(e,p,0,gsl_matrix_get(e_orthonormal,p,0)/e_orthSq1);
           gsl_matrix_set(e,p,1,gsl_matrix_get(e_orthonormal,p,1)/e_orthSq2);
           gsl_matrix_set(e,p,2,gsl_matrix_get(e_orthonormal,p,2)/e_orthSq3);
       }
       //Store the eigen values of M in Me vector...
       //Now that e has changed we have to remake M the same way as above...
       gsl_matrix* tmp_M;
       //We have to deallocate and reallocate here because ddC's are no longer the
       //correct dimensions for the memcpy operation below...
       gsl_matrix_free(ddCe1s);
       gsl_matrix_free(ddCe2s);
       gsl_matrix_free(ddCe3s);
       eS1 = gsl_matrix_submatrix(e,0,0,numBins,1); 
       eS2 = gsl_matrix_submatrix(e,0,1,numBins,1);
       eS3 = gsl_matrix_submatrix(e,0,2,numBins,1);
       MatrixMul(ddC,CblasNoTrans,&eS1.matrix,CblasNoTrans,&ddCe1s);
       MatrixMul(ddC,CblasNoTrans,&eS2.matrix,CblasNoTrans,&ddCe2s);
       MatrixMul(ddC,CblasNoTrans,&eS3.matrix,CblasNoTrans,&ddCe3s);
       for (int p=0;p<numBins;p++) {
           gsl_matrix_set(ddCe,p,0,gsl_matrix_get(ddCe1s,p,0));
           gsl_matrix_set(ddCe,p,1,gsl_matrix_get(ddCe2s,p,0));
           gsl_matrix_set(ddCe,p,2,gsl_matrix_get(ddCe3s,p,0));
       }
       MatrixMul(e,CblasTrans,ddCe,CblasNoTrans,&tmp_M);
       gsl_vector_view Me = gsl_matrix_diagonal(tmp_M);
       //printf("Me(1,1) = %f\n",gsl_vector_get(&Me.vector,0));
       //printf("Me(2,2) = %f\n",gsl_vector_get(&Me.vector,1));
       //printf("Me(3,3) = %f\n",gsl_vector_get(&Me.vector,2));
       gsl_matrix *Su,*Cu;
       MatrixMul(e,CblasTrans,dS,CblasNoTrans,&Su);
       MatrixMul(e,CblasTrans,dC,CblasNoTrans,&Cu);
       gsl_vector* C_min = gsl_vector_calloc(3);
       gsl_vector* C_aim = gsl_vector_calloc(3);
       for (int p=0;p<3;p++) {
           gsl_vector_set(C_min,p,C0-(0.5*gsl_pow_2(gsl_matrix_get(Cu,p,0))/gsl_vector_get(&Me.vector,p)));
       }
       gsl_vector_scale(C_min,2/3);
       gsl_vector_add_constant(C_min,1/3);
       for (int p=0;p<3;p++) {
           gsl_vector_set(C_aim,p,GSL_MAX_DBL(gsl_vector_get(C_min,p),C_exp));
       }
       //start micro-iteration...
       double gde_max = gsl_vector_max(&Me.vector);
       gsl_vector* xu = gsl_vector_calloc(3);
       if (flag == 1) {
           gsl_vector* Cp = gsl_vector_calloc(3);
           gsl_vector* Cc = gsl_vector_calloc(3);
           double a1 = a_min;
           double a2 = a_max;
           for (int u=0;u<3;u++) {
               double p_chop_finished=0;
               //Only deal with significant eigen values...
               if ((gsl_vector_get(&gde.vector,u)/gde_max) > EigenVectorTol) {
                   double p=0; //distance penalty, Q = a S - C - p l^2...
                   double a=(gsl_matrix_get(Cu,u,0)+gsl_vector_get(&Me.vector,u)*l0+p*l0)/(gsl_matrix_get(Su,u,0)-l0); 
                   //don't confuse the number 10 with l0 here.  This is the initial guess for a...
                   if ( a <= a_min) {
                       a = a_min;
                   } else if (a >= a_max) {
                       a = a_max;
                   }
                   while (!p_chop_finished) {
                       double a_chop_finished = 0;
                       double a_chop_success = 0;
                       while (!a_chop_finished) {
                           gsl_vector_set(xu,u,(a*gsl_matrix_get(Su,u,0)*gsl_matrix_get(Cu,u,0))/(p+gsl_vector_get(&Me.vector,u)+a));
                           gsl_vector_set(Cp,u,C0+gsl_matrix_get(Cu,u,0)*gsl_vector_get(xu,u)+0.5*(p+gsl_vector_get(&Me.vector,u))*gsl_pow_2(gsl_vector_get(xu,u)));
                           gsl_vector_set(Cc,u,C0+gsl_matrix_get(Cu,u,0)*gsl_vector_get(xu,u)+0.5*gsl_vector_get(&Me.vector,u)*gsl_pow_2(gsl_vector_get(xu,u)));
                           //Could add print out for debugging purposes here
                           //algorithm as displayed in figure 1.3
                           if (gsl_vector_get(Cp,u) > C0) {
                               a2 = a;
                               a = 0.5*(a+a1);
                           } else if (gsl_vector_get(Cc,u) < gsl_vector_get(C_aim,u)) {
                               a1 = a;
                               a = 0.5*(a+a2);
                           } else if (gsl_pow_2(gsl_vector_get(xu,u)) > gsl_pow_2(l0)) {
                               a1 = a;
                               a = 0.5*(a+a2);
                           } else {
                               a2 = a;
                               a = 0.5*(a+a1);
                               a_chop_success = 1;
                               a_chop_finished = 1;
                               p_chop_finished = 1;
                           }
                           //Criteria for a chop finished...
                           if ((fabs(a1-a2)/GSL_MAX_DBL(a1,a2) < 0.001) || (fabs(a-a_max)/a_max <0.001) || (fabs(a-a_min)/a <0.001)) {
                               a_chop_finished = 1;
                               p_chop_finished = 1;
                           }
                       } //end of while a_chop_finished...   
                   }//end of whil p_chop_finished...
               } else {
                   //The u-th direction is insignificant...
                   gsl_vector_set(xu,u,0);
               }
           }
           gsl_vector_free(Cc);
           gsl_vector_free(Cp);
       } else if (flag==2) {
           gsl_vector_set(xu,0,gsl_matrix_get(Cu,0,0));
           gsl_vector_set(xu,1,gsl_matrix_get(Cu,1,0));
           gsl_vector_set(xu,2,gsl_matrix_get(Cu,2,0));
           gsl_vector_div(xu,&Me.vector);
           gsl_vector_scale(xu,-1.0);
       }
       //printf("xu1 = %f\n",gsl_vector_get(xu,0));
       //printf("xu2 = %f\n",gsl_vector_get(xu,1));
       //printf("xu3 = %f\n",gsl_vector_get(xu,2));
       //I am just going to do this with a for loop and change to matrix operations later...
       //f_new is the deconvoluted histogram...
       //F_new is the reconvoluted histogram...
       gsl_matrix* f_new = gsl_matrix_calloc(numBins,1); 
       gsl_matrix* eXxu;
       gsl_matrix* xuM = gsl_matrix_alloc(3,1);
       gsl_matrix_set(xuM,0,0,gsl_vector_get(xu,0));
       gsl_matrix_set(xuM,1,0,gsl_vector_get(xu,1));
       gsl_matrix_set(xuM,2,0,gsl_vector_get(xu,2));
       MatrixMul(e,CblasNoTrans,xuM,CblasNoTrans,&eXxu);
       //PrintMatrix(eXxu,"eXxu.txt");
       for (int p=0;p<numBins;p++) {
           gsl_matrix_set(f_new,p,0,gsl_vector_get(hdc,p+v1)+gsl_matrix_get(eXxu,p,0));
           if (gsl_matrix_get(f_new,p,0)<2.2*gsl_pow_int(10,-16)) {
               gsl_matrix_set(f_new,p,0,2.2*gsl_pow_int(10,-16));
           }
       }
       //PrintMatrix(f_new,"f_new.txt");
       //now we will construct F_new which is the reconvoluted histogram and will eventually
       //be fed into the global vector thy... 
       gsl_matrix* F_new;
       MatrixMul(R,CblasNoTrans,f_new,CblasNoTrans,&F_new);
       //Now I will feed F_new and f_new into thy and hdc respectively which will add back the ends...
       for (int p=0;p<numBins;p++) {
               gsl_vector_set(hdc,p+v1,gsl_matrix_get(f_new,p,0));
               gsl_vector_set(hyt,p+v1,gsl_matrix_get(F_new,p,0));
       }
       //Construct the chiSq parameter for the return value...
       chiSq = 0;
       for(int i=0;i<numBins;i++){
           chiSq+=gsl_pow_2((gsl_vector_get(hy,v1+i)-gsl_vector_get(hyt,v1+i))/gsl_pow_2(gsl_vector_get(sigma,v1+i)));
       }
       //Construct the TEST parameter proposed by Skilling-Bryan, Eq. 37 in Ref. 1. 
       double sumdS = 0;
       double sumdC = 0;
       for (int p=0;p<numBins;p++) {
           sumdS += gsl_pow_2(gsl_matrix_get(dS,p,0));
           sumdC += gsl_pow_2(gsl_matrix_get(dC,p,0));
       }
       sumdS = sqrt(sumdS);
       sumdC = sqrt(sumdC);
       //testparam is a global variable for TEST described above...
       testparam = 0;
       double sumTest = 0;
       for (int p=0;p<numBins;p++) {
           testparam += 0.5*gsl_pow_2(gsl_matrix_get(dS,p,0)/sumdS-gsl_matrix_get(dC,p,0)/sumdC);
           sumTest += gsl_matrix_get(dS,p,0)/sumdS-gsl_matrix_get(dC,p,0)/sumdC;
       }
       if(sumTest == 0) {
           testparam *= 0;
       } else if (sumTest < 0) {
           testparam *= -1;
       }
       
       //Need to deallocate all the matrices I have setup before exiting...
       //It is possible I could basically make them all static and just make sure I am 
       //resetting all there element during every iteration...
       gsl_matrix_free(e);
       gsl_matrix_free(tmp_pdC);
       gsl_matrix_free(dS);
       gsl_matrix_free(dC);
       gsl_matrix_free(f);
       gsl_matrix_free(D);
       gsl_matrix_free(F);
       gsl_matrix_free(ddCe1);
       gsl_matrix_free(ddCe2);
       gsl_matrix_free(g);
       gsl_matrix_free(M);
       gsl_matrix_free(ddCe1s);
       gsl_matrix_free(ddCe2s);
       gsl_matrix_free(ddCe3s);
       gsl_matrix_free(ddCe);
       gsl_matrix_free(d);
       gsl_matrix_free(W);
       gsl_matrix_free(tmp_g);
       gsl_matrix_free(tmp_W);
       gsl_matrix_free(e_orthonormal);
       gsl_matrix_free(tmp_M);
       gsl_matrix_free(Su);
       gsl_matrix_free(Cu);
       gsl_vector_free(C_min);
       gsl_vector_free(C_aim);
       gsl_vector_free(xu);
       gsl_matrix_free(f_new);
       gsl_matrix_free(eXxu);
       gsl_matrix_free(F_new);
    }
    return chiSq;
}

void hist::MatrixMul(gsl_matrix* A,CBLAS_TRANSPOSE trA,gsl_matrix* B,CBLAS_TRANSPOSE trB,gsl_matrix** S) {
   //Dimensions of S should be the number of rows of the first one by the number of columns of the second one.
   //They can't be multiplied if the number of columns of the first one does not match the number of rows of the second one...
    if ((trA == CblasNoTrans) && (trB == CblasNoTrans)) { 
       if (!((*A).size2 == (*B).size1)) {
           printf("Number of Columns of Matrix A does not match number of rows of Matrix B.\n");
           return;
       }
       *S = gsl_matrix_calloc((*A).size1,(*B).size2); 
       gsl_blas_dgemm( trA, trB, 1.0, A, B, 0.0, *S );
    } else if ((trA == CblasTrans) && (trB == CblasNoTrans)) {
        if (!((*A).size1 ==(*B).size1)) {
            printf("Number of Columns of Transpose of Matrix A does not match number of rows of Matrix B.\n");
            return;
        }
        *S = gsl_matrix_calloc((*A).size2,(*B).size2);
        gsl_blas_dgemm( trA, trB, 1.0, A, B, 0.0, *S );
    } else if ((trA == CblasNoTrans) && (trB == CblasTrans)) {
        if (!((*A).size2 ==(*B).size2)) {
            printf("Number of Columns of Matrix A does not match number of rows of the Transpose of Matrix B.\n");
            return;
        }
        *S = gsl_matrix_calloc((*A).size1,(*B).size1);
        gsl_blas_dgemm( trA, trB, 1.0, A, B, 0.0, *S );
    } else if ((trA == CblasTrans) && (trB == CblasTrans)) {
        if (!((*A).size1 == (*B).size2)) {
            printf("Number of Columns of the transpose of matrix A does not match the number of rows of transpose matrix B.\n");
            return;
        }
        *S = gsl_matrix_calloc((*A).size2,(*B).size1);
        gsl_blas_dgemm( trA, trB, 1.0, A, B, 0.0, *S );
    }
}

void hist::buildResponse() {
    //This method constructs the Reponse matrix using the mean std in the input data
    //Which should be alpha if using constant information binning and working with position only
    Response = gsl_matrix_calloc(n,n);

    double tx = 0;
    double tv = 0;
    double bs = gsl_vector_get(binsize,0);
    double xmi = gsl_vector_get(minx,0);
    double sd = gsl_vector_get(meanAlpha,0)*sqrt(2.0);

    for(int p=0;p<n;p++) {
        //tv is the location where the gaussian is centered... 
        tx = xmi+p*bs;
        for(int i=0;i<n;i++){
            //tx is the location of the bin we are filling... 
            tv = xmi+i*bs;
            double hkern = gsl_sf_exp(-(gsl_pow_2(tv-tx)/2)/gsl_pow_2(sd))/(sqrt(2.0*M_PI)/sd);
            gsl_matrix_set(Response,p,i,hkern);                    
        }
    }
    //    printf("dx: %f\n",dx);
    //    gsl_matrix_scale(Response,dx);
    //PrintMatrix(Response,"Response.txt");
}

void hist::BootPrint(vector<double>* hxv, vector< vector<double> >* hyv, vector< vector<double> >* hdcv, vector< vector<double> >* hytv) {
    // Not properly implemented to work with velocity calculations..... 
    double bins = gsl_vector_get(binsize,0);
    double tothdc = 0;
    double tothy = 0;
    double tothyt = 0;
    for (int i=0;i<n;i++) {
        tothdc += gsl_vector_get(hdc,i);
        tothy += gsl_vector_get(hy,i);
        tothyt += gsl_vector_get(hyt,i);
    }
    tothdc *= bins;
    tothy *= bins;
    tothyt *= bins;
    (*hxv).clear();
    vector<double> hyvT, sigmavT, hdcvT, hytvT;
    for(int i=0;i<gsl_pow_int(n,dim);i++){
        for(int j=0;j<dim;j++){
            (*hxv).push_back(gsl_vector_get(hx[j],i));
        }
        hyvT.push_back(gsl_vector_get(hy,i)/tothy);
        hdcvT.push_back(gsl_vector_get(hdc,i)/tothdc);
        hytvT.push_back(gsl_vector_get(hyt,i)/tothyt);
    }
    (*hyv).push_back(hyvT);
    (*hdcv).push_back(hdcvT);
    (*hytv).push_back(hytvT);
}

void hist::Print(const char* mem, bool bmem){
    //This function will output in probability density units
    double bins = 1;
	for(int k=0;k<dim;k++) {
		bins*=gsl_vector_get(binsize,k);
	}
    double tothdc = 0;
    double tothy = 0;
    double tothyt = 0;
    for (int i=0;i<gsl_pow_int(n,dim);i++) {
        tothdc += gsl_vector_get(hdc,i);
        tothy += gsl_vector_get(hy,i);
        tothyt += gsl_vector_get(hyt,i);
    }
    tothdc *= bins;
    tothy *= bins;
    tothyt *= bins;
    ofstream mout(mem);
    mout << "%Histogram generated by boots." << endl;
    mout << scientific;
    double curSig = 0;
    for(int i=0;i<gsl_pow_int(n,dim);i++){
        for(int j=0;j<dim;j++){
            mout << gsl_vector_get(hx[j],i) << '\t';
        }
        mout << gsl_vector_get(hy,i)/tothy << '\t';
        mout << gsl_vector_get(hdc,i)/tothdc << '\t';
        if (bootInit) {
            curSig = gsl_vector_get(sigma,i);
        }
        mout << curSig/tothdc << '\t';
        mout << gsl_vector_get(hyt,i)/tothyt << endl;
    }
}

void hist::Print_MEM(){
    ofstream mout("mem.out");
    mout << scientific;
    double curSig = 0;
    for(int i=0;i<gsl_pow_int(n,dim);i++){
        for(int j=0;j<dim;j++){
            mout << gsl_vector_get(hx[j],i) << '\t';
        }
        if (bootInit) {
            curSig = gsl_vector_get(sigma,i);
        }
        mout << gsl_vector_get(hy,i) << '\t' << gsl_vector_get(hdc,i) << '\t' << curSig << '\t' << gsl_vector_get(hyt,i) << '\t' << endl;
    }
}
