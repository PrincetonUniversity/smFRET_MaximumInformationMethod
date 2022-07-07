#include "xcorrlib.h"
//#include "fretdata.h"
//#include "libsmurf.h"
#include <fstream>
using namespace std;

xcorr::xcorr(s_fretdata** f, int n, int d){
	farr=f;
	Ntraj=n;
	dim=d;
}

void xcorr::set_data(s_fretdata** d, int nin){
	farr=d;
	Ntraj=nin;
}

void xcorr::set_times(double T, double t, int N){
	// (HY-20060712) modify it so that N can be changed by user
	Tcorr=T;
	tcorr=t;
	Ncorr=N;
	//Ncorr=int(5*log10(T/t)+1);
}

void xcorr::set_times2(double T, double t, int Nt, int nx){
	// (HY-20060827) set constants for private scope
	Tcorr=T;
	tcorr=t;
	Ncorr=Nt;
	Nxcorr=nx;
}

/*double xcorr::xcorr_calc(double tc[], double xct[], double x1, double x2){

	bool tcc=false, xctc=false;
	if(!tc){
		tc=new double[Ncorr];
		tcc=true;
	}

	if(!xct){
		xct=new double[Ncorr];
		xctc=true;
	}

	int count=0;
	int corc=0;
	s_x xt;

	double xc[Ncorr];
	double Te[Ncorr];

	for(int i=0;i<Ncorr;i++) xct[i]=0;
    for(int i=0;i<Ncorr;i++) tc[i]=i*tcorr;
    for(int i=0;i<Ncorr;i++) Te[i]=0;

	for(int i=0;i<Ntraj;i++){
		if(!farr[i]->isvalid())
			continue;

		count=0;
		xt=farr[i]->mipstep(0,alpha,dim);
		while(xt.tmax< farr[i]->b.ta_bleach){
            xt=farr[i]->mipstep(xt.tmax,alpha,dim);
            count++;
        }

        if(xt.tmax<Tcorr){
            cout << "Skipping: too short" << endl;
            continue;
        }

		if((x1<0)||(x2<0)){
			s_x* xd=new s_x[count];
	
	        xd[0]=farr[i]->mipstep(0,alpha,dim);
	        for(int j=1;j<count;j++){
	            xd[j]=farr[i]->mipstep(xd[j-1].tmax,alpha,dim);
	        }
	
	        xcc(tc,xc,Ncorr,xd,xd,count,dim);
	        delete[] xd;
		} else {
			s_x xdt;
			s_x* xd1=new s_x[count];
			s_x* xd2=new s_x[count];
	
	        xdt=farr[i]->mipstep(0,alpha,dim);
			xd1[0]=xdt;
			gsl_vector_set(xd1[0].c,0,Gaussian(gsl_vector_get(xdt.c,0)-x1,alpha));
			xd2[0]=xdt;
			gsl_vector_set(xd2[0].c,0,Gaussian(gsl_vector_get(xdt.c,0)-x2,alpha));
	        for(int j=1;j<count;j++){
	            xdt=farr[i]->mipstep(xd1[j-1].tmax,alpha,dim);
				xd1[j]=xdt;
				gsl_vector_set(xd1[j].c,0,Gaussian(gsl_vector_get(xdt.c,0)-x1,alpha));
				xd2[j]=xdt;
				gsl_vector_set(xd2[j].c,0,Gaussian(gsl_vector_get(xdt.c,0)-x2,alpha));
	        }
	
	        xcc(tc,xc,Ncorr,xd1,xd2,count,dim);
	        delete[] xd1;
	        delete[] xd2;
		}

        for(int i=0;i<Ncorr;i++){ // if(tc[i]>xt.tmax) break;
            xct[i]+=xc[i];
            Te[i]+=xt.tmax-tc[i];
        }
        corc++;

    }

    for(int i=0;i<Ncorr;i++){
        xct[i]/=corc;
    }

	double rtime;
	if((x1<0)||(x2<0))
		rtime=xcorr_e(tc,xct,Ncorr);
	else{
		rtime=xcorr_int(tc,xct,Ncorr);
		ofstream corout("xctmp.xcor");
		for(int i=0;i<Ncorr;i++) corout << tc[i] << '\t' << xct[i] << endl;
	}

	if(tcc) delete[] tc;
	if(xctc) delete[] xct;

	return rtime;
}
*/

double xcorr::xcorr_calc(double tc[], double xct[], double x1, double x2){

	bool tcc=false, xctc=false;
	if(!tc){
		tc=new double[Ncorr];
		tcc=true;
	}

	if(!xct){
		xct=new double[Ncorr];
		xctc=true;
	}

	int count=0;
	int corc=0;
	double totT=0;
	s_x xt;

	double xc[Ncorr];
	double Te[Ncorr];

	for(int i=0;i<Ncorr;i++) xct[i]=0;
    for(int i=0;i<Ncorr;i++) tc[i]=i*tcorr;
    for(int i=0;i<Ncorr;i++) Te[i]=0;

	for(int i=0;i<Ntraj;i++){
		if(!farr[i]->isvalid())
			continue;

		count=0;
		xt=farr[i]->mipstep(0,alpha,dim);
		while(xt.tmax< farr[i]->b.ta_bleach){
            xt=farr[i]->mipstep(xt.tmax,alpha,dim);
            count++;
        }

        if(xt.tmax<Tcorr){
            cout << "Skipping: too short" << endl;
            continue;
        }

		s_x xdt;
		s_x* xd1=new s_x[count];
		s_x* xd2=new s_x[count];
		totT+=farr[i]->b.ta_bleach;

        xdt=farr[i]->mipstep(0,alpha,dim);
		xd1[0]=xdt;
		gsl_vector_set(xd1[0].c,0,Gaussian(gsl_vector_get(xdt.c,0)-x1,alpha));
		xd2[0]=xdt;
		gsl_vector_set(xd2[0].c,0,Gaussian(gsl_vector_get(xdt.c,0)-x2,alpha));
        for(int j=1;j<count;j++){
            xdt=farr[i]->mipstep(xd1[j-1].tmax,alpha,dim);
			xd1[j]=xdt;
			gsl_vector_set(xd1[j].c,0,Gaussian(gsl_vector_get(xdt.c,0)-x1,alpha));
			xd2[j]=xdt;
			gsl_vector_set(xd2[j].c,0,Gaussian(gsl_vector_get(xdt.c,0)-x2,alpha));
        }

        xcc(tc,xc,Ncorr,xd1,xd2,count,dim);
        delete[] xd1;
        delete[] xd2;

        for(int i=0;i<Ncorr;i++){ // if(tc[i]>xt.tmax) break;
            xct[i]+=xc[i];
            Te[i]+=xt.tmax-tc[i];
        }
        corc++;

    }

    for(int i=0;i<Ncorr;i++){
        xct[i]/=corc;
    }

	double rtime;
	rtime=xcorr_int(tc,xct,Ncorr);
	ofstream corout("xctmp.xcor");
	for(int i=0;i<Ncorr;i++) corout << tc[i] << '\t' << xct[i] << endl;

	if(tcc) delete[] tc;
	if(xctc) delete[] xct;

	return totT;
}

// This subroutine computes constant-information binning correlation 
double xcorr::xcorr_calc(double tc[], double xct[]){

	bool tcc=false, xctc=false;
	if(!tc){
		tc=new double[Ncorr];
		tcc=true;
	}

	if(!xct){
		xct=new double[Ncorr];
		xctc=true;
	}

	int count=0;
	int corc=0;
	double totT=0;
	s_x xt;

	double xc[Ncorr];
	double Te[Ncorr];

	for(int i=0;i<Ncorr;i++) xct[i]=0;
    for(int i=0;i<Ncorr;i++) tc[i]=i*tcorr;
    for(int i=0;i<Ncorr;i++) Te[i]=0;

	for(int i=0;i<Ntraj;i++){
		if(!farr[i]->isvalid())
			continue;

		count=0;
		xt=farr[i]->mipstep(0,alpha,dim);
		while(xt.tmax< farr[i]->b.ta_bleach){
            xt=farr[i]->mipstep(xt.tmax,alpha,dim);
            count++;
        }

        if(xt.tmax<Tcorr){
            cout << "Skipping: too short" << endl;
            continue;
        }

		s_x* xd=new s_x[count];
		totT+=farr[i]->b.ta_bleach;

        xd[0]=farr[i]->mipstep(0,alpha,dim);
        for(int j=1;j<count;j++){
            xd[j]=farr[i]->mipstep(xd[j-1].tmax,alpha,dim);
        }

        xcc(tc,xc,Ncorr,xd,xd,count,dim);
        delete[] xd;

        for(int i=0;i<Ncorr;i++){ // if(tc[i]>xt.tmax) break;
            xct[i]+=xc[i];
            Te[i]+=xt.tmax-tc[i];
        }
        corc++;

    }

    for(int i=0;i<Ncorr;i++){
        xct[i]/=corc;
    }

	double rtime;
	rtime=xcorr_e(tc,xct,Ncorr);

	if(tcc) delete[] tc;
	if(xctc) delete[] xct;

	return totT;
}

/*double xcorr::xcorr_calc_ct(double tc[], double xct[], double x1, double x2){

	bool tcc=false, xctc=false;
	if(!tc){
		tc=new double[Ncorr];
		tcc=true;
	}

	if(!xct){
		xct=new double[Ncorr];
		xctc=true;
	}

	int count=0;
	int corc=0;
	s_x xt;

	double xc[Ncorr];
	double Te[Ncorr];

	for(int i=0;i<Ncorr;i++) xct[i]=0;
    for(int i=0;i<Ncorr;i++) tc[i]=i*tcorr;
    for(int i=0;i<Ncorr;i++) Te[i]=0;

	for(int i=0;i<Ntraj;i++){
		if(!farr[i]->isvalid())
			continue;

        if(farr[i]->b.ta_bleach<Tcorr){
            cout << "Skipping: too short" << endl;
            continue;
        }

		s_x xdt;
		s_x* xd1=new s_x[(int)(farr[i]->b.ta_bleach/alpha)];
		s_x* xd2=new s_x[(int)(farr[i]->b.ta_bleach/alpha)];

		int j=0;
        for(double t1=0, t2=alpha;t2 < farr[i]->b.ta_bleach;t1+=alpha, t2+=alpha, j++){
            xdt=farr[i]->box(t1,t2);
			xd1[j]=xdt;
			gsl_vector_set(xd1[j].c,0,Gaussian(gsl_vector_get(xdt.c,0)-x1,alpha));
			xd2[j]=xdt;
			gsl_vector_set(xd2[j].c,0,Gaussian(gsl_vector_get(xdt.c,0)-x2,alpha));
        }

        ctc(tc,xc,Ncorr,xd1,xd2,count,dim);
        delete[] xd1;
        delete[] xd2;

        for(int i=0;i<Ncorr;i++){ // if(tc[i]>xt.tmax) break;
            xct[i]+=xc[i];
            Te[i]+=xt.tmax-tc[i];
        }
        corc++;

    }

    for(int i=0;i<Ncorr;i++){
        xct[i]/=corc;
    }

	double rtime;
	rtime=xcorr_int(tc,xct,Ncorr);
	ofstream corout("xctmp.xcor");
	for(int i=0;i<Ncorr;i++) corout << tc[i] << '\t' << xct[i] << endl;

	if(tcc) delete[] tc;
	if(xctc) delete[] xct;

	return rtime;
}
*/

/********** xcorr filename -x1 -t0.001 -T5 **********/
double xcorr::xcorr_calc_ct(double tc[], double xct[], double xce[]){
  // (HY-20060711) This is the subroutine that computes distance time 
  //               correlation function using xcorr -filename -t -T -x1

	bool tcc=false, xctc=false;
	if(!tc){
		tc=new double[Ncorr];
		tcc=true;
	}

	if(!xct){
		xct=new double[Ncorr];
		xctc=true;
	}

	int count=0;
	s_x xt;

	int Ncorri = (int) floor(Tcorr/tcorr);
	int corc[Ncorr];
	int j, k;
	double xc[Ncorri]; // (HY-20060712), was xc[Ncorr]
	double xc_error[Ncorri]; // (HY-20060717), add
	double xc_n[Ncorri]; // (HY-20060717), add
	double xctn[Ncorr];
	double xcte[Ncorr];
	double Te[Ncorri]; // (HY-20060712), was Te[Ncorr]
	double totT[Ncorri]; // (HY-20060712), was totT[Ncorr]
	double dlogtc; //(HY-20060712)
	double logt0=log10(tcorr), logtN=log10(Tcorr), t1, t2;

	for(int i=0;i<Ncorr;i++) corc[i]=0;
	for(int i=0;i<Ncorr;i++) xct[i]=0;
	// (HY-20060712)
	// The log-average code does not make sense. Rewrite.
	dlogtc = (logtN-logt0) / (double) Ncorr;
	t1 = tcorr;
	for(int i=0;i<Ncorr;i++) {
		t2 = pow(10,logt0+i*dlogtc);
		if ( t2 > t1 ) {
			tc[i] = t2;
		} else {
			tc[i] = t1;
		}
		t1 += tcorr;
		//tc[i]=pow(10,i/5.0)*tcorr;
	}
	//
        for(int i=0;i<Ncorr;i++) Te[i]=0;
        for(int i=0;i<Ncorr;i++) totT[i]=0;

	s_x* xd;
	for(int i=0;i<Ntraj;i++){
		// (HY) loop through trajectories
		count=(int)(farr[i]->b.ta_bleach/tcorr);
		xd=new s_x[count]; // xd is the distance as a function of time
		count=ct_fill(farr[i],xd,tcorr,Tcorr);
		if(count==-1) continue;
		if(count<Tcorr/tcorr){
			//Ncorri=int(5*log10(count)+1); (HY-20060712)
			//cout << "Too short... ncorri = " << Ncorri << endl;
			cout << "Too short... T_traj = " << count*tcorr << endl;
			continue;
		}
        	ctc(xc,xc_n,xc_error,count,xd,xd,count,dim); // (HY) costant-time corr
				// xc[1..count], constant-time correlation function
				// xcn[1..count], number of correlations in each element
				// xce[1..count], error as estimated by <xt^2x0^2>-<xtx0>^2
        	delete[] xd;
		// (HY) The following for loop does log averaging
		// (HY-20060712) Rewrite this section.
		k=0;
		xct[k] += xc[k];
		xcte[k] += xc_error[k]*xc_n[k];
		xctn[k] += xc_n[k];
		corc[k]++;
		t1 = 2*tcorr;
		k++;
		j=1;
		while (j<Ncorri && k < (Ncorr-1)) {
			if ( ( tc[k-1] < t1 ) && ( t1 <= tc[k+1] ) ) {
				xct[k] += xc[j];
				xctn[k] += xc_n[j];
				xcte[k] += xc_error[j]*xc_n[k];
				corc[k] ++;
			} else if (t1 > tc[k+1]) {
				k++;
				xct[k] += xc[j];
				xctn[k] += xc_n[j];
				xcte[k] += xc_error[j]*xc_n[k];
				corc[k]++;
			}
			t1 += tcorr;
			j++;
		}
		while ( (j<Ncorri) && (k==Ncorr-1) ) {
			xct[k]+=xc[j];
			xctn[k] += xc_n[j];
			xcte[k] += xc_error[j]*xc_n[k];
			corc[k]++;
			j++;
		}
    } // end of looping through trajectories
	for(int i=0;i<Ncorr;i++){
		xct[i]/=corc[i];
		xce[i]=sqrt(xcte[i]/xctn[i]/xctn[i]);
	}
	double rtime;
	rtime=xcorr_e(tc,xct,Ncorr);
	//double int2=xcorr_int2(tc,xct,Ncorr);
	// (HY-20060717) comment out incorrect error estimation
	//for(int i=0;i<Ncorr;i++){
	//	xce[i]=sqrt(2*int2/totT[i]);
	//}
	if(tcc) delete[] tc;
	if(xctc) delete[] xct;
	return 0;
}

double xcorr::xcorr_calc_ct_xx(double tc[], double *xxct, double *xxce, double x1, double x2){
  // (HY-20060827) This is the subroutine that computes distance-distance 
  //               time orrelation function using x2corr -filename -t -T -x1

    int count=0;
    s_x xt;

    int Ncorri = (int) floor(Tcorr/tcorr); // mumber of equal-time indices
    int corc[Ncorr];
    int j, k;
    double *xxc, *xxc_error, *xxc_n, *xxctn, *xxcte;
    double Te[Ncorri]; // (HY-20060712), was Te[Ncorr]
    double totT[Ncorri]; // (HY-20060712), was totT[Ncorr]
    double dlogtc; //(HY-20060712)
    double logt0=log10(tcorr), logtN=log10(Tcorr), t1, t2;

    // create work space for various variables
    xxc = (double *) calloc( Ncorri*Nxcorr*Nxcorr, sizeof(double) );
    xxc_error = (double *) calloc( Ncorri*Nxcorr*Nxcorr, sizeof(double) );
    xxc_n = (double *) calloc( Ncorri, sizeof(double) );
    xxctn = (double *) calloc( Ncorr, sizeof(double) );
    xxcte = (double *) calloc( Ncorr*Nxcorr*Nxcorr, sizeof(double) );

    for(int i=0;i<Ncorr;i++) corc[i]=0;
    // (HY-20060712)
    // The log-average code does not make sense. Rewrite.
    dlogtc = (logtN-logt0) / (double) Ncorr;
    t1 = tcorr;
    tc[0] = 0.0;
    for(int i=1;i<Ncorr;i++) {
        t2 = pow(10,logt0+i*dlogtc);
        if ( t2 > t1 ) {
            tc[i] = t2;
        } else {
            tc[i] = t1;
        }
        t1 += tcorr;
        //tc[i]=pow(10,i/5.0)*tcorr;
    }
    //
        for(int i=0;i<Ncorr;i++) Te[i]=0;
        for(int i=0;i<Ncorr;i++) totT[i]=0;

    s_x* xd;
    for(int i=0;i<Ntraj;i++){
        // (HY) loop through trajectories
        count=(int)(farr[i]->b.ta_bleach/tcorr);
        xd=new s_x[count]; // xd is the distance as a function of time
        count=ct_fill(farr[i],xd,tcorr,Tcorr);
        if(count==-1) continue;
        if(count<Tcorr/tcorr){
            //Ncorri=int(5*log10(count)+1); (HY-20060712)
            //cout << "Too short... ncorri = " << Ncorri << endl;
            cout << "Too short... T_traj = " << count*tcorr << endl;
            continue;
        }
        ctc_xx(xxc,xxc_n,xxc_error,count,xd,xd,count,dim,x1,x2); // (HY) costant-time corr
            // count, number of x measures (constant time) along the trajectory
            // xc[1..count], constant-time correlation function
            // xc_n[1..count], number of correlations in each element
            // xc_error[1..count], error as estimated by <xt^2x0^2>-<xtx0>^2
        delete[] xd;

        // (HY) The following for loop does log averaging
        k=0; // the k-th time delay on the log time scale
        for (int ii=0; ii<Nxcorr; ii++) // copy the k-th time map from xxc to xxct
            for (int jj=0; jj<Nxcorr; jj++)
                xxct[ii+Nxcorr*(jj+Nxcorr*k)] += xxc[ii+Nxcorr*(jj+Nxcorr*k)];
        //xcte[k] += xc_error[k]*xc_n[k]; // (HY-20060827) error function, deal with it later
        xxctn[k] += xxc_n[k]; // total number of x measures at the k-th time delay
        corc[k]++; // how many times the k-th delay has been accumulated, may vary from traj to traj
        t1 = tcorr;
        k++;
        j=1;
        while (j<Ncorri && k < (Ncorr-1)) {
            if ( ( tc[k-1] < t1 ) && ( t1 <= tc[k+1] ) ) {
                for (int ii=0; ii<Nxcorr; ii++) // copy the k-th time map from xxc to xxct
                    for (int jj=0; jj<Nxcorr; jj++)
                        xxct[ii+Nxcorr*(jj+Nxcorr*k)] += xxc[ii+Nxcorr*(jj+Nxcorr*j)];
                //xct[k] += xc[j];
                xxctn[k] += xxc_n[j];
                //xxcte[k] += xxc_error[j]*xxc_n[k];
                corc[k] ++;
            } else if (t1 > tc[k+1]) {
                k++;
                for (int ii=0; ii<Nxcorr; ii++) // copy the k-th time map from xxc to xxct
                    for (int jj=0; jj<Nxcorr; jj++)
                        xxct[ii+Nxcorr*(jj+Nxcorr*k)] += xxc[ii+Nxcorr*(jj+Nxcorr*j)];
                //xct[k] += xc[j];
                xxctn[k] += xxc_n[j];
                //xcte[k] += xc_error[j]*xc_n[k];
                corc[k]++;
            }
            t1 += tcorr;
            j++;
        }
        while ( (j<Ncorri) && (k==Ncorr-1) ) {
            for (int ii=0; ii<Nxcorr; ii++) // copy the k-th time map from xxc to xxct
                for (int jj=0; jj<Nxcorr; jj++)
                    xxct[ii+Nxcorr*(jj+Nxcorr*k)] += xxc[ii+Nxcorr*(jj+Nxcorr*j)];
            //xct[k]+=xc[j];
            xxctn[k] += xxc_n[j];
            //xcte[k] += xc_error[j]*xc_n[k];
            corc[k]++;
            j++;
        }

    } // end of looping through trajectories

    for(int i=0;i<Ncorr;i++){
        for (int ii=0; ii<Nxcorr; ii++) // copy the k-th time map from xxc to xxct
            for (int jj=0; jj<Nxcorr; jj++)
                xxct[ii+Nxcorr*(jj+Nxcorr*i)] /= corc[i];
        //xce[i]=sqrt(xcte[i]/xctn[i]/xctn[i]);
    }


    // free up workspace
    free(xxc);
    free(xxc_error);
    free(xxc_n);
    free(xxctn);
    free(xxcte);

    return 0;
}

// (HY) 20060815: constant-time binning with Gaussians centering at x1, x2
double xcorr::xcorr_calc_ct(double tc[], double xct[], double x1, double x2){
	bool tcc=false, xctc=false;

	if(!tc){
		tc=new double[Ncorr];
		tcc=true;
	}

	if(!xct){
		xct=new double[Ncorr];
		xctc=true;
	}

	int count=0;
	int corc=0;
	s_x xt;

	double xc[Ncorr];
	double Te[Ncorr];

	for(int i=0;i<Ncorr;i++) xct[i]=0;
        for(int i=0;i<Ncorr;i++) tc[i]=i*tcorr;
        for(int i=0;i<Ncorr;i++) Te[i]=0;

	for(int i=0;i<Ntraj;i++){
		if(!farr[i]->isvalid())
			continue;

        	if(farr[i]->b.ta_bleach<Tcorr){
            		cout << "Skipping: too short" << endl;
            		continue;
        	}	

		s_x xdt;
		count=(int)(farr[i]->b.ta_bleach/tcorr);
		s_x* xd1=new s_x[count];
		s_x* xd2=new s_x[count];

		int j=0;
        	for(double t1=0, t2=tcorr;t2 < farr[i]->b.ta_bleach;t1+=tcorr, t2+=tcorr, j++){
            		xdt=farr[i]->box(t1,t2);
			xd1[j]=xdt;
			gsl_vector_set(xd1[j].c,0,Gaussian(gsl_vector_get(xdt.c,0)-x1,alpha));
			//cout << Gaussian(gsl_vector_get(xdt.c,0)-x1,alpha) << endl;
			xd2[j]=xdt;
			gsl_vector_set(xd2[j].c,0,Gaussian(gsl_vector_get(xdt.c,0)-x2,alpha));
        	}

//        	ctc(tc,xc,Ncorr,xd1,xd2,count,dim); // (HY-20060717) this line has to be revised to match the new ctc() function
        	delete[] xd1;
        	delete[] xd2;

        	for(int i=0;i<Ncorr;i++){ // if(tc[i]>xt.tmax) break;
            		xct[i]+=xc[i];
            		Te[i]+=xt.tmax-tc[i];
        	}
        	corc++; // (HY) number of trajectories included
    	} // (HY) end of for (i,0,Ntraj)
/***************/

    	for(int i=0;i<Ncorr;i++){
        	xct[i]/=corc; // (HY) average over total number of trajectories
    	}

	double rtime;
	rtime=xcorr_e(tc,xct,Ncorr);

	if(tcc) delete[] tc;
	if(xctc) delete[] xct;

	return rtime;
}

double xcorr::xcorr_calc_ct_pos(double tc[], double xct[], double X1, double X2){
	// (HY-20060711)
	// X1 and X2, range that defines valid x[] such that X1 <= x_valid <= X2

	bool tcc=false, xctc=false;
	if(!tc){
		tc=new double[Ncorr];
		tcc=true;
	}

	if(!xct){
		xct=new double[Ncorr];
		xctc=true;
	}

	int count=0;
	int corc=0;
	s_x xt;

	double xc[Ncorr];
	double Te[Ncorr];

	for(int i=0;i<Ncorr;i++) xct[i]=0;
    for(int i=0;i<Ncorr;i++) tc[i]=i*tcorr;
    for(int i=0;i<Ncorr;i++) Te[i]=0;

	s_x* xd;
	for(int i=0;i<Ntraj;i++){ 
		// (HY) find out how many time points are there in a traj
		count=(int)(farr[i]->b.ta_bleach/tcorr);
		// (HY) create a new array to store them
		xd=new s_x[count];
		// (HY) convert photon counting data to x[t]
		count=ct_fill(farr[i],xd,tcorr,Tcorr);
		if(count==-1) continue;

        ctc_pos(tc,xc,Ncorr,xd,xd,count,dim,X1,X2);
	  // tc[], not used in ctc_pos()
	  // xc[], returns distance fluctuation correlation function
        delete[] xd;

        for(int k=0;k<Ncorr;k++){ // if(tc[i]>xt.tmax) break;
	    // (HY) was for int i=0, ... this should be a bug that may cause errors
            xct[k]+=xc[k];
            Te[k]+=xt.tmax-tc[k]; // ???? tc[] was never computed in ctc_pos
        }
        corc++; // (HY) number of trajectories in the total correlation function
    }

    for(int i=0;i<Ncorr;i++){
        xct[i]/=corc;
    }

	double rtime;
	rtime=xcorr_e(tc,xct,Ncorr);

	if(tcc) delete[] tc;
	if(xctc) delete[] xct;

	return rtime;
}

int xcorr::ct_fill(s_fretdata* fd, s_x* xd, double tcorr, double Tcorr){
	// (HY) constant-time correation function to generate x(t)
	if(!fd->isvalid())
		return -1;

	int count=(int)(fd->b.ta_bleach/tcorr); // (HY) number of time points

	int j=0;
        // (HY) xd[.] is the array that stores distances in x = R/R0
	// (HY) tcorr is the time increment in correlation functions
	for(double t1=0,t2=tcorr; t2<fd->b.ta_bleach;t1+=tcorr, t2+=tcorr, j++){
		xd[j]=fd->box(t1,t2);
	}

	return count;
}

/***** (HY) Subroutine that actually does correlation function *****/
/***** (HY-20060712) The original routine was wrong. Re-code. *****/
void xcorr::ctc(double xc[], double xcn[], double xce[], int Nt,const s_x xd1[],const s_x xd2[],int N,int dim){
	// (HY-20060717)
	// delete tau[], which was never used in the routine.
	// add xcn[], size of Ncorri, number of elements
	// add xce[], size of Ncorri, mean standard deviation
	// (HY-20060712)
	// tau[] is tc[], with a size of Ncorr
	// xc[] has a size of Ncorri
	int Ncorri = (int) floor(Tcorr/tcorr); // (HY-20060712)
	int xn[Ncorri];    // (HY-20060712), was int xn[Ncorr]
	double xm1[Ncorri];
	double xm2[Ncorri];
	double xt[Ncorri]; // (HY-20060712), was xt[Ncorr]
	double xte[Ncorri]; // (HY-20060717), errors
	for(int i=0;i<Ncorri;i++) {
		xc[i]=xce[i]=xcn[i]=0.0;
		xm1[i]=xm2[i]=0.0;
		xn[i]=0;
		xt[i]=xte[i]=0.0;
	}

	int k;
    for(int i=0;i<Nt-1;i++){
		spinner(i,Nt-1,25); // (HY) This is the progress bar on the screen
		for(int j=i+1;j<N && j<(Ncorri+i);j++){
			if(!isfinite(gsl_vector_get(xd1[i].c,0))) continue;
			if(!isfinite(gsl_vector_get(xd2[j].c,0))) continue;
			k = j-i-1;
			xn[k]++;
			xt[k]+=gsl_vector_get(xd1[i].c,0)*gsl_vector_get(xd2[j].c,0);
			xte[k]+= gsl_vector_get(xd1[i].c,0)*gsl_vector_get(xd2[j].c,0)
				   * gsl_vector_get(xd1[i].c,0)*gsl_vector_get(xd2[j].c,0);
			xm1[k]+=gsl_vector_get(xd1[i].c,0);
			xm2[k]+=gsl_vector_get(xd2[j].c,0);
		}
	}

	for(int i=0;i<Ncorri;i++){
		// (HY) This average scheme is compensated normalization.
		//      It is powerful for non-ergodic traces.
		if(xn[i]>0) {
			xc[i]=xt[i]/xn[i]-xm1[i]*xm2[i]/xn[i]/xn[i];
			xcn[i] = (double) xn[i]; // (HY-20060717) add
			xce[i] = xte[i]/xn[i]-xt[i]*xt[i]/xn[i]/xn[i]; // (HY-20060717) add
		}
	}

}

/***** (HY-20060827) Subroutine that computes the <x1(0)x2(t)> correlation function *****/
void xcorr::ctc_xx(double xc[], double xcn[], double xce[], int Nt,
    const s_x xd1[],const s_x xd2[],int N,int dim, double x1, double x2){
    // xcn[], size of Ncorri*Nxcorr*Nxcorr, number of elements
    // xce[], size of Ncorri*Nxcorr*Nxcorr, mean standard deviation
    // xc[] has a size of Ncorri*Nxcorr*Nxcorr
    int Ncorri = (int) floor(Tcorr/tcorr); // (HY-20060712)
    int xn[Ncorri];    // number of correlated points in xn[k]
    double xm1[Ncorri];
    double xm2[Ncorri];
    double dx = (x2-x1)/(double)Nxcorr; // resolution of grids on the x1-x2 plane
    double *xt, *xte;
    //double xt[Ncorri]; // (HY-20060712), was xt[Ncorr]
    //double xte[Ncorri]; // (HY-20060717), errors
    for(int i=0;i<Ncorri;i++) {
        xc[i]=xce[i]=xcn[i]=0.0;
        xm1[i]=xm2[i]=0.0;
        xn[i]=0;
    }
    xt = (double *) calloc( Ncorri*Nxcorr*Nxcorr, sizeof(double) );
    xte = (double *) calloc( Ncorri*Nxcorr*Nxcorr, sizeof(double) );

    int ii,jj,kk; // ii and jj will be indices for x1 and x2
    for(int i=0;i<Nt-1;i++){
        spinner(i,Nt-1,25); // (HY) This is the progress bar on the screen
        for(int j=i;j<N && j<(Ncorri+i);j++){
            if(!isfinite(gsl_vector_get(xd1[i].c,0))) continue;
            if(!isfinite(gsl_vector_get(xd2[j].c,0))) continue;
            ii = (int) floor((gsl_vector_get(xd1[i].c,0)-x1)/dx);
            jj = (int) floor((gsl_vector_get(xd2[j].c,0)-x1)/dx);
            kk = j-i; // time point along the Ncorri elements
            xn[kk]++;
            if ( ii >= Nxcorr || jj >= Nxcorr || ii < 0 || jj < 0 ) continue;
            //xt[kk]+=gsl_vector_get(xd1[i].c,0)*gsl_vector_get(xd2[j].c,0);
            xt[ii+Nxcorr*(jj+Nxcorr*kk)]++;
            /* HY-20060827, have not figured out how to compute error bars. Comment out
            xte[kk]+= gsl_vector_get(xd1[i].c,0)*gsl_vector_get(xd2[j].c,0)
                   * gsl_vector_get(xd1[i].c,0)*gsl_vector_get(xd2[j].c,0);
            xm1[kk]+=gsl_vector_get(xd1[i].c,0);
            xm2[kk]+=gsl_vector_get(xd2[j].c,0);
            */
        }
    }

    for(int i=0;i<Ncorri;i++){
        // for the i-th time delay in <x1(0)x2(i)> ...
        if(xn[i]>0) { // we have data for this time delay
            for (int j=0; j<Nxcorr; j++)
                for (int k=0; k<Nxcorr; k++)
                    xc[k+Nxcorr*(j+Nxcorr*i)]=(double) xt[k+Nxcorr*(j+Nxcorr*i)]/xn[i];
            xcn[i] = (double) xn[i]; // (HY-20060717) add
            //xce[i] = xte[i]/xn[i]-xt[i]*xt[i]/xn[i]/xn[i]; // HY-20060827, error function yet defined
        }
    }
    
    free(xt);
    free(xte);

}

void xcorr::ctc_pos(double tau[],double xc[],int Nt,const s_x xd1[],const s_x xd2[],int N,int dim, double X1, double X2){
	// (HY-20060711)
	// Nt, number of data points in the correlation function, given by Ncorr
	// xd1[], distance time trajectory 1
	// xd2[], distance time trajectory 2
	// N, number of data points in xd1[] and xd2[]
	// xc[], x fluctuation correlation function
	// tau[], not used
	// X1 and X2, range that defines a valid x so that X1 <= x_valid <= X2
	
	// initialization
	for(int i=0;i<Nt;i++) xc[i]=0;

	/***** (HY) Next look computes mean distance, <x> *****/
	double xt=0,xm1=0,xm2=0,xn=0;
	for(int i=0;i<N;i++){
		if(!isfinite(gsl_vector_get(xd1[i].c,0))) continue;
		if(!isfinite(gsl_vector_get(xd2[i].c,0))) continue;
		if(!((gsl_vector_get(xd1[i].c,0)>X1)&&(gsl_vector_get(xd1[i].c,0)<X2))) continue;
		xn++; // (HY) number of valud x points in xd1
		xm1+=gsl_vector_get(xd1[i].c,0); // (HY) sum of all valid x points
	}
	xm1/=xn; // (HY) mean value of x, <x>

	/***** (HY) Next i-j loop computes correlation <x1 * x2> *****/
	for(int i=0;i<Nt;i++){
		spinner(i,Nt,25);
		xt=0; // (HY) xt is the element in time correaltion function for the i-th time point
		      //      on the correlation time axis
		xn=0;
		xm2=0;
		for(int j=0;j<N-i;j++){
			if(!isfinite(gsl_vector_get(xd1[j].c,0))) continue;
			if(!isfinite(gsl_vector_get(xd2[j+i].c,0))) continue;
			if(!((gsl_vector_get(xd1[j].c,0)>X1)&&(gsl_vector_get(xd1[j].c,0)<X2))) continue;
			xn++; // (HY) number of valid data entry
			xt+=gsl_vector_get(xd1[j].c,0)*gsl_vector_get(xd2[j+i].c,0);
			xm2+=gsl_vector_get(xd2[j+i].c,0);
		}
		xm2/=xn;
		xc[i]=xt/xn-xm1*xm2; // (HY) xc[] stores the resulting fluctuation correlation
	}
}

double xcorr::xcc(double tau[],double xc[],int Nt,const s_x xd1[],const s_x xd2[],int N,int dim){
    for(int i=0;i<Nt;i++) xc[i]=0;
	//for(int i=0;i<Nt;i++) exc[i]=0;

    for(int i=0;i<N;i++){
        spinner(i,N,25);
        for(int j=i;j<N;j++){
//			if(xd1[i].tmin-xd2[j].tmax>tau[Nt-1]) continue;
//			if(xd2[j].tmin-xd1[i].tmax>tau[Nt-1]) break;
            xci(xc,tau,Nt,xd1[i],xd2[j],dim);
        }
    }

    for(int i=0;i<Nt;i++) xc[i]/=xd1[N-1].tmax-tau[i];

    double m1=0,m2=0;
    for(int i=0;i<N;i++){
        m1+=(xd1[i].tmax-xd1[i].tmin)*(gsl_vector_get(xd1[i].c,dim-1));
        m2+=(xd2[i].tmax-xd2[i].tmin)*(gsl_vector_get(xd2[i].c,dim-1));
    }
    m1/=xd1[N-1].tmax;
    m2/=xd2[N-1].tmax;

    for(int i=0;i<Nt;i++) xc[i]-=m1*m2;

    return 0;
}

void xcorr::xci(double xc[],double tau[],int Nt, s_x x1, s_x x2, int dim){
    s_x *a,*b;
    bool switched;

    if( x1.tmax-x1.tmin < x2.tmax-x2.tmin ){
        a=&x2;
        b=&x1;
        switched=true;
    } else{
        a=&x1;
        b=&x2;
        switched=false;
    }

    double cf=gsl_vector_get(a->c,dim-1)*gsl_vector_get(b->c,dim-1);

    for(int i=0;i<Nt;i++){
        if(!switched) tau[i]=-tau[i];
        if ( tau[i] <= a->tmin-b->tmax ) xc[i]+=0;
        else if( tau[i] <= a->tmin-b->tmin ) xc[i]+=cf*(b->tmax+tau[i]-a->tmin);
        else if( tau[i] <= a->tmax-b->tmax ) xc[i]+=cf*(b->tmax-b->tmin);
        else if( tau[i] <= a->tmax-b->tmin ) xc[i]+=cf*(a->tmax-b->tmin-tau[i]);
        if(!switched) tau[i]=-tau[i];
    }
}

double xcorr::xcorr_e(double tc[], double xct[], int N){
	double max=xct[1];

	for(int i=0;i<N;i++)
		if(xct[i]<max/exp(1.0)) return tc[i];

	return tc[N-1];
}

double xcorr::xcorr_int(double tc[], double xct[], int N){
	double max=xct[0];

	double ii=0;
	for(int i=0;xct[i]/max>1e-3;i++)
		ii+=xct[i];

	return ii*(tc[1]-tc[0]);
}

double xcorr::xcorr_int2(double tc[], double xct[], int N){
    double ii=0;

	ii+=gsl_pow_2(xct[0])*tc[0];
    for(int i=1;i<N;i++){
		if(isnan(xct[i])) break;
        ii+=gsl_pow_2(xct[i])*(tc[i]-tc[i-1]);
	}

    return ii;
}

inline double xcorr::Gaussian(double x, double s){
	return exp(-x*x/(2*s*s))/(s*pow(2*M_PI,0.5));
}
