#include "fretdata.h"
#include <unistd.h> // HY-20060821

s_fretdata::s_fretdata(string datafile, string procfile){
    rawfilename=datafile;
    procfilename=procfile; // procfilename is declared public
    init();
}

s_fretdata::s_fretdata(string datafile){
    rawfilename=datafile;
    procfilename=datafile+".proc";
    init();
}

void s_fretdata::init(){
    //    processed=tmpfile();
	//lines below were removed so that temp files are put in local
	//directory on the current computer not on nfs...
    //char cwd[512]; // current working directory
    //getcwd(cwd, 512); // HY-20060821
    procfilename=tempnam("/tmp",NULL); // HY-20060821
    //procfilename=tmpnam(NULL);
    processed.open(procfilename.c_str(),ios::binary|ios::out|ios::in|ios::trunc);
    b=smurf_process_fret(rawfilename,processed);

    rnd = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rnd,time(NULL));
    
    CONFIDENCE1 = 0.95;
    CONFIDENCE2 = 0.95;
        
    coxOakesI = -1;
    coxOakesII = -1;
    coxOakesIII=-1;
    coxOakesIV = -1;
//    processed.open(procfilename.c_str(),ios::binary|ios::in);
}

s_fretdata::~s_fretdata(){
    processed.close();
    unlink(procfilename.c_str());
}

void s_fretdata::cross(double E){
    b=smurf_fret_cross(processed,b.ntl,E);
}

s_fr s_fretdata::get_photon(uint32_t p){
    return read_num(processed,p);
}

//Apply the Maximum Information Method to a trajectory contained in the file
//open in processed, with calibration information b and desired standard deviation a.
void s_fretdata::mim(double** times, uint32_t& ltimes, double a){
//    times[0]=0;

//    for(int ltimes=1;times[ltimes-1]<b.ta_bleach;ltimes++){
//        times[ltimes]=mistep(processed,times[ltimes-1],a,b);
//    }
}

//Given a start time t1, find the end time t2 given std a and data processed and b.
s_x s_fretdata::mistep(double t1, double a){
    uint32_t cnum=::findtime(processed,b.nt-1,t1);
    s_fr fp=read_num(processed,cnum);
    cnum++;
    s_fr lp=read_num(processed,cnum);
    s_fr lpm;
    double x,std,tt;
    s_x rv;
    while(lp.time<=b.ta_bleach){
        x=smurf_MLE2B(lp.cpcd-fp.cpcd,lp.cpca-fp.cpca,b);
        std=smurf_STD2B(x,lp.time-t1,b);

        if(std<a){
            cnum--;
            lpm=read_num(processed,cnum);
            x=smurf_MLE2B(lpm.cpcd-fp.cpcd,lpm.cpca-fp.cpca,b);
            std=smurf_STD2B(x,lp.time-t1,b);

            if((std>a)||(isnan(x))||(isnan(std))){
                rv.x=smurf_MLE2B(lp.cpcd-fp.cpcd,lp.cpca-fp.cpca,b);
                rv.std=smurf_STD2B(rv.x,lp.time-t1,b);
                rv.tmin=t1;
                rv.tmax=lp.time;
                return rv;
            }

            tt=t1+pow(smurf_STD2B(x,1,b)/a,2); //Information Rate
//            if(tt>lp.time) return -1;
            if(tt>b.ta_bleach) break;

            rv.tmin=t1;
            rv.tmax=tt;
            rv.x=smurf_MLE2B(lpm.cpcd-fp.cpcd,lpm.cpca-fp.cpca,b);
            rv.std=smurf_STD2B(rv.x,tt-t1,b);
            return rv;
        }

        cnum++;
        lp=read_num(processed,cnum);
    }
        
    rv.tmin=t1;
    rv.tmax=b.ta_bleach;
    rv.x=smurf_MLE2B(lpm.cpcd-fp.cpcd,lpm.cpca-fp.cpca,b);
    rv.std=smurf_STD2B(rv.x,b.ta_bleach-t1,b);
    gsl_vector_set(rv.c,0,rv.x);
    gsl_matrix_set(rv.cS,0,0,rv.std);
    return rv;
}

s_x s_fretdata::mistepe(double t1, double a){
    uint32_t cnum=::findtime(processed,b.nt-1,t1);
    s_fr fp=read_num(processed,cnum);
    cnum++;
    s_fr lp=read_num(processed,cnum);
    s_fr lpm;
    double x,std,tt;
    s_x rv;

    while(lp.time<=b.ta_bleach){
        x=smurf_MLEEB(lp.cpcd-fp.cpcd,lp.cpca-fp.cpca,b);
        std=smurf_STDEB(x,lp.time-t1,b);

        if(std<a){
            cnum--;
            lpm=read_num(processed,cnum);
            x=smurf_MLEEB(lpm.cpcd-fp.cpcd,lpm.cpca-fp.cpca,b);
            std=smurf_STDEB(x,lpm.time-t1,b);

            if((std>a)||(x<0)||(isnan(std))){
                rv.x=smurf_MLEEB(lp.cpcd-fp.cpcd,lp.cpca-fp.cpca,b);
                rv.std=smurf_STDEB(rv.x,lp.time-t1,b);
                gsl_vector_set(rv.c,0,rv.x);
                gsl_matrix_set(rv.cS,0,0,pow(rv.std,2.0));
                rv.tmin=t1;
                rv.tmax=lp.time;
                return rv;
            }

            tt=t1+pow(smurf_STDEB(x,1,b)/a,2); //Information Rate
//            if(tt>lp.time) return -1;
            if(tt>b.ta_bleach) break;

            rv.tmin=t1;
            rv.tmax=tt;
            rv.x=smurf_MLEEB(lpm.cpcd-fp.cpcd,lpm.cpca-fp.cpca,b);
            rv.std=smurf_STDEB(rv.x,tt-t1,b);
            gsl_vector_set(rv.c,0,rv.x);
            gsl_matrix_set(rv.cS,0,0,pow(rv.std,2.0));
            return rv;
        }

        cnum++;
        lp=read_num(processed,cnum);
    }
    rv.tmin=t1;
    rv.tmax=b.ta_bleach;
    rv.x=smurf_MLEEB(lpm.cpcd-fp.cpcd,lpm.cpca-fp.cpca,b);
    rv.std=smurf_STDEB(rv.x,b.ta_bleach-t1,b);
    gsl_vector_set(rv.c,0,rv.x);
    gsl_matrix_set(rv.cS,0,0,pow(rv.std,2.0));
    return rv;
}

//Given a middle time tm, find the end times t1 and t2 given std a and
//data processed and b.
s_x s_fretdata::mistepm(double tm, double a){
    uint32_t cnumm=::findtime(processed,b.nt-1,tm);
    uint32_t cnumf=cnumm-1;
    uint32_t cnuml=cnumm+1;
    s_fr fp=read_num(processed,cnumf);
    s_fr lp=read_num(processed,cnuml);
    s_fr fpc,lpc;
    s_fr fpm,lpm;
    double x,std,tt;
    s_x rv;

    while(lp.time<=b.ta_bleach && fp.time>=0 ){
        x=smurf_MLE2B(lp.cpcd-fp.cpcd,lp.cpca-fp.cpca,b);
        std=smurf_STD2B(x,lp.time-fp.time,b);

        if(std<a){
//            cnum--;
            if( lp.time-tm > tm-fp.time ){
                cnuml--;
                lpm=read_num(processed,cnuml);
                x=smurf_MLE2B(lpm.cpcd-fp.cpcd,lpm.cpca-fp.cpca,b);
                std=smurf_STD2B(x,lp.time-fp.time,b);

                if((std>a)||(isnan(x))||(isnan(std))){
                    rv.x=smurf_MLE2B(lp.cpcd-fp.cpcd,lp.cpca-fp.cpca,b);
                    rv.std=smurf_STD2B(rv.x,2*(lp.time-tm),b);
                    rv.tmin=2*tm-lp.time;
                    rv.tmax=lp.time;
                    return rv;
                }

                tt=pow(smurf_STD2B(x,1,b)/a,2); //Information Rate

                rv.tmin=tm-tt/2;
                rv.tmax=tm+tt/2;
                if( rv.tmax > b.ta_bleach || rv.tmin < 0 ){
                    cout << "Surprise!" << endl;
                    break;
                }
                rv.x=smurf_MLE2B(lpm.cpcd-fp.cpcd,lpm.cpca-fp.cpca,b);
                rv.std=smurf_STD2B(rv.x,tt,b);
                return rv;
            } else { //lp.time-tm < tm-fp.time
                cnumf++;
                fpm=read_num(processed,cnumf);
                x=smurf_MLE2B(lp.cpcd-fpm.cpcd,lp.cpca-fpm.cpca,b);
                std=smurf_STD2B(x,lp.time-fp.time,b);

                if((std>a)||(isnan(x))||(isnan(std))){
                    rv.x=smurf_MLE2B(lp.cpcd-fp.cpcd,lp.cpca-fp.cpca,b);
                    rv.std=smurf_STD2B(rv.x,2*(tm-fp.time),b);
                    rv.tmin=fp.time;
                    rv.tmax=2*tm-fp.time;
                    return rv;
                }

                tt=pow(smurf_STD2B(x,1,b)/a,2); //total time from information rate

                rv.tmin=tm-tt/2;
                rv.tmax=tm+tt/2;
                if( rv.tmax > b.ta_bleach || rv.tmin < 0 ){
                    cout << "Surprise!" << endl;
                    break;
                }
                rv.x=smurf_MLE2B(lp.cpcd-fpm.cpcd,lp.cpca-fpm.cpca,b);
                rv.std=smurf_STD2B(rv.x,tt,b);
                return rv;
            }
        }

        lpc=read_num(processed,cnuml+1);
        fpc=read_num(processed,cnumf-1);
        if(lpc.time-tm < tm-fpc.time){
            if(cnuml==b.ta_bleach){
                rv.tmax=b.ta_bleach+1;
                break;
            }
            cnuml++;
            lp=lpc;
        } else {
            if(cnumf==0){
                rv.tmin=-1;
                break;
            }
            cnumf--;
            fp=fpc;
        }
    } //while

    if(rv.tmin<0){
        rv.tmin=0;
        rv.tmax=2*tm;
        fp=read_num(processed,0);
        rv.x=smurf_MLE2B(lp.cpcd-fp.cpcd,lp.cpca-fp.cpca,b);
        rv.std=smurf_STD2B(rv.x,2*tm,b);
    } else if(rv.tmax>b.ta_bleach){
        rv.tmax=b.ta_bleach;
        rv.tmin=2*tm-b.ta_bleach;
        lp=read_num(processed,cnumf);
        rv.x=smurf_MLE2B(lp.cpcd-fp.cpcd,lp.cpca-fp.cpca,b);
        rv.std=smurf_STD2B(rv.x,2*tm,b);
    }

    return rv;
}

s_x s_fretdata::mipstep(double t1, double a, int n){
    // t1 is the start time for a constant-infomration measurement
    // a is the relative error
    // n is the dimensionality. For distance FRET, n = 1
    s_x ts=mistep(t1,a); // Given a start time t1, find the end time t2 given std a
                         // This routine uses the original maximum-information formula
    s_x xp=xpgen(ts.tmin,ts.tmax,n); // Calculate x and x-prime with standard deviations for a given time interval
                                     // This subroutine uses Taylor expanded series

    gsl_vector_set(xp.c,0,ts.x);
    gsl_matrix_set(xp.cS,0,0,pow(ts.std,2));
    return xp;
}

s_x s_fretdata::miestep(double t1, double a, int n){
    s_x ts=mistepe(t1,a);
//    s_x xp=xpgen(ts.tmin,ts.tmax,n);
    return ts;
//    gsl_vector_set(xp.c,0,ts.x);
//    gsl_matrix_set(xp.cS,0,0,pow(ts.std,2));
//    return xp;
}

s_x s_fretdata::box(double t1, double t2){
    if(t1<0) t1=0;
    if(t2>b.ta_bleach) t2=b.ta_bleach;
    uint32_t n1=::findtime(processed, b.nt-1,t1);
    uint32_t n2=::findtime(processed, b.nt-1,t2);

    s_fr p1=read_num(processed,n1);
    s_fr p2=read_num(processed,n2);

    s_x rv;
    rv.tmin=t1;
    rv.tmax=t2;
    rv.x=smurf_MLE2B(p2.cpcd-p1.cpcd,p2.cpca-p1.cpca,b);
    gsl_vector_set(rv.c,0,rv.x);
    rv.std=smurf_STD2B(rv.x,t2-t1,b);
    gsl_matrix_set(rv.cS,0,0,pow(rv.std,2));
    return rv;
}

s_x s_fretdata::boxe(double t1, double t2){
    if(t1<0) t1=0;
    if(t2>b.ta_bleach) t2=b.ta_bleach;
    uint32_t n1=::findtime(processed, b.nt-1,t1);
    uint32_t n2=::findtime(processed, b.nt-1,t2);

    s_fr p1=read_num(processed,n1);
    s_fr p2=read_num(processed,n2);

    s_x rv;
    rv.tmin=t1;
    rv.tmax=t2;
    rv.x=smurf_MLEEB(p2.cpcd-p1.cpcd,p2.cpca-p1.cpca,b);
    rv.std=smurf_STDEB(rv.x,t2-t1,b);
    return rv;
}

//Calculate x and x-prime with standard deviations for a given time interval
//Taylor expansion of the x^MLE is used for this calculation so that 
//the derivatives can be readily calculated
s_x s_fretdata::xpgen(double t1, double t2,int n){
    int mlestatus=GSL_CONTINUE;

    //Indices of the first and last photons in the given time interval
    uint32_t n1=::findtime(processed, b.nt-1,t1);
    uint32_t n2=::findtime(processed, b.nt-1,t2);

    //Set up photon data set for analysis
    s_fr tp;
    xpparams p(n,n2-n1,t2-t1,b);
    for(uint i=n1;i<n2;i++){
        tp=read_num(processed,i);
        tp.time-=t1;
        p.ts[i-n1]=tp;
    }

    //Calculation of parameters
    //Declare solver space
    const gsl_multiroot_fsolver_type *T;
    gsl_multiroot_fsolver *s;
    T = gsl_multiroot_fsolver_hybrid;
    s = gsl_multiroot_fsolver_alloc(T, n);
    gsl_multiroot_function f = {&drdt_f, n, &p};
    double x_init[n];
    gsl_vector *x = gsl_vector_alloc(n);

    //Iterate until a good solution is achieved
    for(double dx=.001;mlestatus!=GSL_SUCCESS && fabs(dx)<3;dx*=-2.2){
        //Set initial guess
        s_x inittst=box(t1,t2);
        x_init[0]=1+dx;
        for(int i=1;i<n;i++) x_init[i]=0;
        for(int i=0;i<n;i++) gsl_vector_set(x,i,x_init[i]);
    
        gsl_multiroot_fsolver_set (s, &f, x);
    
        mlestatus=GSL_CONTINUE;
        for(int i=0;i<1000 && mlestatus==GSL_CONTINUE;i++){
            mlestatus = gsl_multiroot_fsolver_iterate (s);
            if (mlestatus) break;   //check if solver is stuck
            mlestatus = gsl_multiroot_test_residual (s->f, 1e-3);
        }
    }
    
    for(uint i=0;i<p.Ndim;i++) p.cs[i]=gsl_vector_get(s->x,i);

    //Calculate Fisher Information Matrix and its Inverse
    gsl_integration_workspace *w=gsl_integration_workspace_alloc(1000);
    double result=0;
    double error=0;
    gsl_function F;
    F.function=&JIntegrand;
    F.params=&p;

    gsl_matrix* J=gsl_matrix_calloc(n,n);
    gsl_matrix* S=gsl_matrix_calloc(n,n);
    for(int i=0;i<2*n-1;i++){
        p.currdim=i;
        gsl_integration_qag(&F,-p.T/2,p.T/2,0,1e-10,1000,6,w,&result,&error);
        for(int j=0;j<n;j++){
            for(int k=0;k<n;k++){
                if(j+k==i) gsl_matrix_set(J,j,k,result);
            }
        }
    }

    gsl_permutation* perm=gsl_permutation_calloc(n);
    int* signum=new int(0);
    gsl_linalg_LU_decomp(J,perm,signum);
    gsl_linalg_LU_invert(J,perm,S);
//gsl_matrix_fprintf(stdout,S,"%e");
//cout << endl;

    //Prepare return value
    s_x rv(n);
    if(mlestatus){
        cout << gsl_strerror(mlestatus) << endl;
        rv.error=mlestatus;
    }
    rv.tmin=t1;
    rv.tmax=t2;
    rv.x=p.cs[0];
    rv.std=sqrt(gsl_matrix_get(S,0,0));
    if(n>1){
        rv.xp=p.cs[1];
        rv.stdp=sqrt(gsl_matrix_get(S,1,1));
    }
    gsl_vector_memcpy(rv.c,s->x);
    gsl_matrix_memcpy(rv.cS,S);

    //Free up memory
    gsl_integration_workspace_free(w);
    gsl_multiroot_fsolver_free(s);
    gsl_vector_free(x);
    gsl_matrix_free(J);
    gsl_matrix_free(S);
    gsl_permutation_free(perm);
    delete signum;

    return rv;
}

double MLEIntegrand(double t, void * params){
    xpparams *p=(xpparams*) params;

    double x=0;
    for(uint i=0;i<p->Ndim;i++)
        x+=p->cs[i]*gsl_pow_int(t,i);

    double x5=pow(x,5);

    double common=6*x5/pow(1+x*x5,2)*pow(t,int(p->currdim));
    double id=(p->b.IdB)*(p->b.bd-1)/(p->b.bd);
    double ia=(p->b.IaB)*(1-p->b.ba)/(p->b.ba);
    return (id+ia)*common;
}

double JIntegrand(double t, void * params){
    xpparams *p=(xpparams*) params;

    double x=0;
    for(uint i=0;i<p->Ndim;i++)
        x+=p->cs[i]*gsl_pow_int(t,i);

    double x5=pow(x,5);
    double x6=x*x5;
    
    double common=36*pow(x5,2)/pow(1+x6,3)*pow(t,int(p->currdim));
    double id=p->b.IdB*pow(p->b.bd-1,2)/(p->b.bd*(1+p->b.bd*x6));
    double ia=p->b.IaB*pow(1-p->b.ba,2)/(p->b.ba*(x6+p->b.ba));

    return common*(ia+id);
}

//Calculates the derivative of the Log-likelihood
int drdt_f (const gsl_vector * dd, void *params, gsl_vector * f){
    xpparams *p=(xpparams*) params;

    for(uint i=0;i<p->Ndim;i++)
        p->cs[i]=gsl_vector_get(dd,i);

    //sum term
    double terms[p->Ndim];
    for(uint i=0;i<p->Ndim;i++) terms[i]=0;

    double x,x5,x6,fact,num,den;
    for(uint i=0;i<p->Np;i++){
        x=0;

        for(uint j=0;j<p->Ndim;j++)
            x+=p->cs[j]*gsl_pow_int(p->ts[i].time-p->T/2,j);

        x5=pow(x,5);
        x6=x5*x;

        if(p->ts[i].type==ACCEPTOR){
            num=6*(1-p->b.ba)*x5;
            den=(1+x6)*(x6+p->b.ba);
            fact=num/den;
        } else if(p->ts[i].type==DONOR){
            num=6*(p->b.bd-1)*x5;
            den=(1+x6)*(p->b.bd*x6+1);
            fact=num/den;
        }

        for(uint j=0;j<p->Ndim;j++)
            terms[j]+=fact*gsl_pow_int(p->ts[i].time-p->T/2,j);
    }

    //integral term
    gsl_integration_workspace *w=gsl_integration_workspace_alloc(1000);
    double result=0;
    double numerror=0;
    gsl_function F;
    p->currdim=0;
    F.function=&MLEIntegrand;
    F.params=p;
    int err;

    //gsl_error_handler_t* gsleh=
    gsl_set_error_handler_off();
    

    for(uint i=0;i<p->Ndim;i++){
        p->currdim=i;
        err=gsl_integration_qag(&F,-p->T/2,p->T/2,0,1e-10,1000,6,w,&result,&numerror);
                if(err>0 && err!=GSL_EROUND){
            cout << __FILE__ << ": " << __LINE__ << " - " << gsl_strerror(err) << endl;
            exit(1);
        }
        terms[i]-=result;
    }

//    gsl_set_error_handler(gsleh);

    for(uint i=0;i<p->Ndim;i++)
        gsl_vector_set(f,i,terms[i]);

    gsl_integration_workspace_free(w);

    return GSL_SUCCESS;
}

bool s_fretdata::isvalid(){
    if (b.Vers == 3) {
        //There is only acceptor data.
            coxOakesI = coxOakesBin(0,b.ta_bleach,ACCEPTOR,CONFIDENCE1,-3,-9.1,500,"1");
            coxOakesII = coxOakesBin(b.ta_bleach,get_photon(b.ntl-1).time,ACCEPTOR,CONFIDENCE2,0,0,500,"2");
            coxOakesIII = coxOakesBin((b.ta_bleach/3),b.ta_bleach,ACCEPTOR,CONFIDENCE1,-3,-9.1,500,"3");
            coxOakesIV = coxOakesBin((b.ta_bleach*0.666),b.ta_bleach,ACCEPTOR,CONFIDENCE1,-3,-9.1,500,"4");
            cout << "coxOakesI: " << coxOakesI << endl;
            cout << "coxOakesII: " << coxOakesII << endl;
            cout << "coxOakesIII: " << coxOakesIII << endl;
            cout << "coxOakesIV: " << coxOakesIV << endl;
            if(b.ba<2 || isnan(b.ba)) {
                cout << "Too Much Background Noise! " << "ba value is: "<< b.ba << endl;
                return false;
            }
            if(b.IaB<0 || isnan(b.IaB)) {
                cout << "IaB Error!" << endl;
                return false; 
            }
            return true;
        } else if (b.Vers == 4) {
        //There is only donor data.
            /*
        coxOakesI = coxOakesBin(0,b.td_bleach,DONOR,CONFIDENCE1,-3,-9.1,500,"1");
            cout << "coxOaksI: " << coxOakesI << endl;
            coxOakesII = coxOakesBin(b.td_bleach,get_photon(b.ntl-1).time,DONOR,CONFIDENCE2,0,0,500,"2");
        cout << "coxOaksII: " << coxOakesII << endl;
            coxOakesIII = coxOakesBin((b.td_bleach/3),b.td_bleach,DONOR,CONFIDENCE1,-3,-9.1,500,"3");
            cout << "coxOaksIII: " << coxOakesIII << endl;
            coxOakesIV = coxOakesBin((b.td_bleach*0.6666),b.td_bleach,DONOR,CONFIDENCE1,-3,-9.1,500,"4");
            cout << "coxOaksIV: " << coxOakesIV << endl;
            */

            if(b.bd<2 || isnan(b.bd)) {
                cout << "Too Much Background Noise! " << "bd value is: "<< b.bd << endl;
                return false; 
            }
            if(b.IdB<0 || isnan(b.IdB)) {
                cout << "IdB Error!" << endl;
                return false; 
            }
            cout << "building coxOakes File..." << endl;
            buildDist(0, b.td_bleach,DONOR,"1");
            //buildDist((b.td_bleach)*(0.6666),b.td_bleach,DONOR,"2");
            return true;
        } else if (b.Vers == 2) {    
        //Normal FRET with both channels.
      //buildDist(0, b.td_bleach,DONOR);
            if(b.ta_bleach+.1>b.td_bleach){cout << 1 <<endl; return false;}
            if(b.IaB<0 || isnan(b.IaB)){cout << 2 <<endl;  return false;}
            if(b.IdB<0 || isnan(b.IdB)){cout << 3 << endl; return false;}
            if(b.ba<2 || isnan(b.ba)){cout << 4 <<endl; return false;}
            if(b.bd<2 || isnan(b.bd)){ cout << 5 << endl; return false;}
            if(b.P<0 || isnan(b.P)){ cout << 6 << endl; return false;}
      cout << "building coxOakes File..." << endl; // (HY) 20060816, disabled for speed
      buildDist(b.ta_bleach,b.td_bleach,DONOR,"1");
        /*
        coxOakesI = coxOakesBin(b.ta_bleach, b.td_bleach, DONOR,CONFIDENCE1,-3,-9.1,500,"1");
        if(!coxOakesI) { 
             cout << "Region I: Failed" << endl;
                 return false;
        }
        cout << "Region I: Passed" << endl;
        coxOakesII = coxOakesBin(b.td_bleach,get_photon(b.ntl-1).time,DONOR,CONFIDENCE2,0,0,500,"2");
        if(!coxOakesII) {
             cout << "Region II: Failed" << endl;
                 return false;
        }
            cout << "Region II: Passed" << endl;
        coxOakesIII = coxOakesBin(b.td_bleach,get_photon(b.ntl-1).time,ACCEPTOR,CONFIDENCE2,0,0,500,"3");
        if(!coxOakesIII) {
          cout << "Region III: Failed" << endl;
          return false;
        }
            cout << "Region III: Passed" << endl;
        */
        return true;
        } else {
            cout << "pass through error" << endl;
            return false;
        }
}

void s_fretdata::print_log(string logfile){
    ofstream lout(logfile.c_str());
    if (b.Vers == 2) {
        //Normal FRET with both channels.    
            lout << b.td_bleach << ' ' << b.ta_bleach << endl;
            lout << b.IdB << ' ' << sqrt(b.vIdB) << endl;
            lout << b.IaB << ' ' << sqrt(b.vIaB) << endl;
            lout << b.bd  << ' ' << sqrt(b.vbd) << endl;
            lout << b.ba  << ' ' << sqrt(b.vba)  << endl;
            lout << coxOakesI  << ' ' << coxOakesII << endl;
            lout << coxOakesIII  << ' ' << coxOakesIV << endl;
        } else if (b.Vers == 3) {
        //There is only acceptor data.
            lout << '0' << ' ' << b.ta_bleach << endl;
            lout << b.IaB << ' ' << sqrt(b.vIaB) << endl;
            lout << b.ba  << ' ' << sqrt(b.vba)  << endl;
            lout << coxOakesI  << ' ' << coxOakesII << endl;
            lout << coxOakesIII  << ' ' << coxOakesIV << endl;
        } else if (b.Vers == 4) {
        //There is only donor data.
            lout << b.td_bleach << ' ' << '0' << endl;
            lout << b.IdB << ' ' << sqrt(b.vIdB) << endl;
            lout << b.bd  << ' ' << sqrt(b.vbd) << endl;
            lout << coxOakesI  << ' ' << coxOakesII << endl;
            lout << coxOakesIII  << ' ' << coxOakesIV << endl;
        }
}

void s_fretdata::print_stats(){
        if (b.Vers == 2) {
        //Normal FRET with both channels.
            cout << rawfilename << " loaded: " << b.nt << " photons." << endl;
            cout << "Acceptor survival time: " << b.ta_bleach << " s." << endl;
            cout << "Donor survival time:    " << b.td_bleach << " s." << endl;
            cout << "IdB =    " << b.IdB << " cps" << endl;
            cout << "IaB =    " << b.IaB << " cps" << endl;
            cout << "beta_d = " << b.bd << endl;
            cout << "beta_a = " << b.ba << endl;
        } else if (b.Vers == 3) {
        //There is only acceptor data.
            cout << rawfilename << " loaded: " << b.nt << " photons." << endl;
            cout << "Acceptor survival time: " << b.ta_bleach << " s." << endl;
            cout << "IaB =    " << b.IaB << " cps" << endl;
            cout << "beta_a = " << b.ba << endl;
        } else if (b.Vers == 4) {
        //There is only donor data.
            cout << rawfilename << " loaded: " << b.nt << " photons." << endl;
            cout << "Donor survival time:    " << b.td_bleach << " s." << endl;
            cout << "IdB =    " << b.IdB << " cps" << endl;
            cout << "beta_d = " << b.bd << endl;
        }
}

void s_fretdata::gen_cache(gsl_vector* t, gsl_vector* x, gsl_vector* s){
    uint32_t Nt=x->size;
    double inc=b.ta_bleach/Nt;
    s_x td;
    for(uint i=0;i<Nt;i++){
        td=mistepm(inc*(0.5+i),gsl_vector_get(s,i));
        gsl_vector_set(t,i,inc*(0.5+i));
        gsl_vector_set(x,i,td.x);
        gsl_vector_set(s,i,td.std);
    }
}

double s_fretdata::LL(trajectory traj){
    double ad,bd,cd,aa,ba,ca;

    s_fr last=read_num(processed,b.nt);
    s_fr td;
    for(uint i=0;i<b.nt;i++){
        td=read_num(processed,i);
        if(td.type==DONOR){
            ad+=log(Id(traj.distance(td.time)));
        } else if(td.type==ACCEPTOR){
            aa+=log(Ia(traj.distance(td.time)));
        }
    }

    double iia=0;
    double iid=0;
    ofstream tout("test.txt");
    for(uint i=0;i<traj.Nt;i++){
        iid+=Id(gsl_vector_get(traj.traj,i));
        tout << gsl_vector_get(traj.time,i) << ' ' << Id(gsl_vector_get(traj.traj,i)) << endl;
        iia+=Ia(gsl_vector_get(traj.traj,i));
    }
    bd=iid*traj.dt;
    ba=iia*traj.dt;

    cd=gsl_sf_lngamma(last.cpcd);
    ca=gsl_sf_lngamma(last.cpca);

    return ad+aa-bd-ba-cd-ca;
}

double s_fretdata::Ia(double d){
    double E=1.0/(1+pow(d,6));
    return (E)*(b.IaB-b.IaB/b.ba)+b.IaB/b.ba;
}

double s_fretdata::Id(double d){
    double E=1.0/(1+pow(d,6));
    return (1-E)*(b.IdB-b.IdB/b.bd)+b.IdB/b.bd;
}

bool s_fretdata::KS(double t1, double t2, int chan, int mul) {
        //The only purpose of this method is as a batch so that we can 
        //actually get back values from the KS test if we want to by calling the method below
        //This method hides the change and makes it so if you old code is still compatible.  
        double currentKS;
        currentKS = KSVal(t1,t2,chan,mul);
        if (currentKS < 1) {
            return true;
        } else {
            return false;
        }
}

double s_fretdata::KSVal(double t1, double t2, int chan, int mul) {
    s_fr p1=get_photon(::findtime(processed, b.ntl,t1));
    s_fr p2=get_photon(::findtime(processed, b.ntl,t2));
    
    //CDF(t)=n/N
    double D=0;
    double T=t2-t1;
    double N=p2.cpc(chan)-p1.cpc(chan)+1;
    double sig=mul*1.95/sqrt(N);

    double Dt=0;
    double F=0;
    double SN=0;
    s_fr p;
    
    for(uint32_t i=p1.cpc();i<p2.cpc();i++){
        p=get_photon(i);
        F=(p.time-t1)/T;
                SN=(p.cpc(chan)-p1.cpc(chan))/N;
                Dt=fabs(F-SN);
        if(Dt>D) D=Dt;//return false;
    }
        
        if (b.Vers == 2) {
            //Normal FRET with both channels.
            cout << D << '\t' << sig << endl;
            return D/sig;
        } else if ((b.Vers == 3) || (b.Vers == 4)) {
            //There is only one channel of data.
            return D*sqrt(N);
        } else {
            //YOU'RE fucked because it shouldn't get here.  There is a warning though because
            //if it does get here nothing is returned
            return 0;
        }
        
}

double s_fretdata::KSExp(double t1,double t2, int chan) {
        s_fr p1=get_photon(::findtime(processed, b.ntl,t1));
    s_fr p2=get_photon(::findtime(processed, b.ntl,t2));
        
    double D=0;
    double Dt=0;
    double F=0;
    double SN=0;
    s_fr m;
        vector<double> pdt;
        
        //Here we will read in all times between photons
        m = get_photon(p1.cpc());
        double previoust = m.time;
        double addthem = 0;
        for (uint32_t w=(p1.cpc()+1);w<p2.cpc();w++) {
            m=get_photon(w);
            //Throw out after pulsing
            if (m.time-previoust > 1e-6) {
                pdt.push_back(m.time-previoust);
                addthem += m.time-previoust;
            }
            previoust = m.time;
        }
        double mean = addthem/pdt.size();
        
        //Here we are sorting from shortest to longest interphoton distances
        sort(pdt.begin(),pdt.end());
        
        for(uint32_t i=0;i<pdt.size();i++){
                F=1-exp(-pdt[i]/mean);
                SN=((double)i)/pdt.size();
                Dt=fabs(F-SN);
                if(Dt>D) D=Dt;//return false;
    }
        
        if ((b.Vers == 3) || (b.Vers == 4)) {
            //There is only one channel of data.
            return D*sqrt((double)pdt.size());
        } else {
            //YOU'RE fucked because it shouldn't get here.  There is a warning though because
            //if it does get here nothing is returned
            return D*sqrt((double)pdt.size());
        }
        return 0;
}

bool s_fretdata::coxOakesBin(double t1, double t2, int chan, double confid, double mean, double mean2, uint32_t photons, string testnum) {
  s_fr p1=get_photon(::findtime(processed, b.ntl,t1));
  s_fr p2=get_photon(::findtime(processed, b.ntl,t2));

  s_fr m;
  uint32_t passChan = 0;
  uint32_t previous1 = p1.cpc();
  uint32_t previous2 = p1.cpc();
  double highestvalue1 = 0;
  double highestvalue2 = 0;
  s_fr fallThrough, s1, s2;
  double variance;
  if (photons == 1000) {
    variance = 3.0376;
  } else if(photons == 500) {
    variance = 2.2815;
  }

  for (uint32_t w=p1.cpc();w<p2.cpc();w++) {
    m=get_photon(w);
    if(m.type == chan) {
      passChan += 1;
      if (passChan%photons == 0) {
    s1 = get_photon(previous1);
    s2 = get_photon(w);
    double coxOakesV = coxOakes(s1.time,s2.time,chan,testnum);
    coxOakesV = coxOakesV-mean;
    if ((double)fabs(coxOakesV) > highestvalue1) {
      highestvalue1 = (double)fabs(coxOakesV);
    }
    previous1 = w;
      } 
      if (passChan%5000 == 0) {
    s1 = get_photon(previous2);
    s2 = get_photon(w);
    double coxOakesV = coxOakes(s1.time,s2.time,chan,testnum);
    coxOakesV = coxOakesV-mean2;
    if ((double)fabs(coxOakesV) > highestvalue2) {
      highestvalue2 = (double)fabs(coxOakesV);
    }
    previous2 = w;
      }
      if ((p2.cpc(chan)-p1.cpc(chan)-passChan) == photons) {
    s1 = get_photon(w);
    s2 = get_photon(p2.cpc()-1);
    double coxOakesV = coxOakes(s1.time,s2.time,chan,testnum);
    coxOakesV = coxOakesV-mean;
    if ((double)fabs(coxOakesV) > highestvalue1) {
      highestvalue1 = (double)fabs(coxOakesV);
    }
    break;
      }
    }
  }
  bool pass;
  //If there was no section with at least 500 photons then highestvalue will still be
  //set to zero so we should do the test again based on the number of photons that are there.
  if (highestvalue1 == 0) {
    highestvalue1 = (double)fabs(coxOakes(t1, t2, chan, testnum)+1.5);
    double sigValue = GaussianProb(confid,2.28);
    if (highestvalue1 < sigValue) {
      pass = true;
    } else {
      pass = false;
    }
  }

  cout << "Highest CoxOakes Value: " << highestvalue1 << endl;
  double numBins = floor((double)(p2.cpc(chan)-p1.cpc(chan))/photons)+1;
  double sigValue = GaussianProb(pow(confid,1/numBins), variance);
  cout << "SigValue: " << sigValue << endl;
  if (highestvalue1 < sigValue) {
    pass = true;
  } else {
    pass = false;
  }
  
  if (highestvalue2 == 0) {
    //This happens when the interval is shorter than 5000 photons
    //we will just use the shorter interval test alone.
    return pass; 
  }
  
  cout << "Highest 5K Value: " << highestvalue2 << endl;
  numBins = floor((double)(p2.cpc(chan)-p1.cpc(chan))/5000);
  sigValue = GaussianProb(pow(confid,1/numBins), 10.31);
  cout << "SigValue: " << sigValue << endl;
  if (highestvalue2 < sigValue) {
    pass = true;
  } else {
    pass = false;
  }
  return pass;
}
double s_fretdata::coxOakes(double t1, double t2, int chan, string testnum) {
        s_fr p1=get_photon(::findtime(processed, b.ntl,t1));
        s_fr p2=get_photon(::findtime(processed, b.ntl,t2));
        
        //We need to calculate the intensity between t1 and t2 based on number of photons.
        double N=p2.cpc(chan)-p1.cpc(chan);
        //cout << "Number of Photons: " << N << endl;
        //double inten = N/(p2.time-p1.time);
        
        //setup output for diagnostic files
        //string base = rawfilename;
        //base.erase(rawfilename.length()-2,2);
        //string dummys = base+testnum+"cox";
        //ofstream coxout(dummys.c_str());
        //dummys = base+testnum+"pdts";
        //ofstream pdtsout(dummys.c_str());
        
        s_fr m;
        vector<double> pdt;
        
        uint32_t start = p1.cpc();
        m = get_photon(start);
        start += 1;
        while(m.type != chan) {
            m = get_photon(start);
            start += 1;
        }
        
        double previoust = m.time;
        double addthem = 0;
        //Here we will read in all times between photons
        for (uint32_t w=(start+1);w<=p2.cpc();w++) {
            m=get_photon(w);
            //Throw out after pulsing
            if ((m.time-previoust > 10e-6) && (m.type == chan)) {
            //if (m.type == chan) {
                pdt.push_back(m.time-previoust-10e-6);
                addthem += m.time-previoust-10e-6;
                //pdtout << (m.time-previoust) << endl;
                previoust = m.time;
            }
        }
        
        double mean = addthem/pdt.size();
        //Here we are sorting from shortest to longest interphoton distances.
        sort(pdt.begin(),pdt.end());
        
        //for(uint32_t i=0;i<pdt.size();i++){
        //    pdtsout << pdt[i] << endl;
    //}
        
        //Here we are throwing out all interphoton distances larger than 3*mean 
        //so that we won't get errors due to the tail
        //vector<double>::iterator begin = pdt.begin();
        //vector<double>::iterator end = pdt.end();
        //for(uint32_t g=0;g<pdt.size();g++) {
        //    if (pdt[g] >= 3*mean) {
        //        pdt.erase(find(begin,end,pdt[g]),pdt.end()); 
        //        break;
        //    }
        //}
        N = pdt.size();
        
    //addthem = 0;

    //for(uint32_t w=0;w<pdt.size();w++) {
      //       addthem += pdt[w];
           //}
    //mean = addthem/pdt.size();

        double coxOaks = N;
        //double thisIter = 0;

        for(uint32_t i=0;i<pdt.size();i++){
       coxOaks += (1-pdt[i]/mean)*log(pdt[i]/mean);
       //thisIter = (1-pdt[i]/mean)*log(pdt[i]/mean);
       //coxout << thisIter << endl;
    }
        
        return pow((6/N),0.5)*(coxOaks/M_PI);
}


double s_fretdata::coxOakes2(double t1, double t2, int chan, string testnum) {
        s_fr p1=get_photon(::findtime(processed, b.ntl,t1));
        s_fr p2=get_photon(::findtime(processed, b.ntl,t2));
        
        //We need to calculate the intensity between t1 and t2 based on number of photons.
        double N=p2.cpc(chan)-p1.cpc(chan);
        //cout << "Number of Photons: " << N << endl;
        //double inten = N/(p2.time-p1.time);
        
        //setup output for diagnostic files
        //string base = rawfilename;
        //base.erase(rawfilename.length()-2,2);
        //string dummys = base+testnum+"pdt";
        //ofstream pdtout(dummys.c_str());
        //dummys = base+testnum+"pdts";
        //ofstream pdtsout(dummys.c_str());
        
    s_fr m;
        vector<double> pdt;
        
        uint32_t start = p1.cpc();
        m = get_photon(start);
        start += 1;
        while(m.type != chan) {
            m = get_photon(start);
            start += 1;
        }
        
        double previoust = m.time;
        double addthem = 0;
        //Here we will read in all times between photons
        for (uint32_t w=(start+1);w<=p2.cpc();w++) {
            m=get_photon(w);
            //Throw out after pulsing
            if ((m.time-previoust > 10e-6) && (m.type == chan)) {
            //if (m.type == chan) {
                pdt.push_back(m.time-previoust-10e-6);
                addthem += m.time-previoust-10e-6;
                //pdtout << (m.time-previoust) << endl;
                previoust = m.time;
            }
        }
        
        double mean = addthem/pdt.size();
        //Here we are sorting from shortest to longest interphoton distances.
        sort(pdt.begin(),pdt.end());
        N = pdt.size();
        
        double coxOaks = N;
        
        for(uint32_t i=0;i<pdt.size();i++){
            coxOaks += (1-pdt[i]/mean)*log(pdt[i]/mean);
    }
        
        return pow((6/N),0.5)*(coxOaks/M_PI);
}

void s_fretdata::buildDist(double t1, double t2, int chan, string testnum) {
        s_fr p1=get_photon(::findtime(processed, b.ntl,t1));
        s_fr p2=get_photon(::findtime(processed, b.ntl,t2));
        
        vector<double> vsam10, vsam50, vsam100, vsam500, vsam1000, vsam5000;

        string base = rawfilename;
        base.erase(rawfilename.length()-2,2);
        string dummys = base+"sam";
        ofstream samOut(dummys.c_str());
        s_fr m;
        uint32_t donChan = 0;
        
        uint32_t ps10 = p1.cpc();
        uint32_t ps50 = p1.cpc();
        uint32_t ps100 = p1.cpc();
        uint32_t ps500 = p1.cpc();
        uint32_t ps1000 = p1.cpc();
        uint32_t ps5000 = p1.cpc();
        
        for (uint32_t w=p1.cpc();w<p2.cpc();w++) {
            m=get_photon(w);
            if(m.type == DONOR) {
                donChan += 1;
                if (donChan%10 == 0) {
                    s_fr s1 = get_photon(ps10);
                    s_fr s2 = get_photon(w);
                    double coxOaks10 = coxOakes(s1.time,s2.time,chan,"1"); 
                    vsam10.push_back(coxOaks10);
                    ps10 = w;
                }
                if (donChan%50 == 0) {
                    s_fr s1 = get_photon(ps50);
                    s_fr s2 = get_photon(w);
                    double coxOaks50 = coxOakes(s1.time,s2.time,chan,"1"); 
                    vsam50.push_back(coxOaks50);
                    ps50 = w;
                }
                if (donChan%100 == 0) {
                    s_fr s1 = get_photon(ps100);
                    s_fr s2 = get_photon(w);
                    double coxOaks100 = coxOakes(s1.time,s2.time,chan,"1"); 
                    vsam100.push_back(coxOaks100);
                    ps100 = w;
                }
                if (donChan%500 == 0) {
                    s_fr s1 = get_photon(ps500);
                    s_fr s2 = get_photon(w);
                    double coxOaks500 = coxOakes(s1.time,s2.time,chan,"1"); 
                    vsam500.push_back(coxOaks500);
                    ps500 = w;
                }
                if (donChan%1000 == 0) {
                    s_fr s1 = get_photon(ps1000);
                    s_fr s2 = get_photon(w);
                    double coxOaks1000 = coxOakes(s1.time,s2.time,chan,"1"); 
                    vsam1000.push_back(coxOaks1000);
                    ps1000 = w;
                }
                if (donChan%5000 == 0) {
                    s_fr s1 = get_photon(ps5000);
                    s_fr s2 = get_photon(w);
                    double coxOaks5000 = coxOakes(s1.time,s2.time,chan,"1"); 
                    vsam5000.push_back(coxOaks5000);
                    ps5000 = w;
                }
            }
        }
        for(uint32_t h=0;h<vsam10.size();h++) {
            samOut << vsam10[h] << endl;
        }
        samOut << '*' << endl;
        for(uint32_t h=0;h<vsam50.size();h++) {
            samOut << vsam50[h] << endl;
        }
        samOut << '*' << endl;
        for(uint32_t h=0;h<vsam100.size();h++) {
            samOut << vsam100[h] << endl;
        }
        samOut << '*' << endl;
        for(uint32_t h=0;h<vsam500.size();h++) {
            samOut << vsam500[h] << endl;
        }
        samOut << '*' << endl;
        for(uint32_t h=0;h<vsam1000.size();h++) {
            samOut << vsam1000[h] << endl;
        }
        samOut << '*' << endl;
        for(uint32_t h=0;h<vsam5000.size();h++) {
            samOut << vsam5000[h] << endl;
        }
}

//double s_fretdata:ChiSqBin(double t1, double t2, int chan, double binSize) {
//    s_fr p1=get_photon(::findtime(processed, b.ntl,t1));
//    s_fr p2=get_photon(::findtime(processed, b.ntl,t2));
//    double curChi = 0;
//    double biggestChi = 0;
//    double endTime = 0;
//    for(int n=0;n<=(p2.time-p1.time)/binSize;n++) {
//        if ((p1.time + (n+1)*binSize) > p2.time) {
//            endTime = p2.time;
//        } else {
//            endTime = p1.time + (n+1)*binSize;
//        }
//        curChi = ChiSq(p1.time + n*binSize, endTime,DONOR);
//        if (curChi > biggestChi) {
//           biggestChi = curChi;
//        }
//    }
//    return biggestChi;
//}

/*
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>

double s_fretdata::PowerSpec(double t1, double t2, ) {
  int i, n = 100;
  double data[n];

  gsl_fft_real_wavetable * real;
  gsl_fft_halfcomplex_wavetable * hc;
  gsl_fft_real_workspace * work;

  for (i = 0; i < n; i++)
    {
      data[i] = 0.0;
    }

  for (i = n / 3; i < 2 * n / 3; i++)
    {
      data[i] = 1.0;
    }

  for (i = 0; i < n; i++)
    {
      printf ("%d: %e\n", i, data[i]);
    }
  printf ("\n");

  work = gsl_fft_real_workspace_alloc (n);
  real = gsl_fft_real_wavetable_alloc (n);

  gsl_fft_real_transform (data, 1, n, 
                          real, work);

  gsl_fft_real_wavetable_free (real);

  for (i = 11; i < n; i++)
    {
      data[i] = 0;
    }

  hc = gsl_fft_halfcomplex_wavetable_alloc (n);

  gsl_fft_halfcomplex_inverse (data, 1, n, 
                               hc, work);
  gsl_fft_halfcomplex_wavetable_free (hc);

  for (i = 0; i < n; i++)
    {
      printf ("%d: %e\n", i, data[i]);
    }

  gsl_fft_real_workspace_free (work);
  return 0;
}

*/

double s_fretdata::ChiSq(double t1, double t2, int chan, string testnum) {
        s_fr p1=get_photon(::findtime(processed, b.ntl,t1));
    s_fr p2=get_photon(::findtime(processed, b.ntl,t2));
        
        //setup output for diagnostic files
        string base = rawfilename;
        base.erase(rawfilename.length()-2,2);
        string dummys = base+testnum+"pdt";
        ofstream pdtout(dummys.c_str());
        dummys = base+testnum+"obs";
        ofstream obsOut(dummys.c_str());
        dummys = base+testnum+"cal";
        ofstream calcOut(dummys.c_str());
        
        //We need to calculate the intensity between t1 and t2 based on number of photons.
        double N=p2.cpc(chan)-p1.cpc(chan);
        //cout << "Number of Photons: " << N << endl;
        double inten = N/(p2.time-p1.time);
        
    s_fr m;
        vector<double> pdt;
        
        //Here we will read in all times between photons
        //Should add some code here so that the first photon from the correct channel.
        m = get_photon(p1.cpc());
        double previoust = m.time;
        double addthem = 0;
        for (uint32_t w=(p1.cpc()+1);w<p2.cpc();w++) {
            m=get_photon(w);
            //Throw out after pulsing
            if ((m.time-previoust > 1e-6) && (m.type == chan)) {
                pdt.push_back(m.time-previoust);
                addthem += m.time-previoust;
                previoust = m.time;
            }
        }
        
        double mean = addthem/pdt.size();
        //Here we are sorting from shortest to longest interphoton distances.
        sort(pdt.begin(),pdt.end());
        
        //Here we are throwing out all interphoton distances larger than 3*mean 
        //so that we won't get errors due to the tail
        vector<double>::iterator begin = pdt.begin();
        vector<double>::iterator end = pdt.end();
        for(uint32_t g=0;g<pdt.size();g++) {
            if (pdt[g] >= 3*mean) {
                pdt.erase(find(begin,end,pdt[g]),pdt.end()); 
                break;
            }
        }
        
        double binW = pdt[pdt.size()-1]/50;//sqrt(pdt.size());
        //cout << "Bin Width: " << binW << endl;
        int numBins = 50;//(int)(round(pdt[pdt.size()-1]/binW));
        //cout << "Number of Bins: " << numBins << endl;
        double obs[numBins];
        
        for (int q=1;q<=numBins;q++) {
            obs[q] = 0;
        }
        
        for(uint32_t i=0;i<pdt.size();i++){
            obs[(int)(round((pdt[i]/binW)+0.5))] += 1;
            pdtout << pdt[i] << endl;
    }
        
        for (int q=1;q<=numBins;q++) {
            obsOut << obs[q] << endl;
        }
        
        double chiSum = 0;
        for(int s=1;s<=numBins;s++) {
            double calc = N*inten*binW*exp(-inten*((s-1)*binW+binW/2));
            calcOut << calc << endl;
            chiSum += pow((obs[s]-calc),2)/calc;
        }
         
        
        //return (chiSum/numBins);
        //(chiSum/(numBins-1)),numBins-1
        double chiSq = chiSum/(numBins-1);
        //cout << "chiSq: " << chiSq << endl;
        //if (chiSq > 5) {
        //    chiSq = 5;
        //}
        //numBins-1
        return chiSq;
        
        //chiSqProb(chiSq,numBins-1);
}

double s_fretdata::Variance(double t1, double t2, int chan) {
    s_fr p1=get_photon(::findtime(processed, b.ntl,t1));
    s_fr p2=get_photon(::findtime(processed, b.ntl,t2));
    
    //We need to calculate the intensity between t1 and t2 based on number of photons.
    double N=p2.cpc(chan)-p1.cpc(chan);
    //cout << "Number of Photons: " << N << endl; 
    double variance = 0; 

    s_fr m;
    double binTime = 0 ;
    double chanCounts = 0;
    double previoust = p1.time;
    double previousCs = 0;
    double binInten = 0;
    double numBins = 0;
    for (uint32_t w=p1.cpc();w<p2.cpc();w++) {
        m=get_photon(w);
        if (m.type == chan) { 
            chanCounts += 1;
            if (chanCounts == 500) {
                binTime = m.time-p1.time;
                binInten = (N*binTime)/(p2.time-p1.time); 
            }
            if ((m.time-previoust > binTime) && (chanCounts >= 500)) {
                variance += fabs(chanCounts-previousCs-1)/pow(binTime,2);
                numBins += 1;
                previoust = m.time;
                previousCs = chanCounts;
            }
        }
    }

    return variance/numBins;
}

double s_fretdata::GaussianProb(double x, double variance) {
  //gsl_integration_workspace *w=gsl_integration_workspace_alloc(10000);
    double result;
    double error;
    size_t numFCalls;
    gsl_function F;
    gasParams p(variance);

    F.function = &gas_f;
    F.params = &p;
    double previous = 10;
    double bestValue = 0; 
    for (int y=150;y<=1000;y++) {
        gsl_integration_qng(&F,-(double)y/100,(double)y/100,1e-10,1e-6,&result,&error,&numFCalls);
    if (((double)fabs(result-x))<previous) {
      bestValue = (double)y/100;
      previous = (double)fabs(result-x);
    }
    }
    //gsl_integration_workspace_free(w);
    return bestValue;
}

double gas_f(double x, void * p) {
    gasParams *params = (gasParams *)p;
    double a = (params->A);
    double gass = (1/sqrt((a)*2*M_PI))*exp(-pow(x,2)/(2*a));
    return gass;
}

uint32_t s_fretdata::findtime(double t, bool before){
    return ::findtime(processed, b.nt-1, t, before);
}

s_fr s_fretdata::get_photon(double t, bool before){
    uint32_t ind=findtime(t, before);
    return get_photon(ind);
}

/*double s_fretdata::LL(s_x te){
    if(te.tmin<0) te.tmin=0;
    if(te.tmax>b.ta_bleach) te.tmax=b.ta_bleach;

    uint32_t n1=::findtime(processed, b.nt-1,te.tmin);
    uint32_t n2=::findtime(processed, b.nt-1,te.tmax);

    s_fr p1=read_num(processed,n1);
    s_fr p2=read_num(processed,n2);

    uint32_t na=p2.cpca-p1.cpca;
    uint32_t nd=p2.cpcd-p1.cpcd;

    double E=1/(1+pow(te.x,6));
    double Es=E-1/(1+pow(te.x+te.std,6));
    double Id=b.IdB*((1-1/b.bd)*(1-E)+1/b.bd);
    double Ids=b.IdB*((1-1/b.bd)*(1-E+Es)+1/b.bd)-Id;
    double Ia=b.IaB*((1-1/b.ba)*E+1/b.ba);
    double Ias=Ia-b.IaB*((1-1/b.ba)*(E-Es)+1/b.ba);
    double T=(te.tmax-te.tmin);

//    double elld=-pow(nd-Id*T,2)/(2*(nd+Ids*Ids*T*T))-.5*log(Ids*Ids*nd*(1/(Ids*Ids)+T*T/nd));
//    double ella=-pow(na-Ia*T,2)/(2*(na+Ias*Ias*T*T))-.5*log(Ias*Ias*na*(1/(Ias*Ias)+T*T/na));

//    cerr << E << ' ' << Es << ' ' << Ia << ' ' << Ias << ' ' << ella+elld << endl;

    double ella1=na*log(Ia*T)-Ia*T;
    double ella2=na*log((Ia+2*Ias)*T)-(Ia+2*Ias)*T;
    double ella3=na*log((Ia-2*Ias)*T)-(Ia-2*Ias)*T;
    double ella=.70*ella1+.15*(ella3+ella2)-gsl_sf_lngamma(na);

    double elld1=nd*log(Id*T)-Id*T;
    double elld2=nd*log((Id+2*Ids)*T)-(Id+2*Ids)*T;
    double elld3=nd*log((Id-2*Ids)*T)-(Id-2*Ids)*T;
    double elld=.70*elld1+.15*(elld3+elld2)-gsl_sf_lngamma(nd);
//cout << Ia << ' ' << Ias << ' ' << Id << ' ' << Ids << ' ' << T << ' ' << na << endl;
//cout << ella1 << ' ' << ella2 << ' ' << ella3 << endl;
    return elld+ella;
//    return na*log(Ia)-Ia*T+nd*log(Id)-Id*T-log(Ias)-log(Ids);
//cout << Id << ' ' << Ids << endl;

cout << .5+na/2 << ' ' << -2*Ia+Ias*Ias*T << ' ' << pow((-2*Ia+Ias*Ias*T)/2/Ias,2) << endl;

    double lla1=-pow(Ia/Ias,2)+log(Ias)-na*log(T-2*Ia/Ias/Ias)+(2+na)*log((-2*Ia+Ias*Ias*T)/Ias/Ias);
    double lla2=-log(4.0)-.5*log(2*M_PI)-3*log(-2*Ia+Ias*Ias*T)-gsl_sf_lngamma(1+na);
    double lla3=0,hg=0;

    if(na%2==0){
        lla3=log(4.0)+na*log(T)+2*log(Ias)+log((-2*Ia+Ias*Ias*T)/Ias)+gsl_sf_lngamma((1+na)/2)+log(gsl_sf_hyperg_1F1(.5+na/2,.5,pow((-2*Ia+Ias*Ias*T)/2/Ias,2)));
    } else{
        lla3=log(4.0)+log((double)na)+2*log(-2*Ia+Ias*Ias*T)+na*log(T)+gsl_sf_lngamma(na/2.0)+log(gsl_sf_hyperg_1F1(1.0+na/2.0,1.5,pow((-2*Ia+Ias*Ias*T)/2/Ias,2)));
    }


cout << .5+nd/2 << ' ' << -2*Id+Ids*Ids*T << ' ' << pow((-2*Id+Ids*Ids*T)/2/Ids,2) << endl;

    double lld1=-pow(Id/Ids,2)+log(Ids)-nd*log(T-2*Id/Ids/Ids)+(2+nd)*log((-2*Id+Ids*Ids*T)/Ids/Ids);
    double lld2=-log(4.0)-.5*log(2*M_PI)-3*log(-2*Id+Ids*Ids*T)-gsl_sf_lngamma(1+nd);
    double lld3=0;

    if(nd%2==0){
        lld3=log(4.0)+nd*log(T)+2*log(Ids)+log((-2*Id+Ids*Ids*T)/Ids)+gsl_sf_lngamma((1+nd)/2)+log(gsl_sf_hyperg_1F1(.5+nd/2,.5,pow((-2*Id+Ids*Ids*T)/2/Ids,2)));
    } else{
        lld3=log(4.0)+log((double)nd)+2*log(-2*Id+Ids*Ids*T)+nd*log(T)+gsl_sf_lngamma(nd/2.0)+log(gsl_sf_hyperg_1F1(1.0+nd/2.0,1.5,pow((-2*Id+Ids*Ids*T)/2/Ids,2)));
    }

//log(gsl_sf_hyperg_1F1(.5+na/2,.5,pow((-2*Ia+Ias*Ias*T)/2/Ias,2)));}
//gsl_sf_laguerre_n(



cout << lla1*lla3/lla2 << ' ' << lld1*lld3/lld2 << endl;
return na*log(T);
}*/
