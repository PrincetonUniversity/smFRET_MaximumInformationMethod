#include <fstream>
#include <gsl/gsl_vector.h>
#include <hist.h>

int main(){
	gsl_vector* x1=gsl_vector_alloc(50*50);
	gsl_vector* x2=gsl_vector_alloc(50*50);
	gsl_vector* y=gsl_vector_alloc(50*50);
	gsl_vector* d=gsl_vector_alloc(50*50);
	gsl_vector_set_all(d,0.1);

	double a;
	ifstream zin("zs.txt");
	for(int i=0;i<50*50;i++){
		zin >> a;
		gsl_vector_set(x1,i,a);
		zin >> a;
		gsl_vector_set(x2,i,a);
		zin >> a;
		gsl_vector_set(y,i,a);
	}

	hist dc(2,50,HIST_SCALE_LINEAR);

	dc.setHistX(x1,0);
	dc.setHistX(x2,1);
	dc.setHistXAveDev(0.1,0);
	dc.setHistXAveDev(0.1,1);
	dc.setHistRaw(y);
	dc.set_mean_ey(1.0);
	dc.set_mean_ey2(1.0);
	dc.setHKernel2d(0.1);
	dc.set_stepsize(1);
	dc.set_tolerance(0.1);

	dc.Print_MEM();

	dc.set_lambda(0);
	dc.DoMEM(0);
}
