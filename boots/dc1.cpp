#include <fstream>
#include <gsl/gsl_vector.h>
#include <hist.h>

int main(){
	gsl_vector* x=gsl_vector_alloc(100);
	gsl_vector* y=gsl_vector_alloc(100);
	gsl_vector* d=gsl_vector_alloc(100);
	gsl_vector_set_all(d,0.1);

	double a;
	ifstream zin("hist.out");
	for(int i=0;i<100;i++){
		zin >> a;
		gsl_vector_set(x,i,a);
		zin >> a >> a;
		gsl_vector_set(y,i,a);
		zin >> a >> a;
	}

	hist dc(1,100,HIST_SCALE_LOG);

	dc.setHistX(x,0);
	dc.setHistXAveDev(0.1,0);
	dc.setHistRaw(y);
	dc.set_mean_ey(0.01);
	dc.set_mean_ey2(0.01);
	dc.set_stepsize(1);
	dc.set_tolerance(0.1);

	dc.Print_MEM();

	dc.set_lambda(.01);
	dc.DoMEM(0);
}
