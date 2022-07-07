#include "friction.h"
#include "ltraj.h"

friction::friction(){
	init();
}

friction::friction(double gamma){
	init();
	setup(1);
	for(uint i=0;i<1;i++){
		set_param(i,gamma);
	}
}

friction::friction(double gamma, double nu1, double nu2, int N){
	init();
	set_constant(gamma,nu1,nu2,N);
	print("fout.out");
}

friction::~friction(){
	free_gdata();
}

friction::friction(const friction& nf){
	init();
	u1=nf.u1;
	du=nf.du;
	setup(nf.length());
	if(nf.gdata)
		gsl_vector_memcpy(gdata,nf.gdata);
}

double friction::g(double x){
	return g(n(x));
}

double friction::g(uint n){
	return get_param(n);
}

double friction::x(uint ni){
	return u(ni);
}

double friction::get_dx(){
	return du;
}

void friction::print(std::ostream& fout){
	for(uint i=0;i<length();i++){
		fout << x(i) << '\t' << g(i) << endl;
	}
}

void friction::update_from_params(){
}

double friction::LLR(ltraj* wrt, uint i1, uint i2){
	return wrt->LLR_Langevin(this,static_cast<friction*>(ltt),i1,i2);
}
	
friction* friction::clone(){
	return new friction(*this);
}

void friction::init(){
	gdata=NULL;
	name="friction";
}

void friction::free_gdata(){
	if(gdata)
		gsl_vector_free(gdata);
	gdata=NULL;
	params=NULL;
}

void friction::setup(int size){
	if(!gdata){
		gdata=gsl_vector_calloc(size);
	} else if(gdata->size!=uint(size)){
		free_gdata();
		gdata=gsl_vector_calloc(size);
	}
	params=gdata;
}
