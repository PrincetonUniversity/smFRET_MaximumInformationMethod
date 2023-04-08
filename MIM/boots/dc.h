#ifndef DC_H
#define DC_H

#include <iostream>
#include <fstream>
#include <string>

class dc{
public:
	dc(int dim=1);
	
	void set_ydata(gsl_vector*,int dim=1);
	void set_xdata(gsl_vector*);
	void set_stddata(gsl_vector*);
	void set_xdc();

}
