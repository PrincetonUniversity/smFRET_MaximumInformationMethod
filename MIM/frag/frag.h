#include <libsmurf.h>
#include <fretdata.h>
#include <fstream>
#include <iostream>
#include <stdint.h>
#include <unistd.h>
#include <string>
using namespace std;

/*class s_xc{
public:
	s_xc(){
		x=0;
		sx=0;
		tmin=0;
		tmax=0;
	}

	double x;
	double sx;
	double tmin;
	double tmax;
}
*/
s_x read_dist(fstream& f, double t1, double t2, s_cal b);
void find_enc(double test[], double enc[], double times[], int n);
void sort_times(double[],int);
void print_times(fstream& tf, ostream& of, s_cal b, double times[], int n);
void find_split(fstream& tin,double tmin,double tmax,double times[],uint32_t& ltimes, s_cal b);
s_cal calibrate(fstream &fin, fstream& tin, uint32_t lt);
void mim(fstream &f, double times[], uint32_t& ltimes, double a, s_cal b);
double mistep(fstream &f, double t1, double a, s_cal b);
uint32_t process_data(fstream& data, fstream& processed);
void usage();
