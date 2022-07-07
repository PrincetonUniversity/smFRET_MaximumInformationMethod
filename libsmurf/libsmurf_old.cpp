#include <algorithm>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include "libsmurf.h"
#include <unistd.h>
#include <stdint.h>
using namespace std;

xt_type SMURF_XT=UNSET;

inline void read_Small_Endian(fstream& file, char* data, int size) {
        file.read(data, size);
    #ifdef __APPLE__
            int i = 0;
            int j = size - 1;
            while (i<j) {
                std::swap(data[i], data[j]);
                i++, j--;
            }
    #endif
}

inline void write_Small_Endian(fstream& file, char* data, int size) {
    #ifdef __APPLE__
            int i = 0;
            int j = size - 1;
            while (i<j) {
                std::swap(data[j], data[i]);
                i++, j--;
            }
    #endif
        file.write(data, size);
    #ifdef __APPLE__
            i = 0;
            j = size - 1;
            while (i<j) {
                std::swap(data[j], data[i]);
                i++, j--;
            }
    #endif
}

uint32_t smurf_loadwhich(const string fname) {
	fstream fin;
	fin.open(fname.c_str(),fstream::binary|fstream::in);

	uint32_t tm;
	uint16_t tt,tv;

	fin.seekg(0);
	read_Small_Endian(fin, (char*)&tm, sizeof(uint32_t));
	if(tm!=MAGIC) return NOT_SMURF;

	read_Small_Endian(fin, (char*)&tv, sizeof(uint16_t));

	read_Small_Endian(fin, (char*)&tt, sizeof(uint16_t));

	fin.close();

	switch(tt){
	case ft_time:
//		smurf_loadtt(fname,(s_tt*)data);
		return ft_time;
		break;
	case ft_fret:
//		smurf_loadfret2(fname,(s_fret2*)data);
		return ft_fret;
		break;
	default:
		return NOT_TT&NOT_FRET;
		break;
	}
}

//Check file type. Returns 0 if not a fret file. Otherwise, returns version number.
int smurf_isfret(const string fname) {
	//Check file type
	fstream fin;
	fin.open(fname.c_str(),fstream::binary|fstream::in);
	fin.seekg(0);
	uint32_t tm;
	uint16_t tv;
	uint16_t tt;

	read_Small_Endian(fin, (char*)&tm, sizeof(uint32_t));
	if(tm!=MAGIC){
                cout << tm << endl;
		cerr << "Not a smurf file! :)"  << endl;
		return 0;
	}

	read_Small_Endian(fin, (char*)&tv, sizeof(uint16_t));
	read_Small_Endian(fin, (char*)&tt, sizeof(uint16_t));
	
        if(tt!=ft_fret) {
		cerr << "Not a FRET file! :)" << endl;
		return 0;
	}

	uint32_t ld;
	uint32_t la;
	fin.seekg(8);
	read_Small_Endian(fin, (char*)&(ld),sizeof(uint32_t));
        read_Small_Endian(fin, (char*)&(la),sizeof(uint32_t));

	if((ld==0) && (la==0)){
		cerr << "Zero photons!" << endl;
		return 0;
	} else if (la == 0) {
                cerr << "File only contains donor channel data." << endl;
                tv = 4;
        } else if (ld == 0) {
                cerr << "File only contains acceptor channel data." << endl;
                tv = 3;
        }
	fin.seekg(1024+ld*sizeof(uint32_t));
	uint32_t td;
	read_Small_Endian(fin, (char*)&(td),sizeof(uint32_t));

	return tv;
}

uint32_t smurf_loadtt(const string fname, s_tt* data){
	fstream fin;
	fin.open(fname.c_str(),fstream::binary|fstream::in);

	uint32_t tm;
	uint16_t tt,tv;

	fin.seekg(0);
	read_Small_Endian(fin, (char*)&tm, sizeof(uint32_t));
	if(tm!=MAGIC) return NOT_SMURF;
	read_Small_Endian(fin, (char*)&tv, sizeof(uint16_t));
	read_Small_Endian(fin, (char*)&tt, sizeof(uint16_t));
	if(tt!=ft_time) return NOT_TT;

	fin.seekg(24);
	read_Small_Endian(fin, (char*)&(data->clock),sizeof(double));
	read_Small_Endian(fin, (char*)&(data->power),sizeof(double));

	fin.seekg(0,ios::end);
	int length=fin.tellg();
	data->l=(length-1024)/4;
	data->data=new uint32_t[data->l];


	fin.seekg(1024);
	read_Small_Endian(fin, (char*)data->data,data->l*sizeof(uint32_t));

	return 0;

}

uint32_t smurf_savett(const string fname, s_tt* data){
	fstream fout;
	fout.open(fname.c_str(),fstream::binary|fstream::out|fstream::trunc);

        fout.seekp(0);
        write_Small_Endian(fout, (char*)&MAGIC, sizeof(uint32_t));
        write_Small_Endian(fout, (char*)&VERSION, sizeof(uint16_t));
        write_Small_Endian(fout, (char*)&ft_time, sizeof(uint16_t));
        fout.seekp(24);
        write_Small_Endian(fout, (char*)&data->clock,sizeof(double));
        fout.seekp(1024);

	write_Small_Endian(fout, (char*)(data->data),data->l*sizeof(uint32_t));

	return 0;
}

uint32_t smurf_loadfret2(const string fname, s_fret2 *data){
	fstream fin;
	fin.open(fname.c_str(),fstream::binary|fstream::in);

	uint32_t tm;
	uint16_t tt,tv;

	fin.seekg(0);
	read_Small_Endian(fin, (char*)&tm, sizeof(uint32_t));
	if(tm!=MAGIC) return NOT_SMURF;

	read_Small_Endian(fin, (char*)&tv, sizeof(uint16_t));

	read_Small_Endian(fin, (char*)&tt, sizeof(uint16_t));
	if(tt!=ft_fret) return NOT_FRET;

	fin.seekg(8);
	read_Small_Endian(fin, (char*)&(data->l1),sizeof(uint32_t));
	read_Small_Endian(fin, (char*)&(data->l2),sizeof(uint32_t));

	fin.seekg(24);
	read_Small_Endian(fin, (char*)&(data->clock),sizeof(double));

	fin.seekg(48);
	read_Small_Endian(fin, (char*)&(data->power),sizeof(double));

	fin.seekg(128);
	read_Small_Endian(fin, (char*)&(data->wc1),sizeof(double));
	read_Small_Endian(fin, (char*)&(data->wc2),sizeof(double));

	fin.seekg(128);
	read_Small_Endian(fin, (char*)&(data->ps1),sizeof(double));
	read_Small_Endian(fin, (char*)&(data->ps2),sizeof(double));

	fin.seekg(1024);
	data->data1=new uint32_t[data->l1];
	data->data2=new uint32_t[data->l2];
	read_Small_Endian(fin, (char*)(data->data1),data->l1*sizeof(uint32_t));
	read_Small_Endian(fin, (char*)(data->data2),data->l2*sizeof(uint32_t));

	fin.close();

	return 0;
}

uint32_t smurf_savefret2(const string fname, s_fret2 *data){
	fstream fout;
	fout.open(fname.c_str(),fstream::binary|fstream::out|fstream::trunc);

	fout.seekp(0);
	write_Small_Endian(fout, (char*)&MAGIC, sizeof(uint32_t));
	write_Small_Endian(fout, (char*)&VERSION, sizeof(uint16_t));
	write_Small_Endian(fout, (char*)&ft_fret, sizeof(uint16_t));
	write_Small_Endian(fout, (char*)&data->l1,sizeof(uint32_t));
	write_Small_Endian(fout, (char*)&data->l2,sizeof(uint32_t));

	fout.seekp(24);
        write_Small_Endian(fout, (char*)&data->clock,sizeof(double));
	write_Small_Endian(fout, (char*)&data->power,sizeof(double));

	fout.seekp(128);
	write_Small_Endian(fout, (char*)&data->wc1,sizeof(double));
	write_Small_Endian(fout, (char*)&data->wc2,sizeof(double));

	fout.seekp(160);
	write_Small_Endian(fout, (char*)&data->ps1,sizeof(double));
	write_Small_Endian(fout, (char*)&data->ps2,sizeof(double));

	fout.seekp(1024);
	write_Small_Endian(fout, (char*)data->data1,data->l1*sizeof(uint32_t));
	write_Small_Endian(fout, (char*)data->data2,data->l2*sizeof(uint32_t));

	return 0;
}

double smurf_STD2B(double x, double bd, double ba, double idbT, double iabT){
	double j1=36*pow(x,10)/pow(1+pow(x,6),3);

	double jd=idbT*(1-1/bd)*(1-1/bd)/((pow(x,6)+1/bd));
	double ja=iabT*(1-1/ba)*(1-1/ba)/((1+pow(x,6)/ba));
	double J2B=j1*(ja+jd);
	return pow(J2B,-.5);
}

double smurf_STDEB(double E, double bd, double ba, double IdbT, double IabT){
	double jd=IdbT*pow(1-1/bd,2.0)/(E*(1-1/bd)+1);
	double ja=IabT*pow(1-1/ba,2.0)/(E*(1-1/ba)+1/ba);

	return pow(ja+jd,-0.5);
}

double smurf_MLE2B(double nd, double na, double bd, double ba, double IdbT, double IabT){
	double xorig=pow((IdbT*na*ba-IabT*nd*ba*bd)/(IabT*nd*bd-IdbT*na*ba*bd),1.0/6.0);
	return xorig;
}

double smurf_MLEEB(double nd, double na, double bd, double ba, double IdbT, double IabT){
	double ee=(IdbT*na-IabT*nd/ba)/(IdbT*na*(1-1/bd)+IabT*nd*(1-1/ba));
	if(ee<0 || ee>1) ee=-1;
	return ee;
}

double smurf_STD2B(double x, double T, s_cal b){
	double j1=36*pow(x,10)/pow(1+pow(x,6),3);

	double jd=b.IdB*T*(1-1/b.bd)*(1-1/b.bd)/((pow(x,6)+1/b.bd));
	double ja=b.IaB*T*(1-1/b.ba)*(1-1/b.ba)/((1+pow(x,6)/b.ba));
	double J2B=j1*(ja+jd);
	return pow(J2B,-.5);
}

double smurf_STDEB(double E, double T, s_cal b){
	double jd=b.IdB*T*pow(1-1/b.bd,2.0)/(E*(1-1/b.bd)+1);
	double ja=b.IaB*T*pow(1-1/b.ba,2.0)/(E*(1-1/b.ba)+1/b.ba);

//cout << b.IdB*T*pow(1-1/b.bd,2.0) << ' ' << (E*(1-1/b.bd)-1) << endl;

	return pow(ja+jd,-0.5);
}

double smurf_MLE2B(double nd, double na, s_cal b){
	double xorig=pow((b.IdB*na*b.ba-b.IaB*nd*b.ba*b.bd)/(b.IaB*nd*b.bd-b.IdB*na*b.ba*b.bd),1.0/6.0);
	return xorig;
}

double smurf_MLEEB(double nd, double na, s_cal b){
	double ee=(b.IdB*na-b.IaB*nd/b.ba)/(b.IdB*na*(1-1/b.bd)+b.IaB*nd*(1-1/b.ba));
	if(ee<0 || ee>1) ee=-1;
	return ee;
}

s_cal smurf_process_fret(string infile, fstream& tout) {
	s_cal b;
        
	int ver=smurf_isfret(infile);
        
        fstream fin;
	fin.open(infile.c_str(),ios::binary|ios::in);

        if (ver==3) {
        //There is only acceptor data.    
            b = smurf_process_freta(infile, tout);
            b.Vers = ver;
            return b;
        } else if (ver==4) {
        //There is only donor data.
            b = smurf_process_fretd(infile, tout);
            b.Vers = ver;
            return b;
        }
        
	if(!ver){
		cerr << "Not a fret file!" << endl;
		b.nt=0;
		return b;
	}

	//Find total photon numbers
	uint32_t ld;
	uint32_t la;

	fin.seekg(8);
	read_Small_Endian(fin, (char*)&(ld),sizeof(uint32_t));
	read_Small_Endian(fin, (char*)&(la),sizeof(uint32_t));

	uint32_t lt=la+ld;

        uint32_t ds=1024;
        uint32_t as=1024+ld*sizeof(uint32_t);

	uint32_t nta=0,ntd=0;
	uint64_t cta=0,ctd=0;
	uint64_t tta=0,ttd=0;

	//Clock rate
	uint32_t cr;
	fin.seekg(24);
	read_Small_Endian(fin, (char*)&cr,sizeof(cr));
	if(cr==0) cr=80000000;
	fin.seekg(ds);
	read_Small_Endian(fin, (char*)&ntd,sizeof(uint32_t));
	fin.seekg(as);
	read_Small_Endian(fin, (char*)&nta,sizeof(uint32_t));
	cta=nta;
	ctd=ntd;

	s_fr np;
    
	uint32_t cpcd=0;
	uint32_t cpca=0;
	uint32_t offa=0, offd=0;
	int nza=0,nzd=0;
	uint64_t p232=(uint64_t)pow(2.,32);
	uint64_t p230=(uint64_t)pow(2.,30);

	if(ver==1){
		while((np.cpca<la-1)&&(np.cpcd<ld-1)){
			if(tta+nta<ttd+ntd){
				tta+=nta;
				np.type=ACCEPTOR;
				np.time=tta/double(cr);
				np.cpca++;
				cpca++;
				fin.seekg(as+(cpca+nza)*sizeof(uint32_t));
				read_Small_Endian(fin, (char*)&nta,sizeof(uint32_t));
				if(nta==0){
					np.cpca--;
					nza++;
					continue;
				}
			} else if(ttd+ntd<=tta+nta){
				ttd+=ntd;
				np.type=DONOR;
				np.time=ttd/double(cr);
				np.cpcd++;
				cpcd++;
				fin.seekg(ds+(cpcd+nzd)*sizeof(uint32_t));
				read_Small_Endian(fin, (char*)&ntd,sizeof(uint32_t));
				if(np.time==0){
					np.cpcd--;
					nzd++;
//					cerr << "Donor dt=0" << endl;
					continue;
				}
			}
			write_Small_Endian(tout, (char*)&np,sizeof(np));
		}
	} else if(ver==2){
		while((np.cpca+nza<la-1)&&(np.cpcd+nzd<ld-1)){
			if(cta<ctd){
				tta=cta;
				np.type=ACCEPTOR;
				np.time=tta/double(cr);
				np.cpca++;
				cpca++;
				fin.seekg(as+(cpca+nza)*sizeof(uint32_t));
				read_Small_Endian(fin, (char*)&nta,sizeof(uint32_t));
				if(nta+offa*p232+p230<tta){
					offa++;
				}
				cta=nta+offa*p232;
				if(nta==0){
					np.cpca--;
					nza++;
					continue;
				}
			} else if(ctd<=cta){
				ttd=ctd;
				np.type=DONOR;
				np.time=ttd/double(cr);
				np.cpcd++;
				cpcd++;
				fin.seekg(ds+(cpcd+nzd)*sizeof(uint32_t));
				read_Small_Endian(fin, (char*)&ntd,sizeof(uint32_t));
				if(ntd+offd*p232+p230<ttd) offd++;
				ctd=ntd+offd*p232;
				if(ntd==0){
					np.cpcd--;
					nzd++;
					continue;
				}
			}
			write_Small_Endian(tout, (char*)&np,sizeof(np));
		}
	}

	lt=np.cpca+np.cpcd;

        b=smurf_fret_calibrate(tout,lt,infile);
	b.ntl=lt;
              
        fin.close();
        b.Vers = 2;
	return b;
}

s_cal smurf_process_freta(string infile, fstream& tout) {
        
        fstream fin;
	fin.open(infile.c_str(),ios::binary|ios::in);

	//Find total photon numbers
	uint32_t ld;
	uint32_t la;

	fin.seekg(8);
	read_Small_Endian(fin, (char*)&(ld),sizeof(uint32_t));
	read_Small_Endian(fin, (char*)&(la),sizeof(uint32_t));
        
	uint32_t lt=la;

	uint32_t as=1024;

	uint32_t nta=0;
	uint64_t cta=0;
	uint64_t tta=0;

	//Clock rate
	uint32_t cr;
	fin.seekg(24);
	read_Small_Endian(fin, (char*)&cr,sizeof(cr));
	if(cr==0) cr=80000000;
	fin.seekg(as);
	read_Small_Endian(fin, (char*)&nta,sizeof(uint32_t));
	cta=nta;

	s_fr np;

	uint32_t cpca=0;
	uint32_t offa=0;
	int nza=0;
	uint64_t p232=(uint64_t)pow(2.,32);
	uint64_t p230=(uint64_t)pow(2.,30);

        while(np.cpca+nza<la-1) {
                tta=cta;
                np.type=ACCEPTOR;
                np.time=tta/double(cr);
                np.cpca++;
                cpca++;
                fin.seekg(as+(cpca+nza)*sizeof(uint32_t));
                read_Small_Endian(fin, (char*)&nta,sizeof(uint32_t));
                if(nta+offa*p232+p230<tta){
                    offa++;
                }
                cta=nta+offa*p232;
            if(nta==0){
                np.cpca--;
                nza++;
                continue;
            }
            write_Small_Endian(tout, (char*)&np,sizeof(np));
        }

	lt=np.cpca;

        //Replaced smurf_fret_calibrate with a more specific needed functions......
        
	s_cal b;
        
        b.ta_bleach=tfind(tout,0,lt-1,ACCEPTOR);
        
        fstream goty;
        goty.open(infile.c_str(),ios::binary|ios::in);
        double accBleach = 0;
        goty.seekg(80);
        
        read_Small_Endian(goty, (char*)&accBleach,sizeof(double));
        

        if (accBleach != 0) {
            cout << "Acceptor Bleach from .fr file is " << accBleach << endl;
            printf("resetting acceptor Bleach!!!!\n");
            b.ta_bleach = accBleach;
        }

        goty.close();
        
	uint32_t nla=findtime(tout,lt-1,b.ta_bleach);
	s_fr last_a=read_num(tout,nla);
	s_fr last=read_num(tout,lt-2);
        
        //the intensity in region one is compared to the intensity in region two

	double iIIa=(last_a.cpca)/pow(last_a.time,1);
	double iIIIa=(last.cpca-last_a.cpca)/pow(last.time-last_a.time,1);
        
        b.ba=iIIa/iIIIa;
        b.IaB=iIIa;

        //remains from full implementation for both channels
	double vIIa=pow(sqrt((last_a.cpca)/pow(last_a.time,2))+100,2);
	double vIIIa=(last.cpca-last_a.cpca)/pow(last.time-last_a.time,2);
        
	double B_a=iIIIa;
	double vB_a=vIIIa;
	//double Ia0=iIIa-iIIIa;
	//double vIa0=vIIa+vIIIa;
	double vIaB=vIIa;

	b.vIaB=vIaB;
	b.vba=b.vIaB/pow(B_a,2)+(vB_a)*pow(b.IaB/pow(B_a,2),2);
	b.nt=nla;	
        b.ntl=lt;
        //end of what remains

        fin.close();
	return b;
        
}

s_cal smurf_process_fretd(string infile, fstream& tout) {
        fstream fin;
	fin.open(infile.c_str(),ios::binary|ios::in);

	//Find total photon numbers
	uint32_t ld;
	uint32_t la;

	fin.seekg(8);
	read_Small_Endian(fin, (char*)&(ld),sizeof(uint32_t));
	read_Small_Endian(fin, (char*)&(la),sizeof(uint32_t));

	uint32_t lt=ld;

	uint32_t ds=1024;

	uint32_t ntd=0;
	uint64_t ctd=0;
	uint64_t ttd=0;

	//Clock rate
	uint32_t cr;
	fin.seekg(24);
	read_Small_Endian(fin, (char*)&cr,sizeof(cr));
	if(cr==0) cr=80000000;
	fin.seekg(ds);
	read_Small_Endian(fin, (char*)&ntd,sizeof(uint32_t));
	ctd=ntd;

	s_fr np;
    
	uint32_t cpcd=0;
	uint32_t offd=0;
	int nzd=0;
	uint64_t p232=(uint64_t)pow(2.,32);
	uint64_t p230=(uint64_t)pow(2.,30);

        while(np.cpcd+nzd<ld-1) {
                ttd=ctd;
                np.type=DONOR;
                np.time=ttd/double(cr);
                np.cpcd++;
                cpcd++;
                fin.seekg(ds+(cpcd+nzd)*sizeof(uint32_t));
                read_Small_Endian(fin, (char*)&ntd,sizeof(uint32_t));
                if(ntd+offd*p232+p230<ttd) offd++;
                ctd=ntd+offd*p232;
                if(ntd==0){
                        np.cpcd--;
                        nzd++;
                        continue;
                }
                write_Small_Endian(tout, (char*)&np,sizeof(np));
        }

	lt=np.cpcd;

        //Replaced smurf_fret_calibrate with a more specific needed functions......
	s_cal b;
        
        b.td_bleach=tfind(tout,0,lt-1,DONOR);
        
        fstream goty;
        goty.open(infile.c_str(),ios::binary|ios::in);

        double donBleach = 0;
	    goty.seekg(72);
        read_Small_Endian(goty, (char*)&donBleach,sizeof(double));

        if (donBleach != 0) {
            cout << "Donor Bleach from .fr file is " << donBleach << endl;
            printf("resetting donorBleach!\n");
            b.td_bleach = donBleach;
        }

        goty.close();
        
	uint32_t nld=findtime(tout,lt-1,b.td_bleach);
	s_fr last_d=read_num(tout,nld);
	s_fr last=read_num(tout,lt-2);

        //the intensity in region one is compared to the intensity in region two

	double iIId=(last_d.cpcd)/pow(last_d.time,1);
	double iIIId=(last.cpcd-last_d.cpcd)/pow(last.time-last_d.time,1);
        
        b.bd=iIId/iIIId;
        b.IdB=iIId;
        
        //remainds from full implementation
        double vIId=pow(sqrt((last_d.cpcd)/pow(last_d.time,2))+100,2);
	double vIIId=(last.cpcd-last_d.cpcd)/pow(last.time-last_d.time,2);
        
	double B_d=iIIId;
	double vB_d=vIIId;
	//double Id0=iIId-iIIId;
	//double vId0=vIId+vIIId;
	double vIdB=vIId;
        
	b.vIdB=vIdB;
	b.vbd=b.vIdB/pow(B_d,2)+(vB_d)*pow(b.IdB/pow(B_d,2),2);
	b.nt=nld;
        b.ntl=lt;
        //end of what remains
        
        
        fin.close();
	return b;
}

double vdet(fstream& tin,uint32_t lt,s_fr beg,s_fr end,int n){
	double T=end.time-beg.time;
	double dT=T/n;
	double I;
	s_fr lp,np;
	lp=beg;
	double xx=0,x=0;

	//Calculate binned intensity trajectory and variance
	for(int i=0;i<n;i++){
		np=read_num(tin,findtime(tin,lt-1,beg.time+(i+1)*dT));
		I=(np.cpcd-lp.cpcd)/(np.time-lp.time);
		xx+=I*I;
		x+=I;
		lp=np;
	}

	return xx/n-x/n*x/n;
}

xt_type smurf_fret_xt_parse(char* xtv){
	xt_type sxt;
	
	switch(atoi(xtv)){
	case 0:
		sxt=NONE;
		break;
	case 1:
		sxt=LPW;
		break;
	case 2:
		sxt=JANDM;
		break;
	case 3:
		sxt=XUN;
		break;
	case 4:
		sxt=CY35;
		break;	
	default:
		sxt=UNSET;
		break;
	}
	return sxt;
}


s_cal smurf_fret_calibrate(fstream& tin, uint32_t lt, string infile) {
	//uint32_t ld;
	//uint32_t la;
	s_cal b;
	double Xd,Xa;
	
	//This line reads in the environmental variable
 	//if(SMURF_XT==UNSET) SMURF_XT=smurf_fret_xt_parse(getenv("SMURF_XT"));

    fstream goty;
    goty.open(infile.c_str(),ios::binary|ios::in);

	goty.seekg(89);
	double microUsed = 0;
    read_Small_Endian(goty, (char*)&microUsed,sizeof(double));
	
	//cout << "The double address is: " << microUsed << endl;
	if (SMURF_XT==UNSET) {
		cout << "Automatically Setting Crosstalk for microsope: 192.168.1." << microUsed << endl;
		if (microUsed == 55) {
			SMURF_XT = LPW;			//microscope 1 - berkeley
		} else if (microUsed == 71) {
			SMURF_XT = JANDM;		//m2 - berkeley
		} else if (microUsed == 130) {
			SMURF_XT = JANDM;		//m3 - berkeley 
		} else if (microUsed == 111) {
			SMURF_XT = LPW;			//m1 - princeton
		} else if (microUsed == 121) {
			SMURF_XT = JANDM;		//m2 - princeton
		} else if (microUsed == 131) {
			SMURF_XT = JANDM;		//m3 - princeton
		} else if (microUsed == 135) {
			SMURF_XT = XUN;		//m4 princeton
		//	cout << "PLEASE CORRECT THE X-TALK VALUES FOR M4!!!" << endl;




		} else {
			cout << "The crosstalk values aren't specified for this microsope...." << endl;
			SMURF_XT = NONE;
		}
	} else {
		cout << "Crosstalk supplied by user." << endl;
	}	

	switch(SMURF_XT){
	case NONE:
		cerr << "Crosstalk: None" << endl;
		Xd=0;
		Xa=0;
		break;
	case LPW:	
		cerr << "Crosstalk: Karl" << endl;
		Xa=0.019;
		Xd=0.142;
		break;
	case JANDM:	
		cerr << "Crosstalk: Jeff" << endl;
		Xa=1.85e-4;
		Xd=0.141;
		break;
        case XUN:
               	cerr << "Crosstalk: Xun" << endl;
                Xa=1.5011e-6;
                Xd=0.0738;
                break;
	case CY35:
		cerr << "Crosstalk: Cy3/Cy5" << endl;
		Xa=0.0351;
		Xd=0.0672;
		break;	
	default:
		cerr << "No Crosstalk specified!" << endl;
		exit(1);
		break;
	}
	
        double accBleach = 0;
        double donBleach = 0;
        b.td_bleach=tfind(tin,0,lt-1,DONOR);
	    goty.seekg(72);
        read_Small_Endian(goty, (char*)&donBleach,sizeof(double));
        if (donBleach != 0) {
            printf("resetting donorBleach!\n");
            b.td_bleach = donBleach;
        }

        uint32_t nld=findtime(tin,lt-1,b.td_bleach);
	s_fr last_d=read_num(tin,nld);
	b.ta_bleach=tfind(tin,0,lt-1,ACCEPTOR);
	goty.seekg(80);
        read_Small_Endian(goty, (char*)&accBleach,sizeof(double));
        if (accBleach != 0) {
            printf("resetting acceptorBleach!\n");
            b.ta_bleach = accBleach;
        }
        if((accBleach != 0) || (donBleach != 0)) { 
		   cout << accBleach << '\t' << donBleach << endl;
	    }
		uint32_t nla=findtime(tin,lt-1,b.ta_bleach);
	s_fr last_a=read_num(tin,nla);

        goty.close();

	s_fr last=read_num(tin,lt-2);

	double iId=last_a.cpcd/last_a.time;
	double vId=last_a.cpcd/pow(last_a.time,2);
	double iIa=last_a.cpca/last_a.time;
	double vIa=last_a.cpca/pow(last_a.time,2);
	double iIId=(last_d.cpcd-last_a.cpcd)/pow(last_d.time-last_a.time,1);
	double vIId=(last_d.cpcd-last_a.cpcd)/pow(last_d.time-last_a.time,2);
	//double vIId=pow(sqrt((last_d.cpcd-last_a.cpcd)/pow(last_d.time-last_a.time,2))+100,2);
	double iIIa=(last_d.cpca-last_a.cpca)/pow(last_d.time-last_a.time,1);
	double vIIa=(last_d.cpca-last_a.cpca)/pow(last_d.time-last_a.time,2);
	double iIIId=(last.cpcd-last_d.cpcd)/pow(last.time-last_d.time,1);
	double vIIId=(last.cpcd-last_d.cpcd)/pow(last.time-last_d.time,2);
	double iIIIa=(last.cpca-last_d.cpca)/pow(last.time-last_d.time,1);
	double vIIIa=(last.cpca-last_d.cpca)/pow(last.time-last_d.time,2);
	double B_d=iIIId;
	double vB_d=vIIId;
	double B_a=iIIIa;
	double vB_a=vIIIa;
	double Id0=0;
	double vId0=0;
	double Ia0=0;
	double vIa0=0;

	if(getenv("SMURF_RATIO")==NULL){
		Id0=iIId-iIIId;
		vId0=vIId+vIIId;
		//	double Ia01=iId*(iIIa-iIIIa) + iIId*(iIIIa-iIa) + iIIId*(iIa-iIIa);
		//	double Ia02=(iId-iIId-iIa*Xa+iIIa*Xa);
		//	double Ia0=Ia01/Ia02;
		//	double Xd=(iIIa-iIIIa)/Id0;
		//cerr << "Xd = " << Xd << endl;
		double P=(iIId-iId)/(iIa-iIIa);
		//	double vP=(vIId+vId)
		Ia0=(Id0+Xd*Id0*P)/(P+Xa);
		double vIa01=pow(iId-iIId+Xa*(iIIa-iIa),4);
		double vIa02=vIIId*pow(iId - iIId + (-iIa + iIIa)*Xa,2)*pow(iIa - iIIa + (-iId + iIId)*Xd,2) + pow(iId - iIId,2)*pow(iIId - iIIId,2)*vIa*pow(-1 +Xa*Xd,2) + pow(iIa - iIIa,2)*pow(iIId - iIIId,2)*vId*pow(-1 + Xa*Xd,2) +pow(iId - iIId,2)*pow(iIId - iIIId,2)*vIIa*pow(-1 + Xa*Xd,2) +vIId*pow((iIa - iIIa)*(-iId + iIIId + (iIa - iIIa)* Xa) + (pow(iId -iIId,2) - (iIa - iIIa)*(iId - 2*iIId +iIIId)*Xa)*Xd,2);
		vIa0=vIa02/vIa01;
	} else{
		double ratio=atof(getenv("SMURF_RATIO"));
		cout << "Using SMURF_RATIO=" << ratio << endl;
		Id0= iId-iIIId-Xa*(iIa-iIIIa)+(iIa-iIIIa-Xd*(iId-iIIId))/ratio;
		vId0=vId+vIIId+Xa*(vIa+vIIIa)+(vIa+vIIIa+Xd*(vId+vIIId))/ratio;
		Ia0=Id0*ratio;
		vIa0=vId0*ratio;
	}
	
	double vIdB=vId0+vIIId;
	double vIaB=vIa0+vIIIa;
		
//	double vIa0=vIId*pow(,2)
//	double Ebar=(iIa-B_a-Xd*Id0)/(Ia0-Xd*Id0);
//	double Xa=(iId-B_d-Id0*(1-Ebar))/(Ia0*Ebar);
//cerr << "Xa = " << Xa << endl;
//	b.P=P;
	b.IdB=Id0+B_d;
	b.vIdB=vIdB;
	b.IaB=Ia0+B_a;
	b.vIaB=vIaB;
	b.bd=b.IdB/(B_d+Xa*Ia0);
	b.vbd=b.vIdB/pow(B_d+Xa*Ia0,2)+(vB_d+Xa*vIa0)*pow(b.IdB/pow(B_d+Xa*Ia0,2),2);
	b.ba=b.IaB/(B_a+Xd*Id0);
	b.vba=b.vIaB/pow(B_a+Xa*Ia0,2)+(vB_d+Xd*vId0)*pow(b.IaB/pow(B_a+Xd*Id0,2),2);
	b.nt=nla;

	double vbd=vIdB/pow(B_d+Xa*Ia0,2)+vB_d*pow(-b.IdB/pow(B_d+Xa*Ia0,2),2)+vIa0*pow(Xa*b.IdB/pow(B_d+Xa*Ia0,2),2);
	double vba=vIaB/pow(B_a+Xd*Id0,2)+vB_a*pow(-b.IaB/pow(B_a+Xd*Id0,2),2)+vId0*pow(Xd*b.IaB/pow(B_a+Xd*Id0,2),2);

	vIdB=2.43e4;
	vIaB=1.06e5;
	vbd=.2;
	vba=.4;

// double vxa=36*b.ba*pow(-b.bd*b.IaB+b.ba*b.IdB*last_a.cpca,3.0)*last_a.cpcd*(b.IdB*last_a.cpca-b.bd*b.IaB*last_a.cpcd)*pow(b.ba*(-b.IdB*last_a.cpca+b.bd*b.IaB*last_a.cpcd)/(-b.bd*b.IaB*last_a.cpcd+b.ba*b.IdB*last_a.cpca*last_a.cpcd),2.0/3.0);
// double vxb=2*pow(b.bd*b.IaB,3.0)*b.IdB*last_a.cpca*last_a.cpcd*vba-pow(b.bd*b.IaB,4.0)*pow(last_a.cpcd,2.0)*vba-pow(b.ba*b.IaB*b.IdB*last_a.cpca*(b.ba*last_a.cpcd-1),2.0)*vbd;
// double vxc=-pow(b.bd*last_a.cpca,2.0)*(pow(b.ba*b.IdB*(b.ba*last_a.cpcd-1),2.0)*vIaB+pow(b.IaB,2.0)*(pow(b.IdB,2.0)*vba+pow(b.ba*(b.ba*last_a.cpcd-1),2.0)*vIdB));
//cerr << sqrt(vbd) << ' ' << sqrt(vba) << ' ' << sqrt((vxb+vxc)/vxa) << endl;
	return b;
}

//Determines calibration values for the FRET trajectory in a processed data file
s_cal smurf_fret_cross(fstream& tin, uint32_t lt,double Ebar){
	s_cal b;

	b.td_bleach=tfind(tin,0,lt-1,DONOR);
	uint32_t nld=findtime(tin,lt-1,b.td_bleach);
	s_fr last_d=read_num(tin,nld);
	b.ta_bleach=tfind(tin,0,lt-1,ACCEPTOR);
	uint32_t nla=findtime(tin,lt-1,b.ta_bleach);
	s_fr last_a=read_num(tin,nla);

	s_fr last=read_num(tin,lt-2);

// 	double iId=last_a.cpcd/last_a.time;
// 	double iIa=last_a.cpca/last_a.time;
// 	double iIId=(last_d.cpcd-last_a.cpcd)/(last_d.time-last_a.time);
// 	double iIIa=(last_d.cpca-last_a.cpca)/(last_d.time-last_a.time);
// 	double iIIId=(last.cpcd-last_d.cpcd)/(last.time-last_d.time);
// 	double iIIIa=(last.cpca-last_d.cpca)/(last.time-last_d.time);
// 	double B_a=iIIIa;
// 	double B_d=iIIId;
// 	double Id0=iIId-iIIId;
	
//	double Xa=(iId-iIId+Ebar*iIId-Ebar*iIIId)/(iIa-iIIa+Ebar*iIIa-Ebar*iIIIa);
//	double Xd=(iIIa-iIIIa)/(iIId-iIIId);
//	double Ia0=(iIa-iIIa+Ebar*iIIa-Ebar*iIIIa)/Ebar;

//cerr << "Xd = " << Xd << endl;
//cerr << "Xa = " << Xa << endl;

//	b.P=P;
//	b.IdB=Id0+B_d;
//	b.IaB=Ia0+B_a;
//	b.bd=b.IdB/(B_d+Xa*Ia0);
//	b.ba=b.IaB/(B_a+Xd*Id0);
//	b.nt=nla;
cerr << "!!!!!!!!!!!" << endl;
	return b;
}

//Returns the index of the photon closest to time t. If before is true, returns
//the photon previous to t. If before is false, returns photon subsequent to t.
uint32_t findtime(fstream& f, uint32_t N, double t,bool before){
	if(t==0) return 0;

	uint32_t na=0;
	uint32_t nb=N-1;
	uint32_t nn=(uint32_t)(PHI*N);

	s_fr pa=read_num(f,na);
	s_fr pb=read_num(f,nb);
	s_fr pn=read_num(f,nn);

	while(nb-na>1){
		if(pn.time>t){
			nb=nn;
			pb=pn;
		} else if(pn.time<t){
			na=nn;
			pa=pn;
		} else{
			na=nn;
			pa=pn;
			nb=nn;
			pb=pn;
		}

		nn=(uint32_t)((1-PHI)*na+PHI*nb);
		pn=read_num(f,nn);
	}

	if(before)
		return na;
	else
		return nb;
}

//Returns the photon with index n
s_fr read_num(fstream& f, long n){
	s_fr r;
	if(n<0) return r;
	f.seekg(n*sizeof(r));
	read_Small_Endian(f, (char*)&r,sizeof(r));
	return r;
}

//shifting routine from numerical recipes in C
inline void S2(double &a,double &b,double c){
	a=b;
	b=c;
}

//shifting routine from numerical recipes in C
inline void S3(uint32_t &a,uint32_t &b,uint32_t &c,uint32_t d){
	a=b;
	b=c;
	c=d;
}

//Find the most likely intensity change point in the given processed
//data file, in the segment that begins with photon ns and is n photons long.
double tfind(fstream &data,uint32_t ns,uint32_t n,int chan){

	uint32_t n0,n1,n2,n3,nm,nmt;
	double llmt,llm;
	int N=4;

	for(int i=0;i<N;i++){
		n0=ns+i*n/N;
		n3=ns+(i+1)*n/N;
		n1=2*n0/3+n3/3;
		n2=n0/3+2*n3/3;

		double ll1=-LLT(data,ns,n,n1,chan);
		double ll2=-LLT(data,ns,n,n2,chan);

		while(n3-n0>4){
			if(ll2<ll1){
				S3(n0,n1,n2,(uint32_t)(PHI*n3+(1-PHI)*n2));
				S2(ll1,ll2,-LLT(data,ns,n,n2,chan));
			} else {
				S3(n3,n2,n1,(uint32_t)(PHI*n0+(1-PHI)*n1));
				S2(ll2,ll1,-LLT(data,ns,n,n1,chan));
			}
		}

		if(ll1<ll2){
			nmt=n1;
			llmt=ll1;
		}else{
			nmt=n2;
			llmt=ll2;
		}

		if((i==0)||(llmt<llm)){
			nm=nmt;
			llm=llmt;
		}
	}

	s_fr cp=read_num(data,nm);
	return cp.time;
}

//Calculate the log-likelihood ratio for there being a change point at the nt-th photon
//in the data segment that starts at ns and is n photons long.
double LLT(fstream &tin, uint32_t ns, uint32_t n, uint32_t nt,int chan){
	s_fr ps=read_num(tin,ns);
	s_fr pm=read_num(tin,ns+nt);
	s_fr pl=read_num(tin,ns+n-5);

	double l0=(pl.cpc(chan)-ps.cpc(chan))/(pl.time-ps.time);
	double l1=(pm.cpc(chan)-ps.cpc(chan))/(pm.time-ps.time);
	double l2=(pl.cpc(chan)-pm.cpc(chan))/(pl.time-pm.time);

	double LL=(pm.cpc(chan)-ps.cpc(chan))*log(l1/l0)+(pl.cpc(chan)-pm.cpc(chan))*log(l2/l0);

	return LL;
}

void spinner(uint32_t ndone, uint32_t ntot, int nskip){
	static char* spin="-\\|/";
	static int perc=0;
	static int spinint=0;

	if(ndone==0){
		cerr << "Progress:  0% -";
		spinint++;
	}

	if(ndone%nskip==0) cerr << '\b' << spin[(1000*ndone/ntot)%4];

	if((ndone)%(ntot/100+1)==0){
		perc=((ndone+1)*100)/ntot;

		cerr << "\b\b\b\b\b\b" << setw(3) << perc << "% " << spin[(1000*ndone/ntot)%4];

		cerr.flush();
    }

	if(ndone==ntot-1){
		cerr << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b";
	}

}

void smurf_getlongopt(){
	//	struct option oo[5];
}

//These are some utility functions added to help in debugging hist.cpp
//Now they are also used in hist.cpp and boots.cpp to print gsl_matrix types to files and terminal...

void PrintMatrix(gsl_matrix* toPrint) {
	for (uint i=0;i<(*toPrint).size1;i++) {
		for (uint p=0;p<(*toPrint).size2;p++) {
			cout << gsl_matrix_get(toPrint,i,p) << '\t';
		}
	    cout << endl;
	}
}

void PrintMatrix(gsl_matrix* toPrint,string filename) {
	ofstream Outfile(filename.c_str());
	for (uint i=0;i<(*toPrint).size1;i++) {
		for (uint p=0;p<(*toPrint).size2;p++) {
			Outfile << gsl_matrix_get(toPrint,i,p) << '\t';
		}
	    Outfile << endl;
	}
}

//This function uses a CLAPACK function to solve for the eigen values and vectors A x = lambda B x for A,B matrices...
//Modeled after eig function from matlab in input/output format... 

extern "C" {


#include "f2c.h"
#include "clapack.h"

void eig(gsl_matrix** W, gsl_matrix** d, gsl_matrix* g, gsl_matrix* M) {
	integer N = 3;
	*W = gsl_matrix_calloc(3,3);
	*d = gsl_matrix_calloc(3,3);
	double *A;
	double *B;
    A = (double*) malloc(N*N*sizeof(double));
	B = (double*) malloc(N*N*sizeof(double));
   
	//first for A assignment...
	A[0] = gsl_matrix_get(g,0,0);
	A[3] = gsl_matrix_get(g,0,1);
	A[6] = gsl_matrix_get(g,0,2);
	A[1] = gsl_matrix_get(g,1,0);
	A[4] = gsl_matrix_get(g,1,1);
	A[7] = gsl_matrix_get(g,1,2);
	A[2] = gsl_matrix_get(g,2,0);
	A[5] = gsl_matrix_get(g,2,1);
	A[8] = gsl_matrix_get(g,2,2);
    
	//Now for B assignment...
	B[0] = gsl_matrix_get(M,0,0);
	B[3] = gsl_matrix_get(M,0,1);
	B[6] = gsl_matrix_get(M,0,2);
	B[1] = gsl_matrix_get(M,1,0);
	B[4] = gsl_matrix_get(M,1,1);
	B[7] = gsl_matrix_get(M,1,2);
	B[2] = gsl_matrix_get(M,2,0);
	B[5] = gsl_matrix_get(M,2,1);
	B[8] = gsl_matrix_get(M,2,2);


	/*
	A[0] = gsl_matrix_get(g,0,0);
	A[2] = gsl_matrix_get(g,0,1);
	A[1] = gsl_matrix_get(g,1,0);
	A[3] = gsl_matrix_get(g,1,1);

	B[0] = gsl_matrix_get(M,0,0);
	B[2] = gsl_matrix_get(M,0,1);
	B[1] = gsl_matrix_get(M,1,0);
	B[3] = gsl_matrix_get(M,1,1);
	*/

    //means return eigen values and eigen vectors...
    double *alphar = (double*) malloc(N*sizeof(double));
    double *alphai = (double*) malloc(N*sizeof(double));
	double *beta = (double*) malloc(N*sizeof(double));
	double *eigvectVL = (double*) malloc(N*N*sizeof(double));
	//Should not be used because JOBVR = 'N' ... but need to send function something...
	double *eigvectVR = (double*) malloc(N*N*sizeof(double));
	integer ldvl = N;
	integer ldvr = N;
	char JOBVL = 'N';
    char JOBVR = 'V';
    integer info;
    integer lda = N;
    integer ldb = N;
	//117 was returned as w[0] for A and B 3 by 3...
	integer worklength = (117)*sizeof(double);
	double *w = (double*) malloc(worklength);
    
    dggev_(&JOBVL,&JOBVR,&N,A,&lda,B,&ldb,alphar,alphai,beta,eigvectVL,&ldvl,eigvectVR,&ldvr,w,&worklength,&info); 

	//	printf("w = %f\n",w[0]);

	//	printf("DGGEV INFO (0 if function exited normally): %d\n", info);
	
	/*
	gsl_matrix_set(*W,0,0,eigvectVR[0]);
	gsl_matrix_set(*W,0,1,eigvectVR[2]);
	gsl_matrix_set(*W,1,0,eigvectVR[1]);
	gsl_matrix_set(*W,1,1,eigvectVR[3]);

	gsl_matrix_set(*d,0,0,(alphar[0])/beta[0]);
	gsl_matrix_set(*d,0,1,0);
	gsl_matrix_set(*d,1,0,0);
	gsl_matrix_set(*d,1,1,(alphar[1]+alphai[1])/beta[1]);
    */

    gsl_matrix_set(*W,0,0,eigvectVR[0]);
    gsl_matrix_set(*W,0,1,eigvectVR[3]);
    gsl_matrix_set(*W,0,2,eigvectVR[6]);
    gsl_matrix_set(*W,1,0,eigvectVR[1]);
    gsl_matrix_set(*W,1,1,eigvectVR[4]);
    gsl_matrix_set(*W,1,2,eigvectVR[7]);
    gsl_matrix_set(*W,2,0,eigvectVR[2]);
    gsl_matrix_set(*W,2,1,eigvectVR[5]);
    gsl_matrix_set(*W,2,2,eigvectVR[8]);
    //Now for eigvals...
    gsl_matrix_set(*d,0,0,(alphar[0])/beta[0]);
    gsl_matrix_set(*d,0,1,0.0);
    gsl_matrix_set(*d,0,2,0.0);
    gsl_matrix_set(*d,1,0,0.0);
    gsl_matrix_set(*d,1,1,(alphar[1]+alphai[1])/beta[1]);
    gsl_matrix_set(*d,1,2,0.0);
    gsl_matrix_set(*d,2,0,0.0);
    gsl_matrix_set(*d,2,1,0.0);
	gsl_matrix_set(*d,2,2,(alphar[2]+(2)*alphai[2])/beta[2]);

	free(A);
	A = NULL;
	free(B);
	B = NULL;
	free(alphar);
	alphar = NULL;
	free(alphai);
	alphai = NULL;
	free(beta);
	beta = NULL;
	free(eigvectVR);
	eigvectVR = NULL;
	free(eigvectVL);
	eigvectVL = NULL;
    free(w);
}
}
