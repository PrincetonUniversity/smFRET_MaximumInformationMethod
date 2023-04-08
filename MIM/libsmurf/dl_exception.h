#ifndef DLE_H
#define DLE_H

class dle{
public:
	enum EType{
		NOT_READY,
		LTRAJ_LLP_SIZE,
		LTRAJ_PHOTON_TYPE,
		LTRAJ_NO_SUCH_MODEL,
		XTRAJ_NAN,
		XTRAJ_LENGTH_MISMATCH,
		XTRAJ_RANGE,
		XTRAJ_NO_DATA,
		LANDSCAPE_RANGE,
		LANDSCAPE_NO_DATA,
		LANDSCAPE_LENGTH_MISMATCH,
		INVALID_ORDER,
		LTRAJ_PARAM_LENGTH_MISMATCH
	};

	dle(EType nt){t=nt;}
	~dle(){}
	int type(){return t;}

private:
	EType t;
};

#endif
