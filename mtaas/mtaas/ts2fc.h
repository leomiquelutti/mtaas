#ifndef TS_2_FC
#define TS_2_FC

#include "stdafx.h"

class SetUpRemoteReferenceConcomitance {
public:

	SetUpRemoteReferenceConcomitance() {};
	~SetUpRemoteReferenceConcomitance() {};	

	SetUpRemoteReferenceConcomitance(const SetUpRemoteReferenceConcomitance& element) {*this = element;};
    SetUpRemoteReferenceConcomitance& operator = (const SetUpRemoteReferenceConcomitance& element);

	static void find_concomitance_for_stations( std::vector<StationBase> &station );
};

class Extract_FCs_FixedWindowLength {
public:

	Extract_FCs_FixedWindowLength( int N, int winLngth = 4096 ):
		factorOfEachDecimationLevel(4),
		overlap(0.25),
		minNumDataPerFreq(20),
		pLRatio(0.05) {
			nPointsTS = N;
			windowLength = winLngth;
			this->set_parameters(); };

	~Extract_FCs_FixedWindowLength() {};	

	Extract_FCs_FixedWindowLength(const Extract_FCs_FixedWindowLength& element) {*this = element;};
    Extract_FCs_FixedWindowLength& operator = (const Extract_FCs_FixedWindowLength& element);

	matCUDA::Array<size_t> *fcDistribution;

	static void get_FCs( std::vector<StationBase> &station );

private:

	// members
	size_t windowLength;
	size_t nDecimationLevel;
	size_t factorOfEachDecimationLevel;
	double overlap;
	size_t minNumSegments;
	size_t minNumDataPerFreq;
	double pLRatio;
	size_t nFreq;
	size_t npw;
	size_t nPointsTS;

	// function
	void set_parameters();
};

class Extract_FCs_VariableWindowLength {
public:
	Extract_FCs_VariableWindowLength( int N ) {};
	~Extract_FCs_VariableWindowLength() {};

	Extract_FCs_VariableWindowLength(const Extract_FCs_VariableWindowLength& element) {*this = element;};
    Extract_FCs_VariableWindowLength& operator = (const Extract_FCs_VariableWindowLength& element);

	static void get_FCs( std::vector<StationBase> &station );
};

#endif