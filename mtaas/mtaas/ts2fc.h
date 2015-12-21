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

	Extract_FCs_FixedWindowLength( bool acquisitionContinuous, int N, int winLngth = 4096 );

	~Extract_FCs_FixedWindowLength() {};	

	Extract_FCs_FixedWindowLength(const Extract_FCs_FixedWindowLength& element) {*this = element;};
    Extract_FCs_FixedWindowLength& operator = (const Extract_FCs_FixedWindowLength& element);

	matCUDA::Array<int> *fcDistribution;

	static void get_all_FCs( std::vector<StationBase> &station );
	void set_corrections( StationFile &ts );

private:

	// members
	size_t windowLength;
	size_t fcDraftMadeForWinLengthOf;
	size_t nDecimationLevel;
	size_t factorOfEachDecimationLevel;
	double overlap;
	size_t minNumDataPerFreq;
	double pLRatio;
	size_t nFreq;
	size_t nArCoeff;
	size_t nPointsTS;
	bool isAcquisitionContinuous;
	size_t firstsLinesToRepeatFromFcDistribution;

	// function
	size_t determine_decimation_level( matCUDA::Array<int> auxFcDistribution );
	matCUDA::Array<int> get_draft_of_fc_distribution();
	void set_parameters();
	matCUDA::Array<int> get_fc_distribution( matCUDA::Array<int> auxFcDistribution );
	//void set_corrections( std::vector<Channel> &ch, double samplingFrequency );
};

class Extract_FCs_VariableWindowLength {
public:
	Extract_FCs_VariableWindowLength( int N ) {};
	~Extract_FCs_VariableWindowLength() {};

	Extract_FCs_VariableWindowLength(const Extract_FCs_VariableWindowLength& element) {*this = element;};
    Extract_FCs_VariableWindowLength& operator = (const Extract_FCs_VariableWindowLength& element);

	static void get_all_FCs( std::vector<StationBase> &station );
};

#endif