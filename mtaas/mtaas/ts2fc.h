#ifndef TS_2_FC
#define TS_2_FC

#include "stdafx.h"

class Extract_FCs_FixedWindowLength {
public:

	Extract_FCs_FixedWindowLength( std::vector<StationBase> &station, int winLngth = 4096 );

	~Extract_FCs_FixedWindowLength() {};	

	Extract_FCs_FixedWindowLength(const Extract_FCs_FixedWindowLength& element) {*this = element;};
    Extract_FCs_FixedWindowLength& operator = (const Extract_FCs_FixedWindowLength& element);

private:

	matCUDA::Array<int> *fcDistribution;

	void get_all_FCs( std::vector<StationBase> &station );

	void set_corrections( std::vector<StationBase> &station );
	void extract_fcs_for_each_combination( std::vector<StationBase> &station, const size_t idxStn, const size_t idxTs, const size_t idxCombination );
	void correct_fcs( matCUDA::Array<ComplexDouble> *data, matCUDA::Array<ComplexDouble> *correction, size_t deciLevel );
	void allocate_memory_for_ch_fc( std::vector<StationBase> &station );
	void deallocate_memory_for_ch_fc( std::vector<StationBase> &station );
	void allocate_memory_for_FrequencyResponses( std::vector<StationBase> &station );

	//void get_in( std::vector<StationBase> &station, const size_t iStnIn, const size_t iTsIn, const size_t iComb, const size_t idxSeg  );
	//void get_out( std::vector<StationBase> &station, const size_t iStnIn, const size_t iTsIn, const size_t iComb, const size_t idxSeg  );
	//void get_in_out_inRR( std::vector<StationBase> &station, const size_t iStnIn, const size_t iTsIn, const size_t iComb, const size_t iStnOut, const size_t iTsOut, const size_t idxSeg  );
	void get_in_out_inRR( std::vector<StationBase> &station, const size_t istn, const size_t its, const size_t icomb, const size_t iseg  );
	
//private:

	// members
	matCUDA::Array<double>	*auxDpss;
	matCUDA::Array<int>		*auxFcDistribution;
	size_t					dpssNW;
	size_t					dpssDegree;
	size_t					factorOfEachDecimationLevel;
	size_t					fcDraftMadeForWinLengthOf;
	size_t					firstsLinesToRepeatFromFcDistribution;
	bool					isAcquisitionContinuous;
	size_t 					minNumDataPerFreq;
	size_t 					nArCoeff;
	size_t 					nFreq;
	size_t 					nPointsTS;
	double 					overlap;
	double 					pLRatio;
	size_t 					windowLength;

	// function
	size_t determine_decimation_level( size_t nPoints );
	void get_draft_of_fc_distribution();
	void set_parameters();
	matCUDA::Array<int> get_fc_distribution( size_t maxDecimationLevel );
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