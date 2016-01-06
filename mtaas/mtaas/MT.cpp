#include "MT.h"
#include "Mtu.h"
#include "ts2fc.h"

class Extract_FCs_FixedWindowLength;
class Extract_FCs_VariableWindowLength;
//#include <sstream>

StationBase& StationBase::operator = (const StationBase& element)
{
	this->fileType = element.fileType;
	this->path = element.path;
	this->position = element.position;
	this->stationName = element.stationName;
	this->ts = element.ts;

	return *this;
}

StationFile& StationFile::operator = (const StationFile& element)
{
	this->date = element.date;
	this->exDipoleLength = element.exDipoleLength;
	this->eyDipoleLength = element.eyDipoleLength;
	this->isContinuous = element.isContinuous;
	this->maxFreqInHertz = element.maxFreqInHertz;
	this->minFreqInHertz = element.minFreqInHertz;
	this->mtu = element.mtu;
	this->samplingFrequency = element.samplingFrequency;
	this->timeVector = element.timeVector;
	this->ch = element.ch;
	this->isThereRR = element.isThereRR;
	this->numberOfCombinations = element.numberOfCombinations;
	this->isAcquisitionContinuous = element.isAcquisitionContinuous;

	return *this;
}

Channel& Channel::operator = (const Channel& element)
{
	this->arCoeff = element.arCoeff;
	this->timeSeries = element.timeSeries;
	this->systemResponse = element.systemResponse;
	this->systemResponseFreqs = element.systemResponseFreqs;
	this->countConversion = element.countConversion;
	this->orientationHor = element.orientationHor;
	this->orientationVer = element.orientationVer;
	this->name = element.name;
	this->gain = element.gain;

	return *this;
}

FrequencyResponses& FrequencyResponses::operator = (const FrequencyResponses& element)
{
	this->frequency = element.frequency;
	this->in = element.in;
	this->out = element.out;

	return *this;
}

Combination& Combination::operator = (const Combination& element)
{
	this->idxStn = element.idxStn;
	this->idxTs = element.idxTs;

	return *this;
}

std::string Date::getDateStr()
{
	return std::to_string(startYear) + "-"
		+ std::to_string(startMonth) + "-"
		+ std::to_string(startDay) + " "
		+ std::to_string(startHour) + ":"
		+ std::to_string(startMinute) + ":"
		+ std::to_string(startSecond);
}

void StationBase::read_time_series( DirectoryProperties *dirInfo )
{
	switch(fileType) 
	{
	case FILE_TYPE_MTU:
		ExtractorMTU::read_time_series( this, dirInfo );
	//case FILE_TYPE_ADU_07:
		//ExtractorADU::read_time_series( this, dirInfo );
	}
}

void StationBase::get_all_FCs( std::vector<StationBase> &station, TS_TO_FFT_TYPE ts2fft_type, RR_OR_SS_TYPE rrorss_type )
{
	switch(rrorss_type) 
	{
	case REMOTE_REFERENCE:
		// responsible for define all possible combinations and setup where to read from
		Combination::define_combinations_for_rr( station );
	case SINGLE_SITE:
		Combination::define_combinations_for_ss( station );
	}

	switch(ts2fft_type) 
	{
	case FIXED_WINDOW_LENGTH:
		Extract_FCs_FixedWindowLength::get_all_FCs( station );
	case VARIABLE_WINDOW_LENGTH:
		Extract_FCs_VariableWindowLength::get_all_FCs( station );
	}
}

void Combination::define_combinations_for_rr( std::vector<StationBase> &station )
{
}

void Combination::define_combinations_for_ss( std::vector<StationBase> &station )
{
	Combination auxRRConc;
	auxRRConc.idxStn.resize(1);
	auxRRConc.idxTs.resize(1);

	// this is just the combination for single site
	for( int istn = 0; istn < station.size(); istn++ ) {
		for( int its = 0; its < station[istn].ts.size(); its++ ) {

			station[istn].ts[its].combination.resize(1);
			station[istn].ts[its].combination[0].idxStn.resize(1);
			station[istn].ts[its].combination[0].idxTs.resize(1);

			station[istn].ts[its].combination[0].numberOfConcomitantTs = 1;
			station[istn].ts[its].combination[0].idxStn[0] = istn;
			station[istn].ts[its].combination[0].idxTs[0] = its;

			station[istn].ts[its].combination[0].idxBgn = new Array<int>( 1, auxRRConc.numberOfConcomitantTs + 1 );
			station[istn].ts[its].combination[0].idxBgn[0](0,0) = station[istn].ts[its].timeVector[0](0,0);
			station[istn].ts[its].combination[0].idxBgn[0](0,1) = station[istn].ts[its].timeVector[0](0,0);

			station[istn].ts[its].combination[0].idxEnd = new Array<int>( 1, auxRRConc.numberOfConcomitantTs + 1 );
			station[istn].ts[its].combination[0].idxEnd[0](0,0) = station[istn].ts[its].timeVector[0](0,1);
			station[istn].ts[its].combination[0].idxEnd[0](0,1) = station[istn].ts[its].timeVector[0](0,1);
		}
	}
}

void Utils::delete_all( std::vector<StationBase> *station )
{
	for( int i = 0; i < station->size(); i++ ) {

	}
}

void Utils::draft_build_time_vector( std::vector<StationBase> &station )
{
	for( int istn = 0; istn < station.size(); istn++ ) {
		for( int its = 0; its < station[istn].ts.size(); its++ ) {
			station[istn].ts[its].timeVector = new Array<int>(1,2);
			station[istn].ts[its].timeVector[0](0,0) = 0;
			station[istn].ts[its].timeVector[0](0,1) = station[istn].ts[its].ch[0].timeSeries->getDim(0) - 1;
		}
	}
}