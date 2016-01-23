#include "MT.h"
#include "Mtu.h"
#include "ts2fc.h"
#include "fc2z.h"

class Extract_FCs_FixedWindowLength;
class Extract_FCs_VariableWindowLength;
class Evaluate_Z_LS;
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
	this->combination = element.combination;
	this->ch = element.ch;
	this->date = element.date;
	this->exDipoleLength = element.exDipoleLength;
	this->eyDipoleLength = element.eyDipoleLength;
	this->isAcquisitionContinuous = element.isAcquisitionContinuous;
	this->isContinuous = element.isContinuous;
	this->isThereRR = element.isThereRR;
	this->maxDecimationLevel = element.maxDecimationLevel;
	this->maxFcDist = element.maxFcDist;
	this->maxFreqInHertz = element.maxFreqInHertz;
	this->minFreqInHertz = element.minFreqInHertz;
	this->mtu = element.mtu;
	this->nChannels = element.nChannels;
	this->numberOfCombinations = element.numberOfCombinations;
	this->samplingFrequency = element.samplingFrequency;
	this->timeVector = element.timeVector;

	return *this;
}

Channel& Channel::operator = (const Channel& element)
{	
	this->arCoeff = element.arCoeff;
	this->channel_orientation = element.channel_orientation;
	this->channel_type = element.channel_type;
	this->correction = element.correction;
	this->countConversion = element.countConversion;
	this->dipoleLength = element.dipoleLength;
	this->fc = element.fc;
	this->gain = element.gain;
	this->name = element.name;
	this->orientationHor = element.orientationHor;
	this->orientationVer = element.orientationVer;
	this->systemResponse = element.systemResponse;
	this->systemResponseFreqs = element.systemResponseFreqs;
	this->timeSeries = element.timeSeries;

	return *this;
}

Combination& Combination::operator = (const Combination& element)
{	
	this->deciLevelEachBlock = element.deciLevelEachBlock;
	this->fcDistribution = element.fcDistribution;
	this->fr = element.fr;
	this->idxBgn = element.idxBgn;
	this->idxEnd = element.idxEnd;
	this->idxStn = element.idxStn;
	this->idxTs = element.idxTs;
	this->lengthTsEachBlock = element.lengthTsEachBlock;
	this->measuredFrequencies = element.measuredFrequencies;
	this->nBlocks = element.nBlocks;
	this->nConcomitantTs = element.nConcomitantTs;
	this->nDeci = element.nDeci;
	this->nFreq = element.nFreq;
	this->samplFreqEachDeciLevel = element.samplFreqEachDeciLevel;

	return *this;
}

FrequencyResponses& FrequencyResponses::operator = (const FrequencyResponses& element)
{
	this->deciLevel = element.deciLevel;
	this->elementCounter = element.elementCounter;
	this->frequency = element.frequency;
	this->in = element.in;
	this->inRR = element.inRR;
	this->lowerLimitBandWidth = element.lowerLimitBandWidth;
	this->nFCsPerSegment = element.nFCsPerSegment;
	this->out = element.out;
	this->param = element.param;
	this->transferTensor = element.transferTensor;
	this->upperLimitBandWidth = element.upperLimitBandWidth;

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

	//// checkouts
	//for( int istn = 0; istn < station.size(); istn++ ) {
	//	for( int its = 0; its < station[istn].ts.size(); its++ ) {
	//		for( int icomb = 0; icomb < station[istn].ts[its].combination.size(); icomb++ ) {
	//			cout << station[istn].ts[its].combination[icomb].idxEnd[0](0,0) - station[istn].ts[its].combination[icomb].idxBgn[0](0,0) + 1 << endl;
	//		}
	//		cout << endl;
	//	}
	//			std::cout << std::endl;
	//}

	switch(ts2fft_type) 
	{
	case FIXED_WINDOW_LENGTH:
		Extract_FCs_FixedWindowLength::get_all_FCs( station );
	case VARIABLE_WINDOW_LENGTH:
		Extract_FCs_VariableWindowLength::get_all_FCs( station );
	}
}

void StationBase::get_Z( std::vector<StationBase> &station, ESTIMATOR_TYPE estimator_type )
{
	switch(estimator_type)
	{
	case LEAST_SQUARES:
		Evaluate_Z_LS::get_Z( station );
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

			station[istn].ts[its].combination[0].nConcomitantTs = 1;
			station[istn].ts[its].combination[0].nBlocks = 1;
			station[istn].ts[its].combination[0].idxStn[0] = istn;
			station[istn].ts[its].combination[0].idxTs[0] = its;

			station[istn].ts[its].combination[0].idxBgn = new Array<int>( 1, auxRRConc.nConcomitantTs + 1 );
			station[istn].ts[its].combination[0].idxBgn[0](0,0) = station[istn].ts[its].timeVector[0](0,0);
			station[istn].ts[its].combination[0].idxBgn[0](0,1) = station[istn].ts[its].timeVector[0](0,0);

			station[istn].ts[its].combination[0].idxEnd = new Array<int>( 1, auxRRConc.nConcomitantTs + 1 );
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

WriteOutputs::WriteOutputs( std::vector<StationBase> &station )
{
	for( int istn = 0; istn < station.size(); istn++ )
		for( int its = 0; its < station[istn].ts.size(); its++ )
			for( int icomb = 0; icomb < station[istn].ts[its].combination.size(); icomb++ ) {

				size_t size = station[istn].ts[its].combination[icomb].fr.size();
				matCUDA::Array<double> freq( size );
				matCUDA::Array<double> rhoapp( size, 4 );
				matCUDA::Array<double> phase( size, 4 );

				for( int ifreq = 0; ifreq < station[istn].ts[its].combination[icomb].fr.size(); ifreq++ ) {
					freq( ifreq ) = station[istn].ts[its].combination[icomb].fr[ifreq].frequency;

					rhoapp( ifreq, 0 ) = station[istn].ts[its].combination[icomb].fr[ifreq].param->apparentResistivity( 0, 0 );
					rhoapp( ifreq, 1 ) = station[istn].ts[its].combination[icomb].fr[ifreq].param->apparentResistivity( 1, 0 );
					rhoapp( ifreq, 2 ) = station[istn].ts[its].combination[icomb].fr[ifreq].param->apparentResistivity( 0, 1 );
					rhoapp( ifreq, 3 ) = station[istn].ts[its].combination[icomb].fr[ifreq].param->apparentResistivity( 1, 1 );

					phase( ifreq, 0 ) = station[istn].ts[its].combination[icomb].fr[ifreq].param->phase( 0, 0 );
					phase( ifreq, 1 ) = station[istn].ts[its].combination[icomb].fr[ifreq].param->phase( 1, 0 );
					phase( ifreq, 2 ) = station[istn].ts[its].combination[icomb].fr[ifreq].param->phase( 0, 1 );
					phase( ifreq, 3 ) = station[istn].ts[its].combination[icomb].fr[ifreq].param->phase( 1, 1 );
				}

				this->write_to_file( freq, rhoapp, phase );
			}
}

void WriteOutputs::write_to_file( matCUDA::Array<double> freq, matCUDA::Array<double> rhoapp, matCUDA::Array<double> phase )
{
}