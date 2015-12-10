#include "MT.h"
#include "Mtu.h"
//#include <sstream>

StationBase& StationBase::operator = (const StationBase& element)
{
	this->amountOfTs = element.amountOfTs;
	this->date = element.date;
	this->fileType = element.fileType;
	this->path = element.path;
	this->position = element.position;
	this->stationName = element.stationName;
	this->ts = element.ts;

	return *this;
}

StationFile& StationFile::operator = (const StationFile& element)
{
	this->correctionToData = element.correctionToData;
	this->exDipoleLength = element.exDipoleLength;
	this->eyDipoleLength = element.eyDipoleLength;
	this->isContinuous = element.isContinuous;
	this->maxFreqInHertz = element.maxFreqInHertz;
	this->minFreqInHertz = element.minFreqInHertz;
	this->mtu = element.mtu;
	this->samplingFrequency = element.samplingFrequency;
	this->timeSeries = element.timeSeries;
	this->timeVector = element.timeVector;

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

//bool StationBase::is_acquisition_continuos( Array<double> *timeVector )
//{
//	// check consistency
//	for( int i = 1; i < timeVector->GetDescriptor().GetDim( 0 ); i++ ) {
//		if( (*timeVector)(i,0) - (*timeVector)(i-1,0) < 1 )
//		{
//			std::cout << "There's a problem in time vector" <<std::endl;
//			std::cout << "timeVector(" << i << ") - timeVector(" << i-1 << ") < 1" << std::endl;
//			exit(-42);
//		}
//	}
//
//	// check continuity
//	for( int i = 1; i < timeVector->GetDescriptor().GetDim( 0 ); i++ ) {
//		if( (*timeVector)(i,0) - (*timeVector)(i-1,0) > 1 )
//			return false;
//	}
//	return true;
//}

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

void Utils::delete_all( std::vector<StationBase> *station )
{
	for( int i = 0; i < station->size(); i++ ) {

	}
}