#include "MT.h"
#include "Mtu.h"
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

	return *this;
}

Channel& Channel::operator = (const Channel& element)
{
	this->arCoeff = element.arCoeff;
	this->timeSeries = element.timeSeries;
	this->correction = element.correction;
	this->countConversion = element.countConversion;
	this->orientationHor = element.orientationHor;
	this->orientationVer = element.orientationVer;
	this->name = element.name;
	this->gain = element.gain;

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

void StationBase::get_FCs()
{
}

void Utils::delete_all( std::vector<StationBase> *station )
{
	for( int i = 0; i < station->size(); i++ ) {

	}
}