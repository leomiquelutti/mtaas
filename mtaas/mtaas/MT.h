/***************************************************************************
 *   Copyright (C) 2013 by Leonardo Miquelutti, Bernardo Sampaio           *
 *   leomiquelutti@gmail.com, b_sampaio@yahoo.com	                       *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef MT_H
#define MT_H

class ExtractorMTU;
class MtuDirInfo;
class Adu07DirInfo;
class DirectoryProperties;

#include "stdafx.h"
#include "Mtu.h"

#include "boost/date_time/posix_time/posix_time.hpp"

#define PI 3.141592653589793238462643
#define pi 3.141592653589793238462643

#ifndef MAXPATH
#define MAXPATH 256
#endif

#define MININUM_OF_SEGMENTS_TO_COHERENCE_PRESELECTION 25

using namespace matCUDA;

// identify file type
typedef enum {
	FILE_TYPE_EMI = 1,
	FILE_TYPE_LMT,
	FILE_TYPE_SIO,
	FILE_TYPE_MTU,
	FILE_TYPE_ADU_07
} FILE_TYPE;

// sets whether the segment length will be variable or fixed (with decimation)
typedef enum {
	FIXED_WINDOW_LENGTH = 1,
	VARIABLE_WINDOW_LENGTH
} TS_TO_FFT;

// sets which kind of estimator will be used to retrieve the transfer tensors
typedef enum {
	LEAST_SQUARES = 1,
	WEIGHTED_LEAST_SQUARES,
	M_ESTIMATOR,
	BOUNDED_INFLUENCE,
	TWO_STEPS_BOUNDED_INFLUENCE
} ESTIMATOR_TYPE;

// responsible for time reference of station
class Date
{
public:
	int startYear;
	int startMonth;
	int startDay;
	int startHour;
	int startMinute;
	int startSecond;	
	boost::posix_time::ptime tsTime;

	std::string getDateStr();
};

// responsible for positioning of station
class Position
{
private:
	std::string latitude_str;
	std::string longitude_str;
	double latitude;
	double longitude;
	double elevation;
	double declination;
	double declinationInput;
};

// class containing info from each TS within station
class StationFile
{
public:

	StationFile() {};
	~StationFile() {};

	StationFile(const StationFile& element) {*this = element;};
    StationFile& operator = (const StationFile& element);

	// raw data
	matCUDA::Array<double>	*timeSeries;
	matCUDA::Array<Complex>	*correctionToData;
	matCUDA::Array<double>	*timeVector;

	// parameters for processing
	size_t			samplingFrequency;
	double			maxFreqInHertz;
	double			minFreqInHertz;
	bool			isContinuous;

	// extractors
	ExtractorMTU mtu;

	// functions
	bool			is_acquisition_continuos( Array<double> *timeVector );

private:
	double			exDipoleLength;
	double			eyDipoleLength;
};

// main class contains all info regarding stations
class StationBase
{
public:
	
	StationBase() :
		stationName("NULL"),
		path("NULL"),
		amountOfTs(0) {};

	StationBase(const StationBase& element) {*this = element;};
    StationBase& operator = (const StationBase& element);

	std::string		stationName;
	std::string		path;
	FILE_TYPE		fileType;
	Date			date;
	Position		position;
	size_t			amountOfTs;

	// one StationFile for each time-series (or analogous) file
	std::vector<StationFile> ts;

	// function to read time-series
	void			read_time_series( DirectoryProperties *dirInfo );
};

// utils
class Utils
{
public:

	static void delete_all( std::vector<StationBase> *station );
};

#endif