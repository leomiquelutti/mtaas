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

// identify file type
typedef enum {
	CHANNEL_TYPE_H = 1,
	CHANNEL_TYPE_E
} CHANNEL_TYPE;

// identify file orientation (or label)
typedef enum {
	CHANNEL_ORIENTATION_X = 1,
	CHANNEL_ORIENTATION_Y,
	CHANNEL_ORIENTATION_Z,
} CHANNEL_ORIENTATION;

// sets whether the segment length will be variable or fixed (with decimation)
typedef enum {
	FIXED_WINDOW_LENGTH = 1,
	VARIABLE_WINDOW_LENGTH
} TS_TO_FFT_TYPE;

// sets which kind of estimator will be used to retrieve the transfer tensors
typedef enum {
	LEAST_SQUARES = 1,
	WEIGHTED_LEAST_SQUARES,
	M_ESTIMATOR,
	BOUNDED_INFLUENCE,
	TWO_STEPS_BOUNDED_INFLUENCE
} ESTIMATOR_TYPE;

// sets which kind of estimator will be used to retrieve the transfer tensors
typedef enum {
	SINGLE_SITE = 1,
	REMOTE_REFERENCE
} RR_OR_SS_TYPE;

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

class Combination
{
public:

	Combination() :
		numberOfConcomitantTs(1) {};
	~Combination() {};	

	Combination(const Combination& element) {*this = element;};
    Combination& operator = (const Combination& element);

	// index to reference for which Station and which Ts there is concomitance
	size_t numberOfConcomitantTs;
	std::vector<size_t> idxStn;
	std::vector<size_t> idxTs;
	matCUDA::Array<int> *idxBgn;
	matCUDA::Array<int> *idxEnd;

	static void define_combinations_for_rr( std::vector<StationBase> &station );
	static void define_combinations_for_ss( std::vector<StationBase> &station );
};

// class containing info from each Channel within each TS
class FrequencyResponses
{
public:
	
	FrequencyResponses() :
		numberOfStationsInvolved(1) {};
	~FrequencyResponses() {};

	FrequencyResponses(const FrequencyResponses& element) {*this = element;};
    FrequencyResponses& operator = (const FrequencyResponses& element);
	
	double frequency;
	size_t numberOfStationsInvolved;
	matCUDA::Array<ComplexDouble> *in;
	matCUDA::Array<ComplexDouble> *inRR;
	matCUDA::Array<ComplexDouble> *out;
};

// class containing info from each Channel within each TS
class Channel
{
public:

	
	Channel() : 
		countConversion(1),
		gain(1),
		dipoleLength(1) {};

	~Channel() {
		//delete arCoeff;
		//delete systemResponse;
		//delete timeSeries;
	}

	Channel(const Channel& element) {*this = element;};
    Channel& operator = (const Channel& element);
	
	matCUDA::Array<double> *timeSeries;
	matCUDA::Array<double> *arCoeff;
	matCUDA::Array<ComplexDouble> *systemResponse;
	matCUDA::Array<double>	*systemResponseFreqs;
	matCUDA::Array<ComplexDouble> *correction;
	//matCUDA::Array<ComplexDouble> *fc;
	double countConversion;
	double orientationHor;
	double orientationVer;
	double dipoleLength;
	std::string name;
	CHANNEL_TYPE channel_type;
	CHANNEL_ORIENTATION channel_orientation;
	double gain;
};

// class containing info from each TS within station
class StationFile
{
public:

	StationFile():
		isThereRR(false),
		numberOfCombinations(1),
		isAcquisitionContinuous(true) {};

	~StationFile() {
		/*delete timeVector; */ };

	StationFile(const StationFile& element) {*this = element;};
    StationFile& operator = (const StationFile& element);

	matCUDA::Array<int>	*timeVector;

	// parameters for processing
	double			samplingFrequency;
	double			maxFreqInHertz;
	double			minFreqInHertz;
	bool			isContinuous;
	size_t			nChannels;
	double			exDipoleLength;
	double			eyDipoleLength;
	Date			date;

	bool			isThereRR;
	size_t			numberOfCombinations;
	bool			isAcquisitionContinuous;

	// extractors
	ExtractorMTU mtu;

	// to store time-series for each channel and its respective fourier coefficients
	std::vector<Channel> ch;

	// to store FCs corrected for each frequency, in order to calculate Z
	// stores the possible different combinations allowed
	std::vector<Combination> combination;
	std::vector<std::vector<FrequencyResponses>> fr;

	// functions
	bool			is_acquisition_continuos( Array<double> *timeVector );

private:
};

// main class contains all info regarding stations
class StationBase
{
public:
	
	StationBase() :
		stationName("NULL"),
		path("NULL") {};

	StationBase(const StationBase& element) {*this = element;};
    StationBase& operator = (const StationBase& element);

	std::string		stationName;
	std::string		path;
	FILE_TYPE		fileType;
	Position		position;

	// one StationFile for each time-series (or analogous) file
	std::vector<StationFile> ts;

	// function to read time-series
	void			read_time_series( DirectoryProperties *dirInfo );

	// function to extract corrected FC coefficients to evaluate transfer functions
	static void		get_all_FCs( std::vector<StationBase> &station, TS_TO_FFT_TYPE ts2fft_type, RR_OR_SS_TYPE rrorss_type );
};

// utils
class Utils
{
public:

	static void delete_all( std::vector<StationBase> *station );
	static void draft_build_time_vector( std::vector<StationBase> &station );
};

#endif