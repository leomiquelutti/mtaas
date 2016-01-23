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

class MtuDirInfo;
class Adu07DirInfo;
class DirectoryProperties;

#include "stdafx.h"
#include "Mtu.h"
#include "parameters.h"

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

class WriteOutputs
{
public:

	WriteOutputs( std::vector<StationBase> &station );

	void write_to_file( matCUDA::Array<double> freq, matCUDA::Array<double> rhoapp, matCUDA::Array<double> phase );
};

// class containing info from each Channel within each TS
class FrequencyResponses
{
public:
	
	FrequencyResponses() :
		elementCounter(0) {};
	~FrequencyResponses() {};

	FrequencyResponses(const FrequencyResponses& element) {*this = element;};
    FrequencyResponses& operator = (const FrequencyResponses& element);
	
	size_t							deciLevel;
	size_t							elementCounter;
	double							frequency;
	matCUDA::Array<ComplexDouble>	*in;
	matCUDA::Array<ComplexDouble>	*inRR;
	size_t							lowerLimitBandWidth;
	size_t							nFCsPerSegment;
	matCUDA::Array<ComplexDouble>	*out;
	Parameters						*param;
	matCUDA::Array<ComplexDouble>	*transferTensor;
	size_t							upperLimitBandWidth;
};

class Combination
{
public:

	Combination() :
		nConcomitantTs(1),
		nBlocks(1) {};
	~Combination() {};	

	Combination(const Combination& element) {*this = element;};
    Combination& operator = (const Combination& element);

	// index to reference for which Station and which Ts there is concomitance
	std::vector<size_t> 			deciLevelEachBlock;
	matCUDA::Array<int> 			*fcDistribution;
	std::vector<FrequencyResponses> fr;
	matCUDA::Array<int>				*idxBgn;
	matCUDA::Array<int> 			*idxEnd;
	std::vector<size_t> 			idxStn;
	std::vector<size_t> 			idxTs;
	std::vector<size_t> 			lengthTsEachBlock;
	std::vector<double>				measuredFrequencies;
	size_t 							nBlocks;
	size_t 							nConcomitantTs;
	size_t 							nDeci;
	size_t 							nFreq;
	std::vector<double> 			samplFreqEachDeciLevel;

	static void define_combinations_for_rr( std::vector<StationBase> &station );
	static void define_combinations_for_ss( std::vector<StationBase> &station );
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
	
	matCUDA::Array<double>			*arCoeff;
	CHANNEL_ORIENTATION				channel_orientation;
	CHANNEL_TYPE					channel_type;
	matCUDA::Array<ComplexDouble>	*correction;
	double							countConversion;
	double							dipoleLength;
	matCUDA::Array<ComplexDouble>	*fc;
	double							gain;
	std::string						name;
	double							orientationHor;
	double							orientationVer;
	matCUDA::Array<ComplexDouble>	*systemResponse;
	matCUDA::Array<double>			*systemResponseFreqs;
	matCUDA::Array<double>			*timeSeries;
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
	
	std::vector<Combination>	combination;
	std::vector<Channel>		ch;
	Date						date;
	double						exDipoleLength;
	double						eyDipoleLength;
	bool						isAcquisitionContinuous;
	bool						isContinuous;
	bool						isThereRR;
	size_t						maxDecimationLevel;
	matCUDA::Array<int>			*maxFcDist;
	double						maxFreqInHertz;
	double						minFreqInHertz;
	ExtractorMTU				mtu;
	size_t						nChannels;
	size_t						numberOfCombinations;
	double						samplingFrequency;
	matCUDA::Array<int>			*timeVector;

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

	// function to extract corrected FC coefficients to evaluate transfer functions
	static void		get_Z( std::vector<StationBase> &station, ESTIMATOR_TYPE estimator_type );
};

// utils
class Utils
{
public:

	static void delete_all( std::vector<StationBase> *station );
	static void draft_build_time_vector( std::vector<StationBase> &station );
};

#endif