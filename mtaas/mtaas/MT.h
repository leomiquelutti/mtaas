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

class StationBase;
class ExtractorMTU;

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


// responsible for extracting stations' info, as name, path and type
class DirectoryProperties
{
public:
	
	DirectoryProperties(std::string aux):
		pathName(aux),
		mtuCounter(0), 
		aduCounter(0) 
	{ this->numberOfFolders = get_number_of_subfolders(); };
	
	std::string						pathName;
	std::vector<StationBase>	initialize_stations( std::string inputPath );
	size_t							get_number_of_tbl_files( const std::string inputSubPath );

private:

	size_t							mtuCounter;
	size_t							aduCounter;
	size_t							numberOfFolders;
		
	void							initialize_station_parameters( std::vector<StationBase> *station );
	void							fill_station_names( std::vector<StationBase> *station );
	void							fill_station_types( std::vector<StationBase> *station );

	void							define_tbl_files( StationBase *station );
	size_t							get_number_of_subfolders();
	void							get_tbl_names( std::string *TBLs, std::string file, size_t size );
	void							fill_in_dir_mtu_info( StationBase *station );
	size_t							count_TSn_files_for_each_tbl_file_and_get_its_name( StationBase *station, size_t idx, std::vector<std::string> **v );
};


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

	// raw data
	matCUDA::Array<double>	*timeSeries;
	matCUDA::Array<Complex>	*correctionToData;
	matCUDA::Array<double>	*timeVector;

	// parameters for processing
	size_t			samplingFrequency;
	double			maxFreqInHertz;
	double			minFreqInHertz;
	bool			isContinuous;
	unsigned int	continuousSamples;

	// functions
	void			initialize_parameters();
	bool			is_acquisition_continuos( Array<double> *timeVector );
	void			read_time_series();

private:
	double			exDipoleLength;
	double			eyDipoleLength;

	// functions
	//ExtractorMTU *read_time_series_mtu();
};


// main class contains all info regarding stations
class StationBase
{
public:
	
	//StationBase() {};
	//StationBase(const StationBase& element) {*this = element;};
 //   StationBase& operator = (const StationBase& element);     //which needs definition
	std::string		fileName;
	std::string		stationName;
	FILE_TYPE		fileType;
	Date			date;
	Position		position;
	size_t			amountOfTs;

	// one StationFile for each time-series (or analogous) file
	std::vector<StationFile> *ts;

	// extractors
	ExtractorMTU *mtu;

	// parameters for processing
	size_t			numberOfChannels;

	// functions
	//void			initialize_parameters();
	//bool			is_acquisition_continuos( Array<double> *timeVector );
	void			read_time_series( DirectoryProperties *dirInfo );

private:
	double			exDipoleLength;
	double			eyDipoleLength;

	// functions
	//ExtractorMTU *read_time_series_mtu();
};

#endif