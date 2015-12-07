#include "MT.h"

#include <iostream>
#include <boost/filesystem.hpp>
#include <boost/lambda/bind.hpp>

using namespace boost::filesystem;
using namespace boost::filesystem::path_traits;

std::vector<StationBase> DirectoryProperties::initialize_stations( std::string inputPath )
{
	// create vector of 'numberOfFolders' StationBases
	std::vector<StationBase> station( this->numberOfFolders );
	
	// initialize stations parameters, as name, path, type and its types derivatives
	this->initialize_station_parameters( &station );

	return station;
}

void DirectoryProperties::initialize_station_parameters( std::vector<StationBase> *station )
{
	this->fill_station_names( station );

	this->fill_station_types( station );
}

void DirectoryProperties::fill_station_names( std::vector<StationBase> *station )
{
	path p( this->pathName );

	try
	{
		if( exists( p ) )
		{
			if( is_directory( p ) )
			{
				typedef std::vector<path> vec;             // store paths,
				vec v;                                // so we can sort them later

				copy(directory_iterator(p), directory_iterator(), back_inserter(v));

				sort(v.begin(), v.end());             // sort, since directory iteration is not ordered on some file systems

				size_t counter = 0;
				for (vec::const_iterator it (v.begin()); it != v.end(); ++it)
					(*station)[counter++].stationName = it->filename().string();
					//station[counter++].stationName = it->filename().string();
			}
			else
			{
				std::cout << p << " exists, but is neither a regular file nor a directory" << std::endl;
			}
		}
		else
		{
			std::cout << p << " does not exist" << std::endl;
		}
	}
	catch( const filesystem_error& ex )
	{ 
		std::cout << ex.what() << '\n'; 
	}
}

void DirectoryProperties::fill_station_types( std::vector<StationBase> *station )
{
	size_t amountOfTbl = 0;
	size_t amountOfTsnPerTbl = 0;
	
	for( int i = 0; i < this->numberOfFolders; i++ ) {

		amountOfTbl = this->get_number_of_tbl_files( (*station)[i].stationName );

		if( amountOfTbl > 0 ) {
			(*station)[i].fileType = FILE_TYPE_MTU;
			fill_in_dir_mtu_info( &(*station)[i] );
		}
	}
}

void DirectoryProperties::fill_in_dir_mtu_info( StationBase *station )
{
	std::vector<string> **aux;
	std::string **auxTSnFiles;
	size_t auxTSnCounter = 0;

	station->mtu = new ExtractorMTU;
	station->mtu->amountOfTbl = this->get_number_of_tbl_files( station->stationName ); 
	*aux = new std::vector<string>[station->mtu->amountOfTbl];
	*auxTSnFiles = new std::string[station->mtu->amountOfTbl];

	// get amount of TSn files for all TBL files
	station->ts = new std::vector<StationFile>[1];
	station->mtu->inputTbl = new std::string[station->mtu->amountOfTbl];
	for( int i = 0; i < station->mtu->amountOfTbl; i++ ) {
		define_tbl_files( station );
		auxTSnCounter = count_TSn_files_for_each_tbl_file_and_get_its_name( station, i, aux );
		station->mtu->amountOfTSn += auxTSnCounter;
	}

	// fill in station->ts with each TSn file
	for( int i = 0; i < station->mtu->amountOfTbl; i++ ) {
		define_tbl_files( station );
		station->mtu->amountOfTSn += count_TSn_files_for_each_tbl_file_and_get_its_name( station, i, aux );
	}
}

size_t DirectoryProperties::count_TSn_files_for_each_tbl_file_and_get_its_name( StationBase *station, size_t idx, std::vector<std::string> **vectorOfTSnFiles )
{
	size_t fileCounter = 0;
	std::string inputSubPath = this->pathName + "\\" + station->stationName;

	path p( inputSubPath );
	try
	{
		if( exists( p ) )
		{
			if( is_directory( p ) )
			{
				typedef std::vector<path> vec;             // store paths,
				vec v;                                // so we can sort them later

				copy(directory_iterator(p), directory_iterator(), back_inserter(v));

				sort(v.begin(), v.end());             // sort, since directory iteration is not ordered on some file systems

				for (vec::const_iterator it (v.begin()); it != v.end(); ++it)
				{
					if( it->filename().string().find( station->mtu->inputTbl[idx] + ".TS" ) != std::string::npos ) {
						//vectorOfTSnFiles[idx][fileCounter++] = "";
					}
				}
			}
			else
			{
				std::cout << p << " exists, but no TSn file could be associated with " << station->mtu->inputTbl[idx] << ".TBL file" << std::endl;
				return 0;
			}
		}
		else
		{
			std::cout << p << " does not exist" << std::endl;
			return 0;
		}
	}
	catch( const filesystem_error& ex )
	{ 
		std::cout << ex.what() << '\n' << std::endl; 
		return 0;
	}

	return fileCounter;
}

size_t DirectoryProperties::get_number_of_subfolders()
{
	using namespace boost::lambda;

	path p( this->pathName );   

	size_t subfoldersCounter = std::count_if(
        directory_iterator(p),
        directory_iterator(),
        boost::lambda::bind( static_cast<bool(*)(const path&)>(is_directory),
        boost::lambda::bind( &directory_entry::path, _1 ) ) );

	return subfoldersCounter;
}

void DirectoryProperties::define_tbl_files( StationBase *station )
{
	std::string inputSubPath = this->pathName + "\\" + station->stationName;
	path p( inputSubPath );

	size_t fileCounter = 0;
	try
	{
		if( exists( p ) )
		{
			if( is_directory( p ) )
			{
				typedef std::vector<path> vec;             // store paths,
				vec v;                                // so we can sort them later

				copy(directory_iterator(p), directory_iterator(), back_inserter(v));

				sort(v.begin(), v.end());             // sort, since directory iteration is not ordered on some file systems

				for (vec::const_iterator it (v.begin()); it != v.end(); ++it)
				{
					if( it->filename().string().find( ".TBL" ) != std::string::npos )
						station->mtu->inputTbl[ fileCounter++ ] = it->filename().string().substr( 0, it->filename().string().find( ".TBL" ) );
				}
			}
			else
			{
				std::cout << p << " exists, but is neither a regular file nor a directory" << std::endl;
				return;
			}
		}
		else
		{
			std::cout << p << " does not exist" << std::endl;
			return;
		}
	}
	catch( const filesystem_error& ex )
	{ 
		std::cout << ex.what() << '\n' << std::endl; 
		return;
	}
}

size_t DirectoryProperties::get_number_of_tbl_files( const std::string inputSubPath )
{
	path p( this->pathName + "\\" + inputSubPath );

	size_t fileCounter = 0;
	try
	{
		if( exists( p ) )
		{
			if( is_directory( p ) )
			{
				typedef std::vector<path> vec;             // store paths,
				vec v;                                // so we can sort them later

				copy(directory_iterator(p), directory_iterator(), back_inserter(v));

				sort(v.begin(), v.end());             // sort, since directory iteration is not ordered on some file systems

				for (vec::const_iterator it (v.begin()); it != v.end(); ++it)
				{
					//std::cout << it->filename().string() << std::endl;
					if( it->filename().string().find( ".TBL" ) != std::string::npos )
						fileCounter++;
				}
			}
			else
			{
				std::cout << p << " exists, but is neither a regular file nor a directory" << std::endl;
				return fileCounter;
			}
		}
		else
		{
			std::cout << p << " does not exist" << std::endl;
			return fileCounter;
		}
	}
	catch( const filesystem_error& ex )
	{ 
		std::cout << ex.what() << '\n'; 
		return fileCounter;
	}

	return fileCounter;
}

void DirectoryProperties::get_tbl_names( std::string *TBLs, std::string file, size_t size )
{
	std::ifstream infile;
	infile.open( file.c_str(), std::ios::in );

	for( int i = 0; i < size; ++i )
		std::getline( infile, TBLs[ i ] );

	infile.close();
}