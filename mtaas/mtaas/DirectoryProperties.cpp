#include "MT.h"

#include <iostream>
#include <boost/filesystem.hpp>
#include <boost/lambda/bind.hpp>

using namespace boost::filesystem;
using namespace boost::filesystem::path_traits;

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

void DirectoryProperties::fill_station_names( StationBase *station )
{
	path p( this->pathName );

	int a = this->numberOfFolders;
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
					station[counter++].stationName = it->filename().string();
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

std::vector<StationBase> DirectoryProperties::initialize_stations( std::string inputPath )
{
	DirectoryProperties dirInfo;
	const size_t numberOfFolders = dirInfo.get_number_of_subfolders( inputPath );

	std::vector<StationBase> station( numberOfFolders );
	
	dirInfo.fill_station_names( &station );

	return station;
}

//StationBase* DirectoryProperties::initialize_stations()
//{
//	this->numberOfFolders = this->get_number_of_subfolders( this->pathName );
//
//	StationBase *station = new StationBase[ this->numberOfFolders ];
//
//	this->fill_station_names( station );
//
//	return station;
//}

size_t DirectoryProperties::get_number_of_subfolders( std::string inputPath )
{
	using namespace boost::lambda;

	path p( inputPath );   

	size_t subfoldersCounter = std::count_if(
        directory_iterator(p),
        directory_iterator(),
        boost::lambda::bind( static_cast<bool(*)(const path&)>(is_directory),
        boost::lambda::bind( &directory_entry::path, _1 ) ) );

	return subfoldersCounter;
}

void DirectoryProperties::define_files( const std::string inputSubPath, size_t numberOfFiles, std::string *files )
{
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
						files[ fileCounter++ ] = it->filename().string().substr( 0, it->filename().string().find( ".TBL" ) );
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
					std::cout << it->filename().string() << std::endl;
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