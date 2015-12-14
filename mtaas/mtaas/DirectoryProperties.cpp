#include "DirectoryProperties.h"
#include "MT.h"
#include "Mtu.h"

#include <iostream>
#include <boost/filesystem.hpp>
#include <boost/lambda/bind.hpp>

using namespace boost::filesystem;
using namespace boost::filesystem::path_traits;

std::vector<StationBase> DirectoryProperties::initialize_stations()
{
	// create vector of 'numberOfFolders' StationBases
	std::vector<StationBase> station( this->numberOfFolders );
	
	// initialize stations parameters, as name, path, type and its types derivatives
	this->initialize_station_parameters( &station );

	return station;
}

void DirectoryProperties::initialize_station_parameters( std::vector<StationBase> *station )
{
	// get stations names
	this->fill_station_names( station );

	// get station type and stores paths for TS and complementary needed files, 
	// accordingly to FILE_TYPE (equipment type)
	this->fill_station_types_and_paths( station );
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
				for (vec::const_iterator it (v.begin()); it != v.end(); ++it) {
					(*station)[counter].stationName = it->filename().string();
					(*station)[counter++].path = this->pathName;
				}
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

void DirectoryProperties::fill_station_types_and_paths( std::vector<StationBase> *station )
{
	size_t nTbl = 0;
		
	for( int i = 0; i < this->numberOfFolders; i++ ) {

		nTbl = this->get_number_of_tbl_files( &(*station)[i] );

		// MTU type
		if( nTbl > 0 ) {
			(*station)[i].fileType = FILE_TYPE_MTU;
			MtuDirInfo::fill_in_dir_mtu_info( &(*station)[i], nTbl );
		}
	}
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

size_t DirectoryProperties::get_number_of_tbl_files( StationBase *station )
{
	path p( this->pathName + "\\" + station->stationName );

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

size_t DirectoryProperties::get_number_of_tsn_files( StationBase *station )
{	
	size_t fileCounter = 0;
	size_t tsCounter = 0;
	StationFile auxTs;

	path p( this->pathName + "\\" + station->stationName );

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

				// get amount of TSn files
				for (vec::const_iterator it (v.begin()); it != v.end(); ++it)
					if( it->filename().string().find( ".TS" ) != std::string::npos )
						fileCounter++;

				//// get TSn and its respective TBL and CTS files
				//for( int i = 0; i < nTbl; i++ ) {
				//	for (vec::const_iterator it (v.begin()); it != v.end(); ++it) {
				//		if( it->filename().string().find( auxTblNames[i].substr( 0, auxTblNames[i].find( ".TBL" ) ) + ".TS" ) != std::string::npos ) {
				//			auxTs.mtu->tblFile = auxTblNames[i];
				//			auxTs.mtu->tsnFile = it->filename().string();
				//			auxTs.mtu->mtuTsBand = phoenixTsBand_t( atoi( &(auxTs.mtu->tsnFile.back()) ) );
				//			auxTs.mtu->ctsFile = it->filename().string().substr( 0, auxTblNames[i].find( ".TBL" ) ) + ".CTS";
				//			//station->ts[tsCounter++].mtu = new ExtractorMTU;
				//			station->ts.push_back(auxTs);
				//			tsCounter++;
				//		}
				//	}
				//}
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

	return fileCounter;
}

void MtuDirInfo::fill_in_dir_mtu_info( StationBase *station, size_t nTbl )
{	
	size_t fileCounter = 0;
	vector<string> auxTblNames(nTbl);
	StationFile auxTs;

	string auxPath = station->path + "\\" + station->stationName;
	path p( auxPath );

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

				// get tbl names
				for (vec::const_iterator it (v.begin()); it != v.end(); ++it)
				{
					if( it->filename().string().find( ".TBL" ) != std::string::npos )
						auxTblNames[ fileCounter++ ] = it->filename().string();
				}

				// get TSn and its respective TBL and CTS files
				size_t tsCounter = 0;
				for( int i = 0; i < nTbl; i++ ) {
					for (vec::const_iterator it (v.begin()); it != v.end(); ++it) {
						if( it->filename().string().find( auxTblNames[i].substr( 0, auxTblNames[i].find( ".TBL" ) ) + ".TS" ) != std::string::npos ) {
							auxTs.mtu.tsnFile = it->string();
							auxTs.mtu.tblFile = auxPath + "\\" + auxTblNames[i];
							auxTs.mtu.mtuTsBand = phoenixTsBand_t( atoi( &(auxTs.mtu.tsnFile.back()) ) );
							auxTs.mtu.ctsFile = auxPath + "\\" + it->filename().string().substr( 0, auxTblNames[i].find( ".TBL" ) ) + ".CTS";

							station->ts.push_back(auxTs);
							tsCounter++;
						}
					}
				}
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

//size_t DirectoryProperties::count_TSn_files_for_each_tbl_file_and_get_its_name( StationBase *station, size_t idx, std::vector<std::string> **vectorOfTSnFiles )
//{
//	//size_t fileCounter = 0;
//	//std::string inputSubPath = this->pathName + "\\" + station->stationName;
//
//	//path p( inputSubPath );
//	//try
//	//{
//	//	if( exists( p ) )
//	//	{
//	//		if( is_directory( p ) )
//	//		{
//	//			typedef std::vector<path> vec;             // store paths,
//	//			vec v;                                // so we can sort them later
//
//	//			copy(directory_iterator(p), directory_iterator(), back_inserter(v));
//
//	//			sort(v.begin(), v.end());             // sort, since directory iteration is not ordered on some file systems
//
//	//			for (vec::const_iterator it (v.begin()); it != v.end(); ++it)
//	//			{
//	//				if( it->filename().string().find( station->mtu->inputTbl[idx] + ".TS" ) != std::string::npos ) {
//	//					//vectorOfTSnFiles[idx][fileCounter++] = "";
//	//				}
//	//			}
//	//		}
//	//		else
//	//		{
//	//			std::cout << p << " exists, but no TSn file could be associated with " << station->mtu->inputTbl[idx] << ".TBL file" << std::endl;
//	//			return 0;
//	//		}
//	//	}
//	//	else
//	//	{
//	//		std::cout << p << " does not exist" << std::endl;
//	//		return 0;
//	//	}
//	//}
//	//catch( const filesystem_error& ex )
//	//{ 
//	//	std::cout << ex.what() << '\n' << std::endl; 
//	//	return 0;
//	//}
//
//	//return fileCounter;
//
//	return 0;
//}
//
//void DirectoryProperties::define_tbl_files( StationBase *station )
//{
//	//std::string inputSubPath = this->pathName + "\\" + station->stationName;
//	//path p( inputSubPath );
//
//	//size_t fileCounter = 0;
//	//try
//	//{
//	//	if( exists( p ) )
//	//	{
//	//		if( is_directory( p ) )
//	//		{
//	//			typedef std::vector<path> vec;             // store paths,
//	//			vec v;                                // so we can sort them later
//
//	//			copy(directory_iterator(p), directory_iterator(), back_inserter(v));
//
//	//			sort(v.begin(), v.end());             // sort, since directory iteration is not ordered on some file systems
//
//	//			for (vec::const_iterator it (v.begin()); it != v.end(); ++it)
//	//			{
//	//				if( it->filename().string().find( ".TBL" ) != std::string::npos )
//	//					station->mtu->inputTbl[ fileCounter++ ] = it->filename().string().substr( 0, it->filename().string().find( ".TBL" ) );
//	//			}
//	//		}
//	//		else
//	//		{
//	//			std::cout << p << " exists, but is neither a regular file nor a directory" << std::endl;
//	//			return;
//	//		}
//	//	}
//	//	else
//	//	{
//	//		std::cout << p << " does not exist" << std::endl;
//	//		return;
//	//	}
//	//}
//	//catch( const filesystem_error& ex )
//	//{ 
//	//	std::cout << ex.what() << '\n' << std::endl; 
//	//	return;
//	//}
//}
//
//void DirectoryProperties::get_tbl_names( std::string *TBLs, std::string file, size_t size )
//{
//	std::ifstream infile;
//	infile.open( file.c_str(), std::ios::in );
//
//	for( int i = 0; i < size; ++i )
//		std::getline( infile, TBLs[ i ] );
//
//	infile.close();
//}