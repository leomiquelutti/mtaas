#ifndef DIRECTORY_PROPERTIES_H
#define DIRECTORY_PROPERTIES_H

#include <vector>
#include <string>

class StationBase;

//// responsible for extracting information from folders
//class DirInfo
//{
//public:
//	size_t nTsFiles;
//
//	virtual void foo() = 0;
//	virtual size_t nTbl() = 0;
//	virtual size_t nXml() = 0;
//
//};
//
//// responsible for extracting information from ADU07 folders
//class Adu07DirInfo : public DirInfo
//{
//public:
//	void foo() {};
//	size_t nXml() { return 1; }
//
//};
//
//// responsible for extracting information from MTU folders
//class MtuDirInfo : public DirInfo
//{
//public:
//	void foo() {};
//	size_t nTbl() { return 1; }
//
//};

// responsible for extracting information from ADU07 folders
class Adu07DirInfo
{
public:
	size_t nXml;
};

// responsible for extracting information from MTU folders
class MtuDirInfo
{
public:
	size_t nTbl;
	static void fill_in_dir_mtu_info( StationBase *station );
};

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
	std::vector<StationBase>	initialize_stations();
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
	size_t							count_TSn_files_for_each_tbl_file_and_get_its_name( StationBase *station, size_t idx, std::vector<std::string> **v );
};

#endif