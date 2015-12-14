#include "MT.h"
#include "Mtu.h"
#include "GetDeclination.h"

#include <cctype>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
//#include <vector>
#include <sstream>
#include <boost\iostreams\device\array.hpp>
#include <boost/exception/all.hpp>

using namespace std;
using namespace matCUDA;

ExtractorMTU& ExtractorMTU::operator = (const ExtractorMTU& element)
{
	this->ctsFile = element.ctsFile;
	this->tblFile = element.tblFile;
	this->tsnFile = element.tsnFile;
	this->FS = element.FS;
	this->mtuTsBand = element.mtuTsBand;
	//this->nChannels = element.nChannels;
	this->nScansPerRecord = element.nScansPerRecord;
	this->numberOfBytes = element.numberOfBytes;
	this->numberOfRecords = element.numberOfRecords;
	this->numberOfSamples = element.numberOfSamples;
	this->recordLength = element.recordLength;
	this->sampleLength = element.sampleLength;
	this->tagLength = element.tagLength;
	this->tagsize = element.tagsize;
	this->tbl = element.tbl;

	return *this;
}

void ExtractorMTU::read_time_series( StationBase *station, DirectoryProperties *dirInfo ) 
{
	StationFile auxTs;
	for( int i = 0; i < station->ts.size(); i++ ) {
		auxTs = station->ts[i];

		// get tbl parameters
		auxTs = station->ts[i].mtu.get_parameters( auxTs );

		station->ts[i].mtu.read_mtu_data( auxTs );

		station->ts[i] = auxTs;
	}
}

StationFile ExtractorMTU::get_parameters( StationFile auxTs )
{
	// read tbl an stores result in table tbl file
	auxTs.mtu.tbl = auxTs.mtu.read_tbl();

	// put important parameters on auxTs to be returned
	auxTs.exDipoleLength = auxTs.mtu.tbl.exln;
	auxTs.eyDipoleLength = auxTs.mtu.tbl.eyln;

	// read TSn first tag - filling of auxTs inside function
	ifstream infile( this->tsnFile, ios::binary );
	if ( infile.good() ) {
		read_TSn_tag( infile, &auxTs, 0 );
		infile.close();
    }
	
	return auxTs;
}

void ExtractorMTU::read_mtu_data( StationFile &auxTs  )
{
	std::string input = auxTs.mtu.tsnFile;

	ifstream infile( input, ios::binary );
	if ( infile.good() )
    {
		read_TSn_tag( infile, &auxTs, 0 );
		infile.close();
    }
	
	auxTs.mtu.numberOfBytes = get_number_of_bytes( input );
	auxTs.mtu.numberOfRecords = auxTs.mtu.numberOfBytes / auxTs.mtu.recordLength;
	auxTs.mtu.numberOfSamples = auxTs.mtu.numberOfBytes / auxTs.mtu.recordLength * auxTs.mtu.nScansPerRecord;

	//Array<double> data( (size_t)auxTs.mtu.numberOfSamples, (const index_t)auxTs.nChannels );
	std::vector<Channel> auxCh(auxTs.nChannels);
	for( int i = 0; i < auxTs.nChannels; i++ )
		auxCh[i].timeSeries = new Array<double>((size_t)auxTs.mtu.numberOfSamples);
		//auxCh[i].timeSeries = new Array<double>(5);
	
	*this = auxTs.mtu;
	//std::cout << auxTs.mtu.numberOfSamples << endl;
	if( infile.good() )
	{
		cout << "Reading time series data from " << input << endl << endl;
		get_data( auxCh );
		infile.close();
	}
	else
		cout << "Could not open " << input << endl << endl;

	auxTs.ch = auxCh;
	//for( int i = 0; i < auxTs.nChannels; i++ )
	//	delete auxCh[i].timeSeries;
	//return auxTs;
}

void ExtractorMTU::read_TSn_tag( ifstream &infile, StationFile *ts, size_t position )
{
	unsigned char *CurrentTag = new unsigned char[this->tagsize];
	double sampledenom, sampleenum;
	char sampleunit;

	if( infile.is_open() == true )
	{
		infile.seekg ( position );
		infile.read((char *) CurrentTag, 32);
		ts->date.startSecond = int( CurrentTag[0] );
		ts->date.startMinute = int( CurrentTag[1] );
		ts->date.startHour  = int( CurrentTag[2] );
		ts->date.startDay = int( CurrentTag[3] );
		ts->date.startMonth = int( CurrentTag[4] );
		ts->date.startYear = int( CurrentTag[5] );
		ts->date.startYear = (2000 + ts->date.startYear)*( ts->date.startYear < 70 ) + ( 1900 + ts->date.startYear)*( ts->date.startYear >= 70 );
		ts->mtu.nScansPerRecord = int( CurrentTag[11] * 256 + CurrentTag[10] );
		ts->nChannels = int( CurrentTag[12] );
		ts->mtu.tagLength = int( CurrentTag[13] );
		if( ts->mtu.tagLength != ts->mtu.tagsize )
			cout << "tagsize != taglength = " << ts->mtu.tagLength << "\n" << 
			"Change tagsize em mtu.h to " << ts->mtu.tagLength << "\n";
		ts->mtu.sampleLength = int( CurrentTag[17] );
		sampledenom = double( CurrentTag[19] * 256
			+ CurrentTag[18] );
		sampleunit = CurrentTag[20];
		sampleenum;
		switch (sampleunit)
		{
		case 0:
			sampleenum = 1.0;
			break;
		case 1:
			sampleenum = 60.0;
			break;
		case 2:
			sampleenum = 3600.0;
			break;
		case 3:
			sampleenum = 3600.0 * 24.0;
			break;
		}
		ts->samplingFrequency = (size_t)(sampledenom / sampleenum);
		ts->mtu.recordLength = (int)(ts->mtu.nScansPerRecord * ts->nChannels * ts->mtu.sampleLength + ts->mtu.tagLength);
	}
    delete[] CurrentTag;
}

void ExtractorMTU::fill_time_vector( Array<double> *timeVector, StationBase MtuBase, StationBase MtuCurrent, size_t idxOfTimeVector )
{
	//MtuBase.date.tsTime = boost::posix_time::time_from_string(MtuBase.date.getDateStr());
	//MtuCurrent.date.tsTime = boost::posix_time::time_from_string(MtuCurrent.date.getDateStr());

	//bool checkTime = MtuBase.date.tsTime == MtuCurrent.date.tsTime;
	//unsigned int val = 0;

	//if( checkTime == false )
	//{
	//	boost::posix_time::time_duration diffInTime = MtuCurrent.date.tsTime - MtuBase.date.tsTime;
	//	val = MtuBase.ts->samplingFrequency*diffInTime.total_seconds();
	//}

	//for( int i = 0; i < MtuCurrent.mtu->nScansPerRecord; i++ )
	//	(*timeVector)( i + idxOfTimeVector) = i + val;
}

Array<double> ExtractorMTU::get_mtu_time_vector( StationBase *station, string input )
{
	//StationBase newMtu = *station;
	//size_t startAt = 0;
	//size_t idxOfTimeVector = 0;

	//ifstream infile( input.c_str(), ios::binary );
	//Array<double> result( newMtu.mtu->numberOfSamples );
	//for( int i = 0; i < newMtu.mtu->numberOfRecords; i++ ) {
	//	read_TSn_tag( infile, &newMtu, startAt );
	//	startAt += newMtu.mtu->recordLength;
	//	fill_time_vector( &result, *station, newMtu, i*newMtu.mtu->nScansPerRecord );
	//}
	//infile.close();
	//return result;

	return Array<double>(1);
}

int ExtractorMTU::get_number_of_bytes( string infile )
{
  streampos begin,end;
  ifstream myfile ( infile, ios::binary );
  begin = myfile.tellg();
  myfile.seekg ( 0, ios::end );
  end = myfile.tellg();
  myfile.close();
  return end - begin;
}

void ExtractorMTU::get_data( std::vector<Channel> &auxCh )
{
	ifstream infile( this->tsnFile.c_str(), ios::binary );
	if ( infile.good() )
    {
		int counter = 0;
		//cout << this->numberOfRecords << endl;
		while (infile.good() && counter < this->numberOfRecords )
		{
			read_TSn_time_series( infile, auxCh, counter );
			counter++;
		}
    }
}

void ExtractorMTU::read_TSn_time_series( std::ifstream &infile, std::vector<Channel> &auxCh, size_t counter )
{
	unsigned char *CurrentTag = new unsigned char[this->tagsize];
	char *buffer, *currentbyte;

	infile.read( (char *) CurrentTag, this->tagsize );
	buffer = new char[ this->recordLength - this->tagLength ];
    infile.read( buffer, this->recordLength - this->tagLength );
    currentbyte = buffer;
	int line;

	// TODO - IMPROVE!

	// for nChannels = 5
			//cout << this->nScansPerRecord << endl;
	if( auxCh.size() == 5 ) {
		for (int i = 0; i < this->nScansPerRecord; i++ ) {
			line = i + counter*this->nScansPerRecord;
			auxCh[3].timeSeries[0](line) = read_value( currentbyte );
			currentbyte += this->sampleLength;
			auxCh[4].timeSeries[0](line) = read_value( currentbyte );
			currentbyte += this->sampleLength;
			auxCh[0].timeSeries[0](line) = read_value( currentbyte );
			currentbyte += this->sampleLength;
			auxCh[1].timeSeries[0](line) = read_value( currentbyte );
			currentbyte += this->sampleLength;
			auxCh[2].timeSeries[0](line) = read_value( currentbyte );
			currentbyte += this->sampleLength;
		}
	}
	// for nChannels = 4 - no Hz
	else if( auxCh.size() == 4 ) {
		for (int i = 0; i < this->nScansPerRecord; i++ ) {
			line = i + counter*this->nScansPerRecord;
			auxCh[2].timeSeries[0](line) = read_value( currentbyte );
			currentbyte += this->sampleLength;
			auxCh[3].timeSeries[0](line) = read_value( currentbyte );
			currentbyte += this->sampleLength;
			auxCh[0].timeSeries[0](line) = read_value( currentbyte );
			currentbyte += this->sampleLength;
			auxCh[1].timeSeries[0](line) = read_value( currentbyte );
			currentbyte += this->sampleLength;
		}
	}
	// for nChannels = 3 - just Hz
	else if( auxCh.size() == 4 ) {
		for (int i = 0; i < this->nScansPerRecord; i++ ) {
			line = i + counter*this->nScansPerRecord;
			auxCh[2].timeSeries[0](line) = read_value( currentbyte );
			currentbyte += this->sampleLength;
			auxCh[0].timeSeries[0](line) = read_value( currentbyte );
			currentbyte += this->sampleLength;
			auxCh[1].timeSeries[0](line) = read_value( currentbyte );
			currentbyte += this->sampleLength;
		}
	}
	else {
		std::cout << "number of channels different from 3, 4 or 5. Terminating the program" << std::endl;
		exit(77);
	}

    delete[] buffer;
    delete[] CurrentTag;
}

int ExtractorMTU::read_value( char *pos )
{
	unsigned char byte1, byte2, byte3;

	byte1 = *pos;
	byte2 = *(pos + 1);
	byte3 = *(pos + 2);
	return byte_to_int(byte1, byte2, byte3);
}

int ExtractorMTU::byte_to_int( unsigned char first, unsigned char second, unsigned char last )
{
	const char highbit = 0x80;
	int value;
	if (last & highbit)
    {
		last = ~last;
		second = ~second;
		first = ~first;
		value = -(last * 65536 + second * 256 + first + 1);
    }
	else
		value = last * 65536 + second * 256 + first;
	return (value);
}

table ExtractorMTU::read_tbl() {
    
	table tbl;
	std::string tblstr = this->tblFile;

    // open file
	ifstream tblf(tblstr);
    if(tblf.fail()) {
        cout<<"can't open the tbl file "<<tblstr<<endl;
        exit(1);
    }
    cout<<"opened tbl file "<<tblstr<<endl;

    // look for serial number of MTU unit
    if(positioning(tblf,"SNUM")) {
        tblf.read(reinterpret_cast <char *> (&tbl.snum),sizeof(tbl.snum));
        if(bswap)
            tbl.snum=bswap_32(tbl.snum);
    }
    else  // will try to find in the MTU-LR file
        tbl.snum=-1;
	
	// open file
	tblf.close();
	tblf.open(tblstr.c_str(), ios::in | ios::binary);
    if(tblf.fail()) {
        cout<<"can't open the tbl file "<<tblstr<<endl;
        exit(1);
    }
    cout<<"opened tbl file "<<tblstr<<endl;

    // look for number of samples per minute
    if(positioning(tblf,"SRPM")) {
        tblf.read(reinterpret_cast <char *> (&tbl.srpm),sizeof(tbl.srpm));
        if(bswap)
            tbl.srpm=bswap_32(tbl.srpm);
    }
    else  // will try to find in the MTU-LR file
        tbl.srpm=-1;

	// open file
	tblf.close();
    tblf.open(tblstr.c_str(), ios::in | ios::binary);
    if(tblf.fail()) {
        cout<<"can't open the tbl file "<<tblstr<<endl;
        exit(1);
    }
    cout<<"opened tbl file "<<tblstr<<endl;
    // read site ID of the measurement
    char  SITE[9];
    if(positioning(tblf, "SITE")) {
        tblf.read(SITE,9);
        tbl.site=SITE;
    }
    else {
        tbl.site="";
        cerr<<"can't find site ID of the measurement in tbl file "<<tblstr<<endl
            <<"will use the output file name as site ID"<<endl<<endl;
    }

	// open file
	tblf.close();
    tblf.open(tblstr.c_str(), ios::in | ios::binary);
    if(tblf.fail()) {
        cout<<"can't open the tbl file "<<tblstr<<endl;
        exit(1);
    }
    cout<<"opened tbl file "<<tblstr<<endl;
    // read survey ID
    char  SRVY[13];
    if(positioning(tblf, "SRVY")) {
        tblf.read(SRVY,13);
        tbl.srvy=SRVY;
    }
    else {
        tbl.srvy="";
        cerr<<"can't find site ID of the measurement in tbl file "<<tblstr<<endl;
    }
	
	// open file
	tblf.close();
    tblf.open(tblstr.c_str(), ios::in | ios::binary);
    if(tblf.fail()) {
        cout<<"can't open the tbl file "<<tblstr<<endl;
        exit(1);
    }
    cout<<"opened tbl file "<<tblstr<<endl;
    // read current data acquisition window
    if(positioning(tblf,"AWIN")) {
        tblf.read(reinterpret_cast <char *> (&tbl.awin),sizeof(tbl.awin));
        if(bswap)
            tbl.awin=bswap_32(tbl.awin);
    }
    else {
        tbl.awin=0;
        cerr<<"can't find data acquisition window in tbl file "<<tblstr<<endl;
    }
	
	// open file
	tblf.close();
    tblf.open(tblstr.c_str(), ios::in | ios::binary);
    if(tblf.fail()) {
        cout<<"can't open the tbl file "<<tblstr<<endl;
        exit(1);
    }
    cout<<"opened tbl file "<<tblstr<<endl;
    // read gain for "E" channels. 3,12,48
    if(positioning(tblf, "EGN")) {
        tblf.read(reinterpret_cast <char *> (&tbl.egn),sizeof(tbl.egn));
        if(bswap)
            tbl.egn=bswap_32(tbl.egn);
    }
    else {
        tbl.egn=3;
        cerr<<"can't find gain for E channels in tbl file "<<tblstr<<endl
            <<"will be set to 3"<<endl;
    }
	
	// open file
	tblf.close();
    tblf.open(tblstr.c_str(), ios::in | ios::binary);
    if(tblf.fail()) {
        cout<<"can't open the tbl file "<<tblstr<<endl;
        exit(1);
    }
    cout<<"opened tbl file "<<tblstr<<endl;
    // read gain for "H" channels. 2, 6, 8
    if(positioning(tblf, "HGN")) {
        tblf.read(reinterpret_cast <char *> (&tbl.hgn),sizeof(tbl.hgn));
        if(bswap)
            tbl.hgn=bswap_32(tbl.hgn);
    }
    else {
        tbl.hgn=2;
        cerr<<"can't find gain for H channels in tbl file "<<tblstr<<endl
            <<"will be set to 2"<<endl;
    }
	
	// open file
	tblf.close();
    tblf.open(tblstr.c_str(), ios::in | ios::binary);
    if(tblf.fail()) {
        cout<<"can't open the tbl file "<<tblstr<<endl;
        exit(1);
    }
    cout<<"opened tbl file "<<tblstr<<endl;
    // read Ex sensor azimuth
    if(positioning(tblf, "EAZM")) {
        tblf.read(reinterpret_cast <char *> (&tbl.eazm),sizeof(tbl.eazm));
        if(bswap)
            tbl.eazm=(float8)bswap_64((int8)tbl.eazm);
    }
    else {
        tbl.eazm=0;
        cerr<<"can't find Ex sensor azimuth in tbl file "<<tblstr<<endl
            <<"will be set to 0"<<endl;
    }
	
	// open file
	tblf.close();
    tblf.open(tblstr.c_str(), ios::in | ios::binary);
    if(tblf.fail()) {
        cout<<"can't open the tbl file "<<tblstr<<endl;
        exit(1);
    }
    cout<<"opened tbl file "<<tblstr<<endl;
    // read Ex dipole length, m
    if(positioning(tblf, "EXLN")) {
        tblf.read(reinterpret_cast <char *> (&tbl.exln),sizeof(tbl.exln));
        if(bswap)
            tbl.exln=(float8)bswap_64((int8)tbl.exln);
    }
    else {
        tbl.exln=1.;
        cerr<<"can't find Ex dipole length in tbl file "<<tblstr<<endl
            <<"will be set to 1 "
            <<"(it will work if using system calibration file)"<<endl;
    }
	
	// open file
	tblf.close();
    tblf.open(tblstr.c_str(), ios::in | ios::binary);
    if(tblf.fail()) {
        cout<<"can't open the tbl file "<<tblstr<<endl;
        exit(1);
    }
    cout<<"opened tbl file "<<tblstr<<endl;
    // read Ey dipole length, m
    if(positioning(tblf, "EYLN")) {
        tblf.read(reinterpret_cast <char *> (&tbl.eyln),sizeof(tbl.eyln));
        if(bswap)
            tbl.eyln=(float8)bswap_64((int8)tbl.eyln);
    }
    else {
        tbl.eyln=1.;
        cerr<<"can't find Ey dipole length in tbl file "<<tblstr<<endl
            <<"will be set to 1 "
            <<"(it will work if using system calibration file)"<<endl;
    }
	
	// open file
	tblf.close();
    tblf.open(tblstr.c_str(), ios::in | ios::binary);
    if(tblf.fail()) {
        cout<<"can't open the tbl file "<<tblstr<<endl;
        exit(1);
    }
    cout<<"opened tbl file "<<tblstr<<endl;
    // read Hx sensor azimuth
    if(positioning(tblf, "HAZM")) {
        tblf.read(reinterpret_cast <char *> (&tbl.hazm),sizeof(tbl.hazm));
        if(bswap)
            tbl.hazm=(float8)bswap_64((int8)tbl.hazm);
    }
    else {
        tbl.hazm=0;
        cerr<<"can't find Hx sensor azimuth in tbl file "<<tblstr<<endl
            <<"will be set to 0"<<endl;
    }
	
	// open file
	tblf.close();
    tblf.open(tblstr.c_str(), ios::in | ios::binary);
    if(tblf.fail()) {
        cout<<"can't open the tbl file "<<tblstr<<endl;
        exit(1);
    }
    cout<<"opened tbl file "<<tblstr<<endl;
    // read MAG sensor base line setting
    if(positioning(tblf, "BSEX")) {
        tblf.read(reinterpret_cast <char *> (&tbl.bsex),sizeof(tbl.bsex));
        if(bswap)
            tbl.bsex=bswap_32(tbl.bsex);
    }
    else {
        tbl.bsex=0;
        cerr<<"can't find Hx base line setting in tbl file "<<tblstr<<endl
            <<"will be set to 0"<<endl;
    }
	
	// open file
	tblf.close();
    tblf.open(tblstr.c_str(), ios::in | ios::binary);
    if(tblf.fail()) {
        cout<<"can't open the tbl file "<<tblstr<<endl;
        exit(1);
    }
    cout<<"opened tbl file "<<tblstr<<endl;
    // read MAG sensor base line setting
    if(positioning(tblf, "BSEY")) {
        tblf.read(reinterpret_cast <char *> (&tbl.bsey),sizeof(tbl.bsey));
        if(bswap)
            tbl.bsey=bswap_32(tbl.bsey);
    }
    else {
        tbl.bsey=0;
        cerr<<"can't find Hy base line setting in tbl file "<<tblstr<<endl
            <<"will be set to 0"<<endl;
    }
	
	// open file
	tblf.close();
    tblf.open(tblstr.c_str(), ios::in | ios::binary);
    if(tblf.fail()) {
        cout<<"can't open the tbl file "<<tblstr<<endl;
        exit(1);
    }
    cout<<"opened tbl file "<<tblstr<<endl;
    // read MAG sensor base line setting
    if(positioning(tblf, "BSEZ")) {
        tblf.read(reinterpret_cast <char *> (&tbl.bsez),sizeof(tbl.bsez));
        if(bswap)
            tbl.bsez=bswap_32(tbl.bsez);
    }
    else {
        tbl.bsez=0;
        cerr<<"can't find Hz base line setting in tbl file "<<tblstr<<endl
            <<"will be set to 0"<<endl;
    }
	
	// open file
	tblf.close();
    tblf.open(tblstr.c_str(), ios::in | ios::binary);
    if(tblf.fail()) {
        cout<<"can't open the tbl file "<<tblstr<<endl;
        exit(1);
    }
    cout<<"opened tbl file "<<tblstr<<endl;
    // read elevation altitude of site [m]
    int  ELEV;
    if(positioning(tblf, "ELEV")) {
        tblf.read(reinterpret_cast <char *> (&tbl.elev),sizeof(tbl.elev));
        if(bswap)
            tbl.elev=bswap_32(tbl.elev);
    }
    else {
        tbl.elev=0;
        cerr<<"can't find elevation altitude of site in tbl file "<<tblstr<<endl
            <<"will be set to 0"<<endl;
    }
	
	// open file
	tblf.close();
    tblf.open(tblstr.c_str(), ios::in | ios::binary);
    if(tblf.fail()) {
        cout<<"can't open the tbl file "<<tblstr<<endl;
        exit(1);
    }
    cout<<"opened tbl file "<<tblstr<<endl;
    // read latitude
    char  LATG[13];
    if(positioning(tblf, "LATG")) {
        tblf.read(LATG,13);
        tbl.latg=conv_lat(LATG);
    }
    else {
        cerr<<"can't find latitude in tbl file "<<tblstr<<endl
            <<"enter latitude in degrees (north positive):"<<endl;
        cin>>tbl.latg;
    }
	
	// open file
	tblf.close();
    tblf.open(tblstr.c_str(), ios::in | ios::binary);
    if(tblf.fail()) {
        cout<<"can't open the tbl file "<<tblstr<<endl;
        exit(1);
    }
    cout<<"opened tbl file "<<tblstr<<endl;
    // read longitude
    char  LNGG[13];
    if(positioning(tblf, "LNGG")) {
        tblf.read(LNGG,13);
        tbl.lngg=conv_lon(LNGG);
    }
    else {
        cerr<<"can't find longitude in tbl file "<<tblstr<<endl
            <<"enter longitude in degrees (east positive):"<<endl;
        cin>>tbl.lngg;
    }
	
	GetDeclination get_decl;
	// open file
	tblf.close();
    tblf.open(tblstr.c_str(), ios::in | ios::binary);
    if(tblf.fail()) {
        cout<<"can't open the tbl file "<<tblstr<<endl;
        exit(1);
    }
    cout<<"opened tbl file "<<tblstr<<endl;
    // calculate declination. Will get STIM as a reference date
    if(positioning(tblf, "STIM")) {
        amxtds  STIM;
        tblf.read(reinterpret_cast <char *> (&STIM),sizeof(amxtds));
        tbl.decl = get_decl.declination(tbl.latg, tbl.lngg, tbl.elev/1e3, STIM.yr+STIM.cen*100, STIM.mn, STIM.dy);
    }
    else {
        int year, month, day;
        cout<<"enter a reference year month day to calculate declination"<<endl
            <<"(e.g.: 2001 08 12)"<<endl;
        cin>>year>>month>>day;
        tbl.decl = get_decl.declination(tbl.latg, tbl.lngg, tbl.elev/1e3, year, month, day);
    }
	if (tbl.decl == 0)
		tbl.decl = get_decl.declination_input;
	
	// open file
	tblf.close();
    tblf.open(tblstr.c_str(), ios::in | ios::binary);
    if(tblf.fail()) {
        cout<<"can't open the tbl file "<<tblstr<<endl;
        exit(1);
    }
    cout<<"opened tbl file "<<tblstr<<endl;
    // read coil nominal gain, mV/nT
    if(positioning(tblf, "HNOM")) {
        tblf.read(reinterpret_cast <char *> (&tbl.hnom),sizeof(tbl.hnom));
        if(bswap)
            tbl.hnom=(float8)bswap_64((int8)tbl.hnom);
    }
    else {
        tbl.hnom=4.;
        cerr<<"can't find coil nominal gain in tbl file "<<tblstr<<endl
            <<"will be set to 4 mV/nT"<<endl;
    }
    tblf.close();

    return tbl;
}

bool ExtractorMTU::correct_types(void) {
    if(sizeof(INT2)==2 && sizeof(INT4)==4 && sizeof(int8)==8 && sizeof(float8)==8)
        return true;
    else {
        cout<<"some type haven't the correct size:"<<endl;
        cout<<"size of int2       "<<setw(4)<<sizeof(INT2)<<endl;
        cout<<"size of int4       "<<setw(4)<<sizeof(INT4)<<endl;
        cout<<"size of float8     "<<setw(4)<<sizeof(float8)<<endl<<endl;

        cout<<"you will need to edit the bytes.h file."<<endl;
        cout<<"the sizes of types in your machine are:"<<endl;
        cout<<"size of short      "<<setw(4)<<sizeof(short)<<endl;
        cout<<"size of int        "<<setw(4)<<sizeof(int)<<endl;
        cout<<"size of long       "<<setw(4)<<sizeof(long)<<endl;
        cout<<"size of long long  "<<setw(4)<<sizeof(long long)<<endl;
        cout<<"size of float      "<<setw(4)<<sizeof(float)<<endl;
        cout<<"size of double     "<<setw(4)<<sizeof(double)<<endl;
        cout<<"size of long double"<<setw(4)<<sizeof(long double)<<endl;
        return false;
    }
}

INT2 ExtractorMTU::bswap_16(INT2 datum) {
    return (((datum & 0xff00) >>  8) | ((datum & 0x00ff) <<  8));
}

INT4 ExtractorMTU::bswap_32(INT4 datum) {
    return (((datum & 0xff000000) >> 24) | ((datum & 0x00ff0000) >>  8) | \
            ((datum & 0x0000ff00) <<  8) | ((datum & 0x000000ff) << 24));
}

int8 ExtractorMTU::bswap_64(int8 datum) {
    INT4 lword= INT4(datum & 0x00000000ffffffff);
    INT4 hword=INT4(datum>>32);
    int8 lw=bswap_32(lword);
    int8 hw=bswap_32(hword);
    return (int8)(hw<<32 | lw);
}

bool ExtractorMTU::positioning(ifstream& tblf, char* label) {
    tblf.seekg(0);

    char labeltmp[5];

    // look for label
    tblf.read(reinterpret_cast <char *> (labeltmp),5);
    while(strcmp(labeltmp,label)!=0 && !tblf.eof()) {
        tblf.seekg(20,ios::cur); // skip to the next label
        tblf.read(labeltmp,5);   // and read it
    }

    // if found label, positioning to get it
    if(strcmp(labeltmp,label)==0) {
        tblf.seekg(7,ios::cur);
        return true;
    }
    else
        return false;
}

double ExtractorMTU::conv_lat(const char* LATG) {
    string latstr(LATG,13);

    string degstr(latstr,0,2);
    int deg=atoi(degstr.c_str());

    string::size_type coma=latstr.find(",");
    string minstr(latstr,2,coma-2);
    double min=atof(minstr.c_str());

    double signal=1.;
    if(latstr.at(coma+1)=='S')
        signal*=-1.;

    return signal*(deg+min/60.);
}

double ExtractorMTU::conv_lon(const char* LNGG) {
    string lonstr(LNGG,13);

    string degstr(lonstr,0,3);
    int deg=atoi(degstr.c_str());

    string::size_type coma=lonstr.find(",");
    string minstr(lonstr,3,coma-3);
    double min=atof(minstr.c_str());

    double signal=1.;
    if(lonstr.at(coma+1)=='W')
        signal*=-1.;

    return signal*(deg+min/60.);
}

Array<Complex> ExtractorMTU::read_cts_file( StationBase *station, string ctsFileName )
{
	//ifstream in(ctsFileName);
	//if(!in) {
 //       cerr<<"Could not open CTS file \""<<ctsFileName<<'\"'<<endl;
 //       exit(1);
 //   }
	//    
	//string line;

	//// detect beginning of data in cts file
 //   while(line.find( '/' ) == string::npos)
	//	getline(in,line);
	//
	//size_t numberOfPoints = 0;
	//while(in.eof() == false) {
	//	getline(in,line);
	//	numberOfPoints++;
	//}
	//in.close();

	//Array<Complex> result( numberOfPoints-1,  station->numberOfChannels + 2 );
	//in.open(ctsFileName);
	//if(!in) {
 //       cerr<<"Could not open CTS file \""<<ctsFileName<<'\"'<<endl;
 //       exit(1);
 //   }

	//// detect beginning of data in cts file
 //   while(line.find( '/' ) == string::npos)
	//	getline(in,line);

	//// fills in result with correction values
	//for( int i = 0; i < result.GetDescriptor().GetDim(0); i++ ) {
	//	getline(in,line);
	//	cts2array( &result, line, i );		
	//}

	//return result;

	return 0;
}

void ExtractorMTU::cts2array( Array<Complex> *data, string line, index_t idx )
{
	std::stringstream istr(line.c_str());
    double freq, real, imag;

    istr >> (*data)( idx, 0 );
	if (istr.peek() == ',')
			istr.ignore();

    istr >> (*data)( idx, 1 );
	if (istr.peek() == ',')
			istr.ignore();

	for( int i = 2; i < data->GetDescriptor().GetDim(1); i++ ) {
		if (istr.peek() == ',')
			istr.ignore();
		istr >> real;
		if (istr.peek() == ',')
			istr.ignore();
		istr >> imag;
		(*data)( idx, i ) = Complex( real, imag );
	}
}

phoenixTsBand_t ExtractorMTU::get_phoenix_TS_band( string fileName )
{
	if( fileName.find("TS2") != std::string::npos )
		return TS2;
	else if( fileName.find("TS3") != std::string::npos )
		return TS3;
	else if( fileName.find("TS4") != std::string::npos )
		return TS4;
	else if( fileName.find("TS5") != std::string::npos )
		return TS5;
	else
		return noTS;
}
