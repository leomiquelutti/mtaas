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

#include "boost/date_time/posix_time/posix_time.hpp"

#define PI 3.141592653589793238462643
#define pi 3.141592653589793238462643

#ifndef MAXPATH
#define MAXPATH 256
#endif

#ifdef WORDS_BIGENDIAN
    const bool bswap=true;
#else
    const bool bswap=false;
#endif

#define MININUM_OF_SEGMENTS_TO_COHERENCE_PRESELECTION 25

typedef short INT2;
typedef int INT4;
typedef long long int8;
typedef double float8;
typedef std::complex<double> Complex;

using namespace matCUDA;

class StationBase;

typedef struct {                   /* AMX Time/Date Structure             */
      unsigned char sec;           /* seconds  (0-59)                     */
      unsigned char min;           /* minutes  (0-59)                     */
      unsigned char hr;            /* hours    (0-23)                     */
      unsigned char dy;            /* day      (1-31)                     */
      unsigned char mn;            /* month    (1-12)                     */
      unsigned char yr;            /* year     (0-99)                     */
      unsigned char ow;            /* day of week (Mon=1 to Sun=7)        */
      unsigned char cen;           /* 0 if time/date is incorrect         */
                                   /* century if time/date is correct     */
} amxtds;

typedef struct {

    int   snum;          // Serial Number of MTU unit
    int   srpm;          // Number of samples per minute
    std::string site;          // Site ID of the measurement
    std::string srvy;          // Survey ID
    int   awin;          // Current data acquisition window
    int   egn;           // Gain for "E" channels. 3,12,48
    int   hgn;           // Gain for "H" channels. 2, 6, 8
    float8 eazm;          // Ex sensor azimuth
    float8 exln;          // Ex dipole length, m
    float8 eyln;          // Ey dipole length, m
    float8 hazm;          // Hx sensor azimuth
    int   bsex;          // MAG sensor base line setting
    int   bsey;          // MAG sensor base line setting
    int   bsez;          // MAG sensor base line setting
    int   elev;          // Elevation altitude of site [m]
    float  decl;          // magnetic declination
    double latg;          // Latitude
    double lngg;          // longitude
    float8 hnom;          // Coil nominal gain, mV/nT
    short  fieldtype;     // =0 for box calibration, =1 for system calibration
    unsigned short nchan; // number of channels
	int   l3ns;		  // number of seconds acquired continuously of band TS3
	int   l4ns;		  // number of seconds acquired continuously of band TS4
	int   v5sr;		  // sampling frequency for bands
} table;

typedef enum{
	noTS = -1,
    TS2  = 2,
    TS3  = 3,
    TS4  = 4,
    TS5  = 5
} phoenixTsBand_t;

typedef enum {
	FILE_TYPE_EMI = 1,
	FILE_TYPE_LMT,
	FILE_TYPE_SIO,
	FILE_TYPE_MTU,
	FILE_TYPE_ADU_07
} FILE_TYPE;

class DirectoryProperties
{
public:
	
	size_t						numberOfFolders;
	std::string					pathName;

	void						define_files( const std::string inputSubPath, size_t numberOfFiles, std::string *files );
	std::vector<StationBase>	initialize_stations();
	//StationBase*				initialize_stations();
	size_t						get_number_of_tbl_files( const std::string inputSubPath );
	void						get_tbl_names( std::string *TBLs, std::string file, size_t size );

private:
	
	void						fill_station_names( StationBase *station );
	void						fill_station_names( std::vector<StationBase> *station );
	size_t						get_number_of_subfolders( std::string inputPath );
};

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

//template<typename TElement> 
class MtuExtractor
	//:public Array<TElement>
{
	//template<typename TElement> friend class Array;
public:
	MtuExtractor():
		empirical_factor(0.1),
		ECV(1.E6*empirical_factor),
		HCV(1.E9),
		FS(0x7FFFFF),
		tagsize(32){}

	// variables
	std::string		inputTbl;
	phoenixTsBand_t MtuTsBand;
	table			tbl;
	size_t			amountOfTbl;

	// functions
	Array<double>	get_mtu_time_vector( StationBase *station, std::string input );
	phoenixTsBand_t get_phoenix_TS_band( std::string fileName );
	Array<Complex>	read_cts_file( StationBase *station, std::string ctsFileName );
	Array<double>	read_mtu_data( StationBase *station, std::string input );
	table			read_tbl(std::string tblstr);

private:

	// variables
	std::string		ctsFileName;
	const double	ECV;// = 1.E6*empirical_factor;    // convert V/m to mV/km
	const double	empirical_factor;// = 0.1;
	const long		FS;// = 0x7FFFFF;  // full scale (normalizing the ADN)
	const double	HCV;// = 1.E9;     // convert T to nT
	int				nScansPerRecord;
	int				numberOfBytes;
	int				numberOfCtsFiles;
	int				nChannels;
	int				numberOfRecords;
	int				numberOfSamples;
	int				numberOfTsFilesForCurrentCts;
	int				recordLength;
	int				sampleLength;
	//double		samplingRate;
	int				tagLength;
	const size_t	tagsize;// = 32;
	
	// functions
	static INT2		bswap_16(INT2 datum);
	static INT4		bswap_32(INT4 datum);
	static int8		bswap_64(int8 datum);
	static int		byte_to_int( unsigned char first, unsigned char second, unsigned char last );
	static double	conv_lat(const char* LATG);
	static double	conv_lon(const char* LNGG);
	static bool		correct_types(void);
	static void		cts2array( Array<Complex> *, std::string, index_t );
	static void		fill_time_vector( Array<double> *timeVector, StationBase MtuBase, StationBase MtuCurrent, size_t idxOfTimeVector );
	static void		get_data( const std::string filename, Array<double> *data, StationBase *station );
	static int		get_number_of_bytes( std::string infile );
	static bool		positioning( std::ifstream &tbl, char* label);
	static void		read_TSn_tag( std::ifstream &filestream, StationBase *Mtu, int position );
	static void		read_TSn_time_series( std::ifstream &infile, int counter, Array<double> *data, StationBase *station );
	static int		read_value( char *pos );
};

class StationBase : public MtuExtractor
{
public:
	
	StationBase() {};
	StationBase(const StationBase& element) {*this = element;};
    StationBase& operator = (const StationBase& element);     //which needs definition
	std::string		fileName;
	std::string		stationName;
	FILE_TYPE		fileType;
	Date			date;
	Position		position;

	// raw data
	Array<double>	*data;
	Array<Complex>	*correctionToData;
	Array<double>	*timeVector;

	// parameters for processing
	size_t			numberOfChannels;
	size_t			samplingRateInHertz;
	double			maxFreqInHertz;
	double			minFreqInHertz;
	bool			isContinuous;
	unsigned int	continuousSamples;

	// functions
	void			initialize_parameters();
	bool			is_acquisition_continuos( Array<double> *timeVector );

private:
	double			exDipoleLength;
	double			eyDipoleLength;
};

#endif

//// USELESS

//class Header
//{
//public:
//	Header();
//	~Header();
//
//public:
//	char			filename[MAXPATH]; // file's name which is being processed
//	FILE_TYPE		filetype; // file's type - EMI, LEMI or SIO
//	long double		datastart; // which time acquisition started
//	int				nfiles; // number of files being processed
//	int				q;
//	char			datdir_name[MAXPATH]; // NOT BEING USED
//	char			sendir_name[MAXPATH]; // NOT BEING USED
//	int				say; // NOT BEING USED
//	int				fhi; // number of frequencies
//	Array<double>	&freq; // frequencies evaluated
//	Array<double>	&freqall; // all frequencies evaluated - for bands and decimations case
//	int				fmax; // maximum number of frequencies
//	int				cmax; // NOT BEING USED
//	double			satmx; // maximum value for saturation checking
//	Array<double>	&vpsmpl; // volt per sample
//	int				ncpblk; // number of crosspower blocks used to stack
//	Array<double>	&par; // calculated parameters
//	Array<double>	&parerr; // errors of calculated parameters
//	Array<double>	&parall; // all calculated parameters - for bands and decimations case
//	Array<double>	&parerrall; // all errors of calculated parameters - for bands and decimations case
//	Array<int>		&chmt; // chanels  identifiers
//	int				ropt; // for rotation options - see 'set_rot' function
//	int				strk; // strike angle information - for use this -> ropt = 4;
//	int				deci; // decimation level
//	char			efrsp[8][32]; // name of sensors files used for hardware response correction
//	Array<double>	&eh; // chanels order information
//	Array<int>		&pos; // NOT BEING USED
//	int				pos2; // used for reading lemi time series
//	int				meanon; // NOT BEING USED
//	int				stdlim; // NOT BEING USED
//	Array<double>	&pltorder; // NOT BEING USED
//	Array<int>		&dateform; // NOT BEING USED
//	int				window; // NOT BEING USED
//	int				unit_type; // NOT BEING USED
//	Array<double>	&tlenght; // NOT BEING USED
//	char			sitenam[128]; // site name
//	char			header[32][128]; // header (char form)
//	char			chid[8][8]; // chanels identifier
//	Array<double>	&indy; // NOT BEING USED
//	int				nosite; // number of sites - NOT BEING USED
//	int				nochnl; // number of chanels
//	int				pwline; // powerline frequency - NOT BEING USED
//	Array<double>	&efpar; // acquisition parameters, such as angle between dipoles, it's lenght, gain, band(HF/LF), filter, resistance and spontaneous potencial
//	Array<int>		&gain; // gain
//	int				tppb; // header information about windows acquired
//	int				tnblk; // header information about windows acquired
//	double			hpfv; // frequency lower value
//	double			lpfv; // frequency higher value
//	double			smplfq; // sample acquisition frequency
//	Array<int>		&filech; // files chanels - NOT	BEING USED
//	Array<int>		&datach; // data chanels - NOT	BEING USED - same as 'filech'
//	Array<double>	&sn; // some header information
//	char			strngs[32][32]; // NOT BEING USED
//	char			cohnam[16][8]; // NOT BEING USED
//	int				nfq; // number of frequencies in decimation/band
//	int				f1; // first frequency in decimation/band
//	Array<double>	&cavg; // crosspower average
//	Array<double>	&navg; // number of crosspower blocks for average
//	Array<double>	&fqs; // NOT BEING USED
//	Array<double>	&bws; // bandwidth
//	Array<double>	&nfcs;
//	Array<double>	&navs; // NOT BEING USED
//	int				lwind;
//	Array<double>	&w;
//	Array<double>	&chsat; // chanels saturation
//	int				lngth; // windows lngth ( number of points)
//	Array<double>	&time; // store time acquisition information
//	Array<double>	&fqavg; // frequency average
//	Array<double>	&bwavg; // bandwidth average
//	Array<int>		&ifq; // index of frequencies
//	Array<double>	&nfc; // number of frequency coefficients
//	Array<double>	&bw; // NOT BEING USED
//	char			acqby[32]; // who acquired
//	double			freqrate; // NOT BEING USED
//	char			chanid[8][8]; // chanel identifier
//	Array<double>	&chnuse; 
//	double			plt_start;
//	int				directory1;
//	double			fs; // SIO header info
//	int				type; // SIO header info
//	int				nch; // number of chanels - SIO
//	char			headername[MAXPATH]; // NOT BEING USED
//};
//
//class FileExtractor : StationBase
//{
//public:
//
//	// UNUSED
//	int readheader_emi(Header *files, int cont);
//	int readheader_lemi(Header *files, int cont, int idec);
//	void readheader_sio(Header *files, int cont);
//	int localiza(Header *files, int cont);
//	int readts_emi(Header *files, int cont,  Array<double> &data, int ijk, int icp, int start);
//	void readts_lemi(Header *files, int cont, Array<double> &data, int icp);
//	int readts_sio(Header *files, int cont, Array<double> &data, int icp);
//	void printresults(Header *files, int cont, Array<double> &freqall, Array<double> &parall, Array<double> &parerrall, int index);
//	void printresults_robust(Header *files, int cont, Array<double> &freqall, Array<double> &parall, Array<double> &parerrall, int index);
//	void edi(Header *files, int cont, int index, Array<double> &freqall, Array<double> &parall, Array<double> &parerrall);
//	void read_sys_cor_emi(Header *files, int cont, Array<double> &ampf, Array<double> &phsf, Array<int> &rspfile);
//
//private:
//
//	// UNUSED
//	void process_line(char *line, Header *files, int cont);
//	void process_data_file(char *line, Header *files, int cont);
//	void process_version_id(char *line, Header *files, int cont);
//	void process_date_time(char *line, Header *files, int cont);
//	void process_survey_id(char *line, Header *files, int cont);
//	void process_survey_co(char *line, Header *files, int cont);
//	void process_client_co(char *line, Header *files, int cont);
//	void process_operator(char *line, Header *files, int cont);
//	void process_area(char *line, Header *files, int cont);
//	void process_remote(char *line, Header *files, int cont);
//	void process_weather(char *line, Header *files, int cont);
//	void process_si(char *line, Header *files, int cont);
//	void process_ch(char *line, Header *files, int cont);
//	void process_ext_clock(char *line, Header *files, int cont);
//	void process_powerline(char *line, Header *files, int cont);
//	void process_site(char *line, Header *files, int cont);
//	void process_chanel(char *line, Header *files, int cont);
//	void process_chaext(char *line, Header *files, int cont);
//	void process_parameter(char *line, Header *files, int cont);
//	void calculate_sn(Header *files, int cont);
//	void print_header(Header *files, int cont);
//	void print_header_lemi(Header *files, int cont);
//	void process_line_bin(char *line, Header *file, int cont);
//	void process_date_time_bin(char *line, Header *files, int cont);
//	void process_chanel_bin(char *line, Header *files, int cont);
//	void process_chaext_bin(char *line, Header *files, int cont);
//	void process_data_type(char *line, Header *files, int cont);
//	void process_parameter_bin(char *line, Header *files, int cont);
//	void process_filebin(FILE *fp, Header *files, int cont);
//	void calculate_sn_bin(Header *files, int cont);
//	void process_chnuse(FILE *fp, Header *files, int cont);
//	void calculate_satmx(Header *files, int cont);
//	void print_header_bin(Header *files, int cont);
//	void read_sys_cor_sio(Header *files, int cont, Array<double> &ampf, Array<double> &phsf, Array<int> &rspfile);
//	void getdata(Header *files, Array<double> &datag, long int startblk, int nblks, int cont);
//};
//
//class TransformerBase
//{
//};
//
//class CrossPowerTransformer : TransformerBase
//{
//public:
//	void cp_add(Header *files, int cont, Array<double> &c);
//	void cpcrrct(Header *files, int cont, int f1, Array<double> &c);
//	void crossm(Header *files, int cont, Array<double> &d, Array<double> &cp); 
//	void fillInFftValues(Header *files, int cont, Array<double> &d, int icp, Array<Complex> &fftValues);
//	void sys_cor(Header *files, int cont, Array<double> &c, Array<double> &ampf, Array<double> &phsf);
//};
//
//class FFTTransformer : TransformerBase
//{
//public:
//	void freqselect(Header *files, int cont);
//	void sat_check(Header *files, int cont, Array<double> &dat); 
//	void fft9(Header *files, int cont,  Array<double> &d); 
//
//private:
//	void hammingm(Array<double> &x, int ltrf, double lwind, Array<double> &w); 
//	void fftm(Array<double> &x, Array<double> &C, int n2);
//	double round(double d);
//};
//
//class RobustTransformer : TransformerBase
//{
//public:
//	void multipleCoherences(Header *files, int cont, Array<Complex> &fftValues, Array<double> &multiCoh, int icp);
//	void robust_processing(Header *files, int cont, Array<Complex> &fftValues, Array<double> &multiCoh,Array<double> &ampf,Array<double> &phsf, int icp);
//
//private:
//	int decide_size_sn_and_vpsmpl(Header *files, int cont);
//	void fills_in_sn_and_vpsmpl(Header *files, int cont, int counter, Array<double> &sn, Array<double> &vpsmpl);
//	void cohbul(Header *files, int cont, Array<double> &multiCoh, Array<int> &wtCoh );
//	void fills_in_wt( Header *files, int cont, Array<int> &wtCoh, Array<double> &wt, int frequencyNumber );
//	void fills_in_aa( Header *files, int cont, Array<Complex> &aa, Array<double> &ampf, Array<double> &phsf, int index);
//	void fills_in_chanel_values(Header *files, int cont, Array<Complex> &aa, Array<Complex> &fftValues, Array<Complex> &chanels, Array<double> &fct, int frequencyNumber );
//	void fills_in_fct(Header *files, int cont, Array<double> &fct);
//	void impedance_and_tipper_calculus(Header *files, int cont, int ex, int ey, int hx, int hy, int hz, int rx, int ry, Array<Complex> &chanels, Array<double> &wt, int frequencyNumber, Array<Complex> &Z );
//	Complex impedance_calculus( Array<Complex> &chanels, Array<double> &wt, int num1, int num2, int num3, int num4, int denom1, int denom2, int denom3, int denom4, int numberOfEvents, int wt_identifier, int icp );
//	void tipper_calculus( Array<Complex> &chanels, Array<double> &wt, int hx, int hy, int hz, int rx, int ry, int numberOfEvents, int icp, Array<Complex> &T );
//	void residual_calculus( Header *files, int cont, Array<Complex> &chanels, Array<double> &wt, int ex, int ey, int hx, int hy, int hz, int frequencyNumber, Array<Complex> &Z, Array<double> &residual, Array<double> &sumResidual );
//	void hat_weight_calculus( Header *files, int cont, Array<Complex> &chanels, Array<double> &wt, Array<double> &wtHatWeight, int ex, int ey, int hx, int hy, int hz, int rx, int ry, int frequencyNumber, Array<Complex> &Z );
//	void huber_weight( Header *files, int cont, Array<Complex> &chanels, Array<double> &wt, Array<double> &wtHatWeight, Array<double> &residual, int frequencyNumber, int ex, int ey, int hx, int hy, int hz, int rx, int ry, Array<Complex> &Z );
//	void fills_in_E_h_adjustFactor( Header *files, int cont, Array<Complex> A, Array<Complex> R, Array<Complex> E, Array<double> h, Array<double> adjustFactor, int frequencyNumber );
//};
//
//class LoaderBase
//{
//};
//
//class FileLoader : LoaderBase
//{
//public:
//	double results(Header *files, int cont, int icp, Array<double> &dis);
//	void results_robust(Header *files, int cont, Array<double> &z);
//
//private:
//	double szcalc(double dfdxr, double dfdxi, double dfdyr, double dfdyi, double dfdan, double c2r, double c2i, double tyr, double tyi, double txr, double txi, Array<double> &a, double com);
//	double fix_dir(double angl, int erot);
//	double minf(double n1, double n2);
//	void set_rot(Array<double> &par, Array<double> &parerr, Array<double> &freq, int nfre, int ropt, int erot, double strk); 
//	void multcoh(double *coh12, double *coh1y, double *coh2y, double *mcoh, Array<double> &c, int C, int ich1, int ich2, int indx);
//	void ff (Header *files, int cont, double zxxr, double zyyr, double zxxi, double zyyi, double zxyr, double zyxr, double zxyi, double zyxi, double an, Array<double> &a, Array<double> &n, Array<double> &nav, int k, Array<double> &f);
//	double fnsum (Array<double> &e, Array<double> &a, Array<double> &g);
//	void zrot(double *rxyr, double *rxyi, double *ryxr, double *ryxi, Array<double> &par, double theta, int ii);
//};
//
//class Util
//{
//public:
//	static void initstr(Header *files, int cont);
//	static void zerarparametros(Header *files, int cont);
//	static int char2int32(char *cbuffer);
//	static int char2int16(char *cbuffer);
//	static double mindbv(Array<double> &vector, int dim);
//	static double maxdbv(Array<double> &vector, int dim);
//	static double mindb (double n1, double n2); 
//	static void copy_ampf_and_phsf(Header *files, int cont, Array<double> &ampf, Array<double> &phsf, Array<double> &ampf_new, Array<double> &phsf_new);
//	static void copyfreqselect(Header *files, int cont);
//	static int finalresults(Header *files, int cont, Array<double> &freqall, Array<double> &parall, Array<double> &parerrall);
//	static void impedance(Array<double> &dis);
//	static Complex conjugate( Complex a );
//	static void invert_2x2_complex_matrix( Array<Complex> &matrix );
//	static void qr_decomposition_via_GramSchimidt( Array<Complex> &A, Array<Complex> &Q, Array<Complex> &R );
//	static int get_icp( Header *files, int cont );
//	static int get_number_of_events( Header *files, int cont, int frequencyNumber );
//	static int getInicialPositionOfEvent( Header *files, int cont, int frequencyNumber );
//	static void print_complex_array( Array<Complex> &A );
//	static void print_double_array( Array<double> &A );
//
//private:
//	double minfv(Array<double> &vector); 
//	double maxfv(Array<double> &vector);
//};