#ifndef MTU_H
#define MTU_H

#include "stdafx.h"

#ifdef WORDS_BIGENDIAN
    const bool bswap=true;
#else
    const bool bswap=false;
#endif

using namespace std;

class StationBase;

typedef short INT2;
typedef int INT4;
typedef long long int8;
typedef double float8;

// auxiliary for ExtractorMTU
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

// auxiliary for ExtractorMTU
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

// auxiliary for ExtractorMTU
typedef enum{
	noTS = -1,
    TS2  = 2,
    TS3  = 3,
    TS4  = 4,
    TS5  = 5
} phoenixTsBand_t;

// class for reading MTU data
class ExtractorMTU
{
public:
	ExtractorMTU():
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
	static void read_time_series( StationBase *station );
	matCUDA::Array<double>	get_mtu_time_vector( StationBase *station, std::string input );
	phoenixTsBand_t get_phoenix_TS_band( std::string fileName );
	matCUDA::Array<Complex>	read_cts_file( StationBase *station, std::string ctsFileName );
	matCUDA::Array<double>	*read_mtu_data( StationBase *station, std::string input );
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
	static void		cts2array( matCUDA::Array<Complex> *, std::string, index_t );
	static void		fill_time_vector( matCUDA::Array<double> *timeVector, StationBase MtuBase, StationBase MtuCurrent, size_t idxOfTimeVector );
	static void		get_data( const std::string filename, matCUDA::Array<double> *data, StationBase *station );
	static int		get_number_of_bytes( std::string infile );
	static bool		positioning( std::ifstream &tbl, char* label);
	static void		read_TSn_tag( std::ifstream &filestream, StationBase *mtu, int position );
	static void		read_TSn_time_series( std::ifstream &infile, int counter, matCUDA::Array<double> *data, StationBase *station );
	static int		read_value( char *pos );
};

#endif // MTU_H