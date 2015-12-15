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
class StationFile;
class Channel;
class DirectoryProperties;

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
		empirical_factor(1/8),
		tagsize(32),
		FS(0x7FFFFF) {
			ECV = 1.E6*empirical_factor;    // convert V/m to mV/km
			HCV = 1.E9;     // convert T to nT
		}

	ExtractorMTU(const ExtractorMTU& element) {*this = element;};

    ExtractorMTU& operator = (const ExtractorMTU& element);     //which needs definition

	// paths to files
	std::string		tblFile;
	std::string		tsnFile;
	std::string		ctsFile;

	// data
	phoenixTsBand_t mtuTsBand;
	table			tbl;

	// functions
	static void				read_time_series( StationBase *station, DirectoryProperties *dirInfo );
	matCUDA::Array<double>	get_mtu_time_vector( StationBase *station, std::string input );
	phoenixTsBand_t			get_phoenix_TS_band( std::string fileName );
	void					read_mtu_data( StationFile &auxTs  );
	StationFile				get_parameters( StationFile auxTs );
	table					read_tbl();

private:

	// variables
	double			empirical_factor; 
	double			ECV;    // convert V/m to mV/km
	double			HCV;     // convert T to nT
	long			FS;  // full scale (normalizing the ADN)
	int				nScansPerRecord;
	int				numberOfBytes;
	int				numberOfRecords;
	int				numberOfSamples;
	int				recordLength;
	int				sampleLength;
	int				tagLength;
	size_t			tagsize;// = 32;
	
	// functions
	static INT2		bswap_16(INT2 datum);
	static INT4		bswap_32(INT4 datum);
	static int8		bswap_64(int8 datum);
	static int		byte_to_int( unsigned char first, unsigned char second, unsigned char last );
	static double	conv_lat(const char* LATG);
	static double	conv_lon(const char* LNGG);
	static bool		correct_types(void);
	static void		cts2array( std::vector<Channel> &, std::string, index_t );
	static void		fill_time_vector( matCUDA::Array<double> *timeVector, StationBase MtuBase, StationBase MtuCurrent, size_t idxOfTimeVector );
	void			get_data( std::vector<Channel> &auxCh );
	static int		get_number_of_bytes( std::string infile );
	static bool		positioning( std::ifstream &tbl, char* label);
	void			read_TSn_tag( ifstream &infile, StationFile *ts, size_t position );
	void 			read_TSn_time_series( std::ifstream &infile, std::vector<Channel> &auxCh, size_t counter );
	static int		read_value( char *pos );
	void			finish_channel_details( std::vector<Channel> &auxCh  );
	void			read_cts_file( std::vector<Channel> &auxCh );
};

#endif // MTU_H