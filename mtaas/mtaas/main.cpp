#ifdef ENGINE_EXPORTS
#define ENGINE_API __declspec(dllexport)
#else
#define ENGINE_API __declspec(dllimport)
#endif

#include "stdafx.h"
#include "Engine.h"
#include "DirectoryProperties.h"
#include "Mtu.h"

#define N_BANDS 4

#define N_INPUT N_BANDS*2
#define N_OUTPUT 2
#define N_FILES N_INPUT

int main()
{
	// INPUTS
	TS_TO_FFT_TYPE ts2fft_type = FIXED_WINDOW_LENGTH;
	ESTIMATOR_TYPE estimator_type = LEAST_SQUARES;
	RR_OR_SS_TYPE rrorss_type = SINGLE_SITE;
	std::string inputPath = "C:\\Users\\Usuario\\Google Drive\\Documentos\\MT data\\mtu ss";
	//std::string inputPath = "C:\\Users\\leonardo\\Google Drive\\Documentos\\MT data\\mtu";

	// object to explore the main path and its subfolders
	DirectoryProperties dirInfo( inputPath );

	// get stations name, path and type
	std::vector<StationBase> station = dirInfo.initialize_stations();

	// read raw time-series
	// TODO - define timeVector properly for RR
	for( int i = 0; i < station.size(); i++ )
		station[i].read_time_series( &dirInfo );

	// build the draft for the time vector, considering all continuous acquisition
	// TODO - do it correctly
	Utils::draft_build_time_vector( station );

	//// checkouts
	//for( int istn = 0; istn < station.size(); istn++ ) {
	//	for( int its = 0; its < station[istn].ts.size(); its++ ) {		
	//		//station[istn].ts[its].timeVector[0].print();
	//		//std::cout << endl << istn << " " << its << endl;
	//		for( int irow = 0; irow < station[istn].ts[its].ch[0].systemResponseFreqs[0].getDim(0); irow++ ) {
	//			//std::cout << station[istn].ts[its].ch[0].systemResponseFreqs[0](irow) << " ";
	//			for( int ich = 0; ich < station[istn].ts[its].ch.size(); ich++ )	{
	//				station[istn].ts[its].ch[ich].gain = 1;
	//				std::cout << station[istn].ts[its].ch[ich].gain << " ";
	//			}
	//			std::cout << std::endl;
	//		}
	//	}
	//			std::cout << std::endl;
	//}

	// extract corrected (to physical units) Fourier coefficients 
	StationBase::get_all_FCs( station, ts2fft_type, rrorss_type );

	// evaluate Z 
	StationBase::get_Z( station, estimator_type );

	//// checkouts
	//for( int istn = 0; istn < station.size(); istn++ ) {
	//	for( int its = 0; its < station[istn].ts.size(); its++ ) {		
	//		cout << station[istn].ts[its].combination[0].idxStn[0] << endl;	
	//		cout << station[istn].ts[its].combination[0].idxTs[0] << endl;		
	//		station[istn].ts[its].combination[0].idxBgn[0].print();
	//		station[istn].ts[its].combination[0].idxEnd[0].print();
	//	}
	//	std::cout << std::endl;
	//}

	//			// checkout
	//for( int istn = 0; istn < station.size(); istn++ )
	//	for( int its = 0; its < station[istn].ts.size(); its++ )
	//		for( int icomb = 0; icomb < station[istn].ts[its].combination.size(); icomb++ ) {
	//			for( int ifreq = 0; ifreq < station[istn].ts[its].combination[icomb].fr.size(); ifreq++ ) 
	//				cout << istn << " " << its << " " << icomb << " " << ifreq << " " << station[istn].ts[its].combination[icomb].fr[ifreq].frequency << endl;
	//			cout << endl;
	//		}

	WriteOutputs::WriteOutputs( station );



	// delete everything allocated
	// TODO
	Utils::delete_all( &station );



	
	//int inputc;
	//int outputc;
	//char *input;
	//char *output;
	//unsigned int type = FILE_TYPE_MTU;
	//int decimation = 5;
	//BOOL robust = true;

	//CEngineBase *engine = NULL;
	//if(robust)
	//	engine = new CEngineRobust(	inputc,
	//								outputc,
	//								&input,
	//								&output,
	//								static_cast<FILE_TYPE>(type),
	//								decimation);
	//else
	//	engine = new CEngine(	inputc,
	//							outputc,
	//							&input,
	//							&output,
	//							static_cast<FILE_TYPE>(type),
	//							decimation);

	//int ret = engine->process();
	//delete engine;

	//return ret;
}