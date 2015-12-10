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
	TS_TO_FFT ts_to_fft = FIXED_WINDOW_LENGTH;
	ESTIMATOR_TYPE estimator_type = LEAST_SQUARES;
	std::string inputPath = "C:\\Users\\Usuario\\Google Drive\\Documentos\\MT data\\mtu";

	// object to explore the main path and its subfolders
	DirectoryProperties dirInfo( inputPath );

	// get stations name, path and type
	std::vector<StationBase> station = dirInfo.initialize_stations();

	// read time-series
	for( int i = 0; i < station.size(); i++ )
		station[i].read_time_series( &dirInfo );

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