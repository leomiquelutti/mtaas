#ifdef ENGINE_EXPORTS
#define ENGINE_API __declspec(dllexport)
#else
#define ENGINE_API __declspec(dllimport)
#endif

#include "stdafx.h"
#include "Engine.h"

#define N_BANDS 4

#define N_INPUT N_BANDS*2
#define N_OUTPUT 2
#define N_FILES N_INPUT

int main()
{
	// INPUTS
	TS_TO_FFT ts_to_fft = FIXED_WINDOW_LENGTH;
	ESTIMATOR_TYPE estimator_type = LEAST_SQUARES;

	// get stations name, path and type
	std::vector<StationBase> station = DirectoryProperties::initialize_stations( std::string("G:\\teste_mtprox") );

	// read time-series
	for( int i = 0; i < station.size(); i++ )
		station[i].read_time_series();
	
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