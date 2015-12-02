// The following ifdef block is the standard way of creating macros which make exporting 
// from a DLL simpler. All files within this DLL are compiled with the ENGINE_EXPORTS
// symbol defined on the command line. This symbol should not be defined on any project
// that uses this DLL. This way any other project whose source files include this file see 
// ENGINE_API functions as being imported from a DLL, whereas this DLL sees symbols
// defined with this macro as being exported.
#ifdef ENGINE_EXPORTS
#define ENGINE_API __declspec(dllexport)
#else
#define ENGINE_API __declspec(dllimport)
#endif

#include "stdafx.h"

#define N_BANDS 4

#define N_INPUT N_BANDS*2
#define N_OUTPUT 2
#define N_FILES N_INPUT

extern "C" ENGINE_API int process(	int inputc,
									int outputc,
									char *input[],
									char *output[],
									unsigned int type,
									int decimation,
									BOOL robust);

extern "C" JNIEXPORT jint JNICALL Java_com_mt_Engine_process(	JNIEnv *env,
																jobject this_object,
																jobjectArray j_input,
																jobjectArray j_output,
																jint j_type,
																jint j_decimation,
																jboolean j_robust);

class CEngineBase {
public:
	CEngineBase(int inputc,
				int outputc,
				char *input[],
				char *output[],
				FILE_TYPE type,
				int decimation);

	virtual ~CEngineBase() {}
	
	int process();

protected:
	virtual int process_emi() = 0;
	virtual int process_lmt() = 0;
	virtual int process_sio() = 0;
	virtual int process_mtu() = 0;
	virtual void reset() = 0;

	int			m_inputc;
	char*		m_input[N_INPUT];
	int			m_outputc;
	char*		m_output[N_OUTPUT];
	FILE_TYPE	m_type;
	int			m_decimation;
	
	//Header m_files[N_FILES];
};

class CEngine : public CEngineBase {
public:
	CEngine(int inputc,
			int outputc,
			char *input[],
			char *output[],
			FILE_TYPE type,
			int decimation) : CEngineBase(	inputc,
											outputc,
											input,
											output,
											type,
											decimation) {}

	virtual ~CEngine() {}

private:
	int process_emi();
	int process_lmt();
	int process_sio();
	int process_mtu();
	void reset() {}

	int process_emi_previous();
	int process_lmt_previous();
	int process_sio_previous();
	
	//FileExtractor			m_extractor;
	//CrossPowerTransformer	m_crosspower;
	//FFTTransformer			m_fft;
	//FileLoader				m_loader;
};

class CEngineRobust : public CEngineBase {
public:
	CEngineRobust(	int inputc,
					int outputc,
					char *input[],
					char *output[],
					FILE_TYPE type,
					int decimation);

	virtual ~CEngineRobust();

private:
	int process_emi();
	int process_lmt();
	int process_sio();
	int process_mtu();
	void reset() {}

	int process_emi_previous();
	
	//FileExtractor			m_extractor;
	//CrossPowerTransformer	m_crosspower;
	//FFTTransformer			m_fft;
	//FileLoader				m_loader;
	//RobustTransformer		m_robustTransformer;

	Array<double> &m_freqall;
	Array<double> &m_parall;
	Array<double> &m_parerrall;
	Array<double> &m_dis;
	Array<double> &m_c;
};

