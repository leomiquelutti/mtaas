//// Engine.cpp : Defines the exported functions for the DLL application.

#include "stdafx.h"
#include "Engine.h"

//ENGINE_API int process(	int inputc,
//						int outputc,
//						char *input[],
//						char *output[],
//						unsigned int type,
//						int decimation,
//						BOOL robust)
//{
//	if(inputc > N_INPUT || outputc > N_OUTPUT || type > FILE_TYPE_ADU_07)
//		return -1;
//
//	CEngineBase *engine = NULL;
//	if(robust)
//		engine = new CEngineRobust(	inputc,
//									outputc,
//									input,
//									output,
//									static_cast<FILE_TYPE>(type),
//									decimation);
//	else
//		engine = new CEngine(	inputc,
//								outputc,
//								input,
//								output,
//								static_cast<FILE_TYPE>(type),
//								decimation);
//
//	int ret = engine->process();
//	delete engine;
//
//	return ret;
//}
//
//JNIEXPORT jint JNICALL Java_com_mt_Engine_process(	JNIEnv *env,
//													jobject this_object,
//													jobjectArray j_input,
//													jobjectArray j_output,
//													jint j_type,
//													jint j_decimation,
//													jboolean j_robust)
//{
//	int inputc = env->GetArrayLength(j_input);
//	int outputc = env->GetArrayLength(j_output);
//
//	if(inputc > N_INPUT || outputc > N_OUTPUT || j_type > FILE_TYPE_SIO)
//		return -1;
//	
//	char *input[N_INPUT];
//	for (int i = 0; i < inputc; i++)
//	{
//		jstring tmp = static_cast<jstring>(env->GetObjectArrayElement(j_input, i));
//		input[i] = const_cast<char*>(env->GetStringUTFChars(tmp, 0));
//	}
//	
//	char *output[N_OUTPUT];
//	for (int i = 0; i < outputc; i++)
//	{
//		jstring tmp = static_cast<jstring>(env->GetObjectArrayElement(j_output, i));
//		output[i] = const_cast<char*>(env->GetStringUTFChars(tmp, 0));
//	}
//
//	CEngineBase *engine = NULL;
//	if(static_cast<BOOL>(j_robust))
//		engine = new CEngineRobust(	inputc,
//									outputc,
//									input,
//									output,
//									static_cast<FILE_TYPE>(j_type),
//									static_cast<int>(j_decimation));
//	else
//		engine = new CEngine(	inputc,
//								outputc,
//								input,
//								output,
//								static_cast<FILE_TYPE>(j_type),
//								static_cast<int>(j_decimation));
//
//	int ret = engine->process();
//	delete engine;
//	
//	return ret;
//}
//

CEngineBase::CEngineBase(	int inputc,
							int outputc,
							char *input[],
							char *output[],
							FILE_TYPE type,
							int decimation) : m_inputc(inputc), m_outputc(outputc), m_type(type), m_decimation(decimation)
{
	memset(m_input, NULL, sizeof(m_input));
	for(int i = 0; i < inputc; i++)
		m_input[i] = input[i];

	memset(m_output, NULL, sizeof(m_output));
	for(int i = 0; i < outputc; i++)
		m_output[i] = output[i];

	std::cout << "Type:" << m_type << std::endl;
	std::cout << "Decimation:" << m_decimation << std::endl;

	std::cout << "Input files:" << std::endl;
	for(int i = 0; i < m_inputc; i++)
		std::cout << m_input[i] << std::endl;

	std::cout << "Output files:" << std::endl;
	for(int i = 0; i < m_outputc; i++)
		std::cout << m_output[i] << std::endl;

}

int CEngine::process_emi()
{
	std::cout << "CEngine::process_emi" << std::endl;

	return 0;
}

int CEngine::process_lmt()
{
	std::cout << "CEngine::process_lmt" << std::endl;

	return 0;
}

int CEngine::process_sio()
{
	std::cout << "CEngine::process_sio" << std::endl;

	return 0;
}

int CEngine::process_mtu()
{
	std::cout << "CEngine::process_mtu" << std::endl;

	return 0;
}

int CEngineBase::process()
{
	switch(m_type) 
	{
	case FILE_TYPE_EMI:
		return process_emi();
	case FILE_TYPE_LMT:
		return process_lmt();
	case FILE_TYPE_SIO:
		return process_sio();
	case FILE_TYPE_MTU:
		return process_mtu();
	default:
		return -1;
	}
}
