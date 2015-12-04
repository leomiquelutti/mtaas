#include "stdafx.h"
#include "Engine.h"

CEngineRobust::CEngineRobust(	int inputc,	int outputc,
								char *input[], char *output[],
								FILE_TYPE type, int decimation) :
CEngineBase(inputc, outputc, input, output, type, decimation),
m_freqall((*new Array<double>(ArrayUtil::GetArrayDescriptor(1, 256), -9999))),
m_parall((*new Array<double>(ArrayUtil::GetArrayDescriptor(2, 256, 64), -9999))),
m_parerrall((*new Array<double>(ArrayUtil::GetArrayDescriptor(2, 256, 64), -9999))),
m_dis((*new Array<double>(ArrayUtil::GetArrayDescriptor(3, 16, 1024, 16)))),
m_c((*new Array<double>(ArrayUtil::GetArrayDescriptor(4, 16, 16, 16, 4)))) {
}

CEngineRobust::~CEngineRobust() {
	delete &m_freqall;
	delete &m_parall;
	delete &m_parerrall;
	delete &m_dis;
	delete &m_c;
}

int CEngineRobust::process_emi()
{
	std::cout << "CEngineRobust::process_emi" << std::endl;

	return 0;
}

int CEngineRobust::process_lmt()
{
	std::cout << "CEngineRobust::process_lmt" << std::endl;

	return 0;
}

int CEngineRobust::process_sio()
{
	std::cout << "CEngineRobust::process_sio" << std::endl;

	return 0;
}

int CEngineRobust::process_mtu()
{
	std::cout << "CEngineRobust::process_mtu" << std::endl;

	std::string inputTbl1 = "G:\\ANP\\parana\\dados_recebidos_141112\\PHASE 5 POINT BY POINT\\02-139\\RAW-BIN\\1327A21B.TBL";
	std::string inputTSn1 = "G:\\ANP\\parana\\dados_recebidos_141112\\PHASE 5 POINT BY POINT\\02-139\\RAW-BIN\\1327A21B.TS5";
	std::string inputCts = "G:\\ANP\\parana\\dados_recebidos_141112\\PHASE 5 POINT BY POINT\\02-139\\RAW-TXT\\1327A21B.CTS";

	std::string inputTbl2 = "G:\\ANP\\parana\\dados_recebidos_141112\\PHASE 5 POINT BY POINT\\02-138\\RAW-BIN\\1327924B.TBL";
	std::string inputTSn2 = "G:\\ANP\\parana\\dados_recebidos_141112\\PHASE 5 POINT BY POINT\\02-138\\RAW-BIN\\1327924B.TS3";
	//Array<double> *data;
	//Array<double> *timeVector;

	// Array<double> taper = dpss<double>( 134, 4 );

	//this->tbl = station[0].read_tbl( inputTbl1 );
	//station[0].tbl = station[0].read_tbl( inputTbl1 );
	//Array<double> data1 = station[0].read_mtu_data( &station[0], inputTSn1 );
	//Array<double> timeVector1 = station[0].get_mtu_time_vector( &station[0], inputTSn1 );
	//Array<Complex> correctionToData1 = station[0].read_cts_file( &station[0], inputCts );
	//
	//station[0].MtuTsBand = station[0].get_phoenix_TS_band( inputTSn1 );
	//station[0].isContinuous = station[0].is_acquisition_continuos( &timeVector1 );

	//data1 = data1.detrend()*taper;

	//station[1].tbl = station[1].read_tbl( inputTbl2 );
	//Array<Complex> data2 = station[1].read_mtu_data( &station[1], inputTSn2 );
	//Array<double> timeVector2 = station[1].get_mtu_time_vector( &station[1], inputTSn2 );

	//delete[] station;

	return 0;
}

