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

//int CEngineRobust::process_emi_previous()
//{
//	int decimation = 1, sc = 0, index = 0;
//	bool rr = false;
//	for (int cont = 0; cont < m_inputc; cont++)
//	{
//		Util::initstr(m_files, cont);
//
//		m_files[cont].filetype = m_type;
//		m_files[cont].deci = m_decimation;
//		strcpy(m_files[cont].filename, m_input[cont]);
//
//		if(cont + 1 < m_inputc && strstr("rr_", m_input[cont + 1]))
//		{
//			rr = true;
//			strcpy(m_files[4].filename, m_input[cont + 1]);
//			m_files[cont].nfiles = 2;
//			m_files[cont].chmt(6) = m_files[cont].chmt(6) + 5;
//			m_files[cont].chmt(7) = m_files[cont].chmt(7) + 5;
//
//			for(int i = 0; i < 8; i++)
//				m_files[4].vpsmpl(i) = 0;
//		}
//
//		// Read header information
//		for(int i = 1; i <= m_files[cont].nfiles; i++)
//		{
//			if (i == 1)
//				m_extractor.readheader_emi(m_files, cont);
//			else 
//				m_extractor.readheader_emi(m_files, 4);
//
//				sc = 1000;
//		}
//
//		// !!!!!!!!!!!!!!!!!!!!!!!!!! OLHAR ESTA SEÇÃO DEPOIS DE PR IMPLEMENTADO
//		for (int i = 0; i < 16; i++)
//		{
//			for (int j = 0; j < 16; j++)
//			{
//				for (int k = 0; k < 16; k++)
//				{
//					for (int l = 0; l < 3; l++)
//					{
//						m_files[cont].cavg(i, j, k, l) = 0;
//						m_c(i, j, k, l) = m_files[cont].cavg(i, j, k, l);
//					}
//				}
//			}
//		}
//		// !!!!!!!!!!!!!!!!!!!!!!!!!!
//
//		// selection of frequencies and time windows to be used
//		// maximum number of frequencies allowed = 40
//		m_fft.freqselect(m_files, cont);
//		if (rr)
//			Util::copyfreqselect(m_files, cont);
//
//		// vector that contains time series data converted
//		Array<double> dat(ArrayUtil::GetArrayDescriptor(1, (m_files[cont].nochnl + 2*(rr) )*m_files[cont].lngth));	
//
//		// vector for checking whether there is response file
//		Array<int> rspfile(ArrayUtil::GetArrayDescriptor(1, 5));
//		for(int i = 0; i < 5; i++)
//			rspfile(i) = 1;
//
//		// output of function read_sys_cor_emi, used for correction of acquired measurements
//		Array<double> ampf(ArrayUtil::GetArrayDescriptor(2, 31, 21));
//		Array<double> phsf(ArrayUtil::GetArrayDescriptor(2, 31, 21));
//		for(int i = 0; i < 31; i++)
//		{
//			for(int j = 0; j < 21; j++)
//			{
//				ampf(i, j)=0;
//				phsf(i, j)=0;
//			}
//		}
//
//		for(int i = 1; i <= m_files[cont].nfiles; i++) {
//			if (i == 1)
//				m_extractor.read_sys_cor_emi(m_files, cont, ampf, phsf, rspfile);
//			else 
//				m_extractor.read_sys_cor_emi(m_files, 4, ampf, phsf, rspfile);			
//		}
//
//		// Stack number
//		int ncpb = 0;
//
//		// index of cross power blocks
//		int cpsay = 1; 
//		int ijk = 1, start = 0, startrr = 0;
//
//		start = m_extractor.localiza(m_files, cont);
//		if (rr)
//			startrr = m_extractor.localiza(m_files, 4);
//
//		for (int i = 0; i <= ( m_files[cont].nochnl + 2*(rr) )*m_files[cont].lngth; i++)
//			dat(i) = 0;	
//
//		// $$$$$$$$$ - vector declarations
//
//		// vector that contains complex crosspower data values
//		Array<Complex> crossPowers(ArrayUtil::GetArrayDescriptor(4, m_files[cont].nochnl + 2*(rr), m_files[cont].nochnl + 2*(rr), m_files[cont].nfq, m_files[cont].tppb*m_files[cont].tnblk ) );
//		// vector that contains complex fft data values
//		Array<Complex> fftValues(ArrayUtil::GetArrayDescriptor(3, m_files[cont].nochnl, (m_files[ cont ].ifq(m_files[ cont ].nfq + 1) - m_files[ cont ].ifq(1)), m_files[cont].tppb*m_files[cont].tnblk ) );
//		// vector that contains real (between 0 and 1) multiple coherence data values
//		Array<double> multiCoh(ArrayUtil::GetArrayDescriptor(2, m_files[cont].tppb*m_files[cont].tnblk, 3 ) );
//
//		// $$$$$$$$$ - end of vector declarations
//
//		// Read time series
//		for (int icp = 1; icp <= m_files[cont].tppb*m_files[cont].tnblk; icp++)
//		{
//			for (int ijk = 1; ijk <= m_files[cont].nfiles; ijk++)
//			{
//				if (ijk == 1)
//					start = m_extractor.readts_emi(m_files, cont, dat, ijk, icp, start);
//				else
//					startrr = m_extractor.readts_emi(m_files, 4, dat, ijk, icp, startrr);
//			}
//
//			/*for( int i = 0; i < 5*files[cont].lngth; i++ ) {
//				printf("dat(%d) = %e\n", i, dat(i) );
//				if( 2*(i + 1)%files[cont].lngth == 0 )
//					printf("\n");
//			}*/
//
//			// time series saturation check
//			m_fft.sat_check(m_files, cont, dat);
//
//			// fast fourier transform (fft)
//			m_fft.fft9(m_files, cont, dat);
//
//			// Calcula auto e cross powers e average frequency range determinado por ifq e grava em "c"
//			m_crosspower.fillInFftValues(m_files, cont, dat, icp, fftValues);
//
//			m_robustTransformer.multipleCoherences(m_files, cont, fftValues, multiCoh, icp);
//
//			// Número de cross power blocks
//			//ncpb = ncpb + 1;
//
//			//// Quando ncpb for igual ao número de cross power blocks, entra nesse if - faz o stack e calcula os resultados.
//			//if (ncpb == files[cont].ncpblk)
//			//{
//			//	// Correção das constantes do equipamento e parâmetros de exploração
//			//	m_crosspower.cpcrrct(files, cont, 1, c);
//
//			//	// Correção da função resposta com valores obtidos em ampf e phsf
//			//	m_crosspower.sys_cor(files, cont, c, ampf, phsf);
//
//			//	// Cross power stack - média das cross powers
//			//	m_crosspower.cp_add(files, cont, c);
//
//			//	// Calcula MT tranfer function e seus derivados
//			//	double erot = m_loader.results(files, cont, icp, dis);
//
//			//	// Dizer quantos blocos foram calculados
//			//	cpsay = cpsay + 1;
//
//			//	// Zerar ncpb a fim de voltar nesse "if"
//			//	ncpb = 0;
//
//			//	for (int i = 0; i <= ( files[cont].nochnl + 2*(rr) )*files[cont].lngth; i++)
//			//		dat(i) = 0;
//
//			//	for (int i = 0; i < 16; i++)
//			//	{
//			//		for (int j = 0; j < 16; j++)
//			//		{
//			//			for (int k = 0; k < 16; k++)
//			//			{
//			//				for (int l = 0; l < 3; l++)
//			//				{
//			//					c(i, j, k, l) = 0;
//			//				}
//			//			}
//			//		}
//			//	}
//			//}
//		}
//
//		int icp = m_files[cont].tppb*m_files[cont].tnblk;
//		m_robustTransformer.robust_processing(m_files, cont, fftValues, multiCoh, ampf, phsf, icp);
//		
//		index = Util::finalresults(m_files, cont, m_freqall, m_parall, m_parerrall);
//
//		// Zera alguns parâmetros importantes caso deci != 1 e files[cont].par, files[cont].parerr e files[cont].freq
//		Util::zerarparametros(m_files, cont);
//
//		// Imprime os valores do tensor de impedancia
//		Util::impedance(m_dis);
//	}
//
//	// Imprime resultados
//	m_extractor.printresults(m_files, m_inputc - 1, m_freqall, m_parall, m_parerrall, index);
//
//	// Gera o arquivo edi
//	m_extractor.edi(m_files, m_inputc - 1, index - 1, m_freqall, m_parall, m_parerrall);	
//
//	return 0;
//}

int CEngineRobust::process_mtu()
{
	std::cout << "CEngineRobust::process_mtu" << std::endl;

	DirectoryProperties dirInfo;
	dirInfo.pathName = std::string( "G:\\teste_mtprox" );

	std::vector<StationBase> station = dirInfo.initialize_stations();

	for( int idxFolder = 0; idxFolder < station.size(); idxFolder++ )
	{
	}

	std::string inputTbl1 = "G:\\ANP\\parana\\dados_recebidos_141112\\PHASE 5 POINT BY POINT\\02-139\\RAW-BIN\\1327A21B.TBL";
	std::string inputTSn1 = "G:\\ANP\\parana\\dados_recebidos_141112\\PHASE 5 POINT BY POINT\\02-139\\RAW-BIN\\1327A21B.TS5";
	std::string inputCts = "G:\\ANP\\parana\\dados_recebidos_141112\\PHASE 5 POINT BY POINT\\02-139\\RAW-TXT\\1327A21B.CTS";

	std::string inputTbl2 = "G:\\ANP\\parana\\dados_recebidos_141112\\PHASE 5 POINT BY POINT\\02-138\\RAW-BIN\\1327924B.TBL";
	std::string inputTSn2 = "G:\\ANP\\parana\\dados_recebidos_141112\\PHASE 5 POINT BY POINT\\02-138\\RAW-BIN\\1327924B.TS3";
	//Array<double> *data;
	//Array<double> *timeVector;

	// Array<double> taper = dpss<double>( 134, 4 );
	
	station[0].tbl = station[0].read_tbl( inputTbl1 );
	Array<double> data1 = station[0].read_mtu_data( &station[0], inputTSn1 );
	Array<double> timeVector1 = station[0].get_mtu_time_vector( &station[0], inputTSn1 );
	Array<Complex> correctionToData1 = station[0].read_cts_file( &station[0], inputCts );
	
	station[0].MtuTsBand = station[0].get_phoenix_TS_band( inputTSn1 );
	station[0].isContinuous = station[0].is_acquisition_continuos( &timeVector1 );

	//data1 = data1.detrend()*taper;

	//station[1].tbl = station[1].read_tbl( inputTbl2 );
	//Array<Complex> data2 = station[1].read_mtu_data( &station[1], inputTSn2 );
	//Array<double> timeVector2 = station[1].get_mtu_time_vector( &station[1], inputTSn2 );

	//delete[] station;

	return 0;
}

