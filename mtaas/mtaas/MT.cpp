#include "MT.h"
#include <sstream>

StationBase& StationBase::operator = (const StationBase& element)
{
	*this = element;
	return *this;
}

void StationBase::initialize_parameters()
{
	maxFreqInHertz = 1;
}

std::string Date::getDateStr()
{
	return std::to_string(startYear) + "-"
		+ std::to_string(startMonth) + "-"
		+ std::to_string(startDay) + " "
		+ std::to_string(startHour) + ":"
		+ std::to_string(startMinute) + ":"
		+ std::to_string(startSecond);
}

bool StationBase::is_acquisition_continuos( Array<double> *timeVector )
{
	// check consistency
	for( int i = 1; i < timeVector->GetDescriptor().GetDim( 0 ); i++ ) {
		if( (*timeVector)(i,0) - (*timeVector)(i-1,0) < 1 )
		{
			std::cout << "There's a problem in time vector" <<std::endl;
			std::cout << "timeVector(" << i << ") - timeVector(" << i-1 << ") < 1" << std::endl;
			exit(-42);
		}
	}

	// check continuity
	for( int i = 1; i < timeVector->GetDescriptor().GetDim( 0 ); i++ ) {
		if( (*timeVector)(i,0) - (*timeVector)(i-1,0) > 1 )
			return false;
	}
	return true;
}

//Header::Header() :
//	freq((*new Array<double>(ArrayUtil::GetArrayDescriptor(1, 128)))),
//	freqall((*new Array<double>(ArrayUtil::GetArrayDescriptor(1, 128)))),
//	vpsmpl((*new Array<double>(ArrayUtil::GetArrayDescriptor(1, 128)))),
//	par((*new Array<double>(ArrayUtil::GetArrayDescriptor(2, 128, 128)))),
//	parerr((*new Array<double>(ArrayUtil::GetArrayDescriptor(2, 128, 128)))),
//	parall((*new Array<double>(ArrayUtil::GetArrayDescriptor(2, 128, 128)))),
//	parerrall((*new Array<double>(ArrayUtil::GetArrayDescriptor(2, 128, 128)))),
//	chmt((*new Array<int>(ArrayUtil::GetArrayDescriptor(1, 8)))),
//	eh((*new Array<double>(ArrayUtil::GetArrayDescriptor(1, 8)))),
//	pos((*new Array<int>(ArrayUtil::GetArrayDescriptor(1, 8)))),
//	pltorder((*new Array<double>(ArrayUtil::GetArrayDescriptor(1, 16)))),
//	dateform((*new Array<int>(ArrayUtil::GetArrayDescriptor(1, 8)))),
//	tlenght((*new Array<double>(ArrayUtil::GetArrayDescriptor(1, 8)))),
//	indy((*new Array<double>(ArrayUtil::GetArrayDescriptor(1, 32)))),
//	efpar((*new Array<double>(ArrayUtil::GetArrayDescriptor(2, 8, 8)))),
//	gain((*new Array<int>(ArrayUtil::GetArrayDescriptor(1, 16)))),
//	filech((*new Array<int>(ArrayUtil::GetArrayDescriptor(1, 8)))),
//	datach((*new Array<int>(ArrayUtil::GetArrayDescriptor(1, 8)))),
//	sn((*new Array<double>(ArrayUtil::GetArrayDescriptor(1, 8)))),
//	cavg((*new Array<double>(ArrayUtil::GetArrayDescriptor(4, 16, 16, 16, 4)))),
//	navg((*new Array<double>(ArrayUtil::GetArrayDescriptor(1, 16)))),
//	fqs((*new Array<double>(ArrayUtil::GetArrayDescriptor(1, 16)))),
//	bws((*new Array<double>(ArrayUtil::GetArrayDescriptor(1, 16)))),
//	nfcs((*new Array<double>(ArrayUtil::GetArrayDescriptor(1, 16)))),
//	navs((*new Array<double>(ArrayUtil::GetArrayDescriptor(1, 16)))),
//	w((*new Array<double>(ArrayUtil::GetArrayDescriptor(1, 512)))),
//	chsat((*new Array<double>(ArrayUtil::GetArrayDescriptor(1, 32)))),
//	time((*new Array<double>(ArrayUtil::GetArrayDescriptor(1, 1000000)))),
//	fqavg((*new Array<double>(ArrayUtil::GetArrayDescriptor(1, 16)))),
//	bwavg((*new Array<double>(ArrayUtil::GetArrayDescriptor(1, 16)))),
//	ifq((*new Array<int>(ArrayUtil::GetArrayDescriptor(1, 32)))),
//	nfc((*new Array<double>(ArrayUtil::GetArrayDescriptor(1, 128)))),
//	bw((*new Array<double>(ArrayUtil::GetArrayDescriptor(1, 16)))),
//	chnuse((*new Array<double>(ArrayUtil::GetArrayDescriptor(1, 8))))
//{
//}
//
//Header::~Header()
//{
//	delete &freq;
//	delete &freqall;
//	delete &vpsmpl;
//	delete &par;
//	delete &parerr;
//	delete &parall;
//	delete &parerrall;
//	delete &chmt;
//	delete &eh;
//	delete &pos;
//	delete &pltorder;
//	delete &dateform;
//	delete &tlenght;
//	delete &indy;
//	delete &efpar;
//	delete &gain;
//	delete &filech;
//	delete &datach;
//	delete &sn;
//	delete &cavg;
//	delete &navg;
//	delete &fqs;
//	delete &bws;
//	delete &nfcs;
//	delete &navs;
//	delete &w;
//	delete &chsat;
//	delete &time;
//	delete &fqavg;
//	delete &bwavg;
//	delete &ifq;
//	delete &nfc;
//	delete &bw;
//	delete &chnuse;
//}
//
//void Util::initstr(Header *files, int cont)
//{
//	files[cont].stdlim = 4;
//	files[cont].say = 1;   
//	files[cont].q = 3;
//	files[cont].fhi = 0;
//	files[cont].fmax = 60;
//	files[cont].cmax = 15;
//	files[cont].satmx = 16380;
//
//	for(int i = 0; i < 8; i++)
//		files[cont].chmt(i) = 0;
//
//	files[cont].eh(1) = 12;
//	files[cont].eh(2) = 13; 
//	files[cont].eh(3) = 14; 
//	files[cont].eh(4) = 15; 
//	files[cont].eh(5) = 7; 
//	files[cont].eh(6) = 8; 
//	files[cont].eh(7) = 9;
//
//	files[cont].chmt(1) = 1;
//	files[cont].chmt(2) = 2;
//	files[cont].chmt(3) = 3;
//	files[cont].chmt(4) = 4;
//	files[cont].chmt(5) = 5;
//	files[cont].chmt(6) = 3;
//	files[cont].chmt(7) = 4;
//	
//	files[cont].ropt = 0;
//	files[cont].strk = 0;
//	files[cont].deci = 1;
//
//	files[cont].nfiles = 1;
//
//	for(int i = 0; i < 8; i++)
//		files[cont].vpsmpl(i) = 0;
//
//	for(int i = 1; i <= 75; i++)
//	{
//		files[cont].freq(i) = -9999;
//		for (int j = 1; j <= 59; j++)
//		{
//			files[cont].par(i, j) = -9999;
//			files[cont].parerr(i, j) = -9999;
//		}
//	}
//
//	files[cont].ncpblk = 1;
//}
//
//void Util::zerarparametros(Header *files, int cont)
//{
//	files[cont].nfq = 1;
//	files[cont].q = 3;
//	files[cont].fhi = 0;
//	files[cont].fmax = 60;
//	files[cont].lngth = 0;
//
//	for(int i = 0; i < 17; i++)
//		files[cont].ifq(i) = 0;
//
//	for(int i = 0; i < 512; i++)
//		files[cont].time(i) = 0;
//
//	for(int i = 0; i < 81; i++)
//		files[cont].nfc(i) = 0;
//
//	for(int i = 0; i < 11; i++)
//	{
//		files[cont].bw(i) = 0;
//		files[cont].bwavg(i) = 0;
//		files[cont].fqavg(i) = 0;
//	}
//
//	for(int i = 1; i <= 75; i++)
//	{
//		files[cont].freq(i) = -9999;
//		for(int j = 1; j <= 59; j++)
//		{
//			files[cont].par(i, j) = -9999;
//			files[cont].parerr(i, j) = -9999;
//		}
//	}
//}
//
//int Util::finalresults(Header *files, int cont, Array<double> &freqall, Array<double> &parall, Array<double> &parerrall)
//{
//	int indice1 = 0, indice2 = 0;
//
//	//Achar índice1 até o qual files[cont].freq, files[cont].par e files[cont].parerr foi preenchida
//	for (int i = 1; i <= 200; i++)
//	{
//		if(files[cont].freq(i) == -9999)
//		{
//			indice1 = i;
//			break;
//		}
//	}
//
//	//Achar índice2 até o qual freqall, parall e parerrall foi preenchida
//	for (int i = 1; i <= 200; i++)
//	{
//		if(freqall(i) == -9999)
//		{
//			indice2 = i;
//			break;
//		}
//	}
//
//	//Preencher freqall, parall e parerrall
//	for (int i = 1; i < indice1; i++)
//	{
//		freqall(i + indice2 - 1) = files[cont].freq(i);
//		for (int j = 1; j <= 59; j++)
//		{
//			parall(i + indice2 - 1, j) = files[cont].par(i, j);
//			parerrall(i + indice2 - 1, j) = files[cont].parerr(i, j);
//
//			if(files[cont].par(i, j) == -9999)
//				parall(i + indice2 - 1, j) = 0;
//
//			if(files[cont].parerr(i, j) == -9999)
//				parerrall(i + indice2 - 1, j) = 0;
//		}
//	}
//
//	//Achar índice2 até o qual freqall, parall e parerrall foi preenchida
//	for (int i = 1; i <= 200; i++)
//	{
//		if(freqall(i) == -9999)
//		{
//			indice2 = i;
//			break;
//		}
//	}
//
//	return indice2;
//}
//
//void Util::copyfreqselect(Header *files, int cont)
//{
//	files[4].nfq = files[cont].nfq;
//	files[4].fhi = files[cont].fhi;
//	files[4].f1 = files[cont].f1;
//	files[4].lwind = files[cont].lwind;
//
//	for (int i = 0; i < 16; i++)
//	{
//		for (int j = 0; j < 16; j++)
//		{
//			for (int k = 0; k < 11; k++)
//			{
//				for (int l = 0; l < 3; l++)
//				{
//					files[4].cavg(i, j, k, l) = files[cont].cavg(i, j, k, l);
//				}
//			}
//		}
//	}
//
//	for (int k = 0; k < 11; k++)
//	{
//		files[4].bw(k) = files[cont].bw(k);
//		files[4].fqavg(k) = files[cont].fqavg(k);
//		files[4].bwavg(k) = files[cont].bwavg(k);
//		files[4].navg(k)  = files[cont].navg(k);
//	}
//
//	for (int j = 0; j < 16; j++)
//	{
//		files[4].fqs(j) = files[cont].fqs(j);
//		files[4].bws(j) = files[cont].bws(j);
//		files[4].nfcs(j) = files[cont].nfcs(j);
//		files[4].navs(j) = files[cont].navs(j);
//	}
//
//	for (int j = 0; j < 17; j++)
//		files[4].ifq(j) = files[cont].ifq(j);
//
//	for (int j = 0; j < 81; j++)
//		files[4].nfc(j) = files[cont].nfc(j);
//
//	for (int j = 0; j < 101; j++)
//		files[4].freq(j) = files[cont].freq(j);
//
//	for (int j = 0; j < 500; j++)
//		files[4].w(j) = files[cont].w(j);
//
//	return;
//}
//	
//double Util::minfv(Array<double> &vector)
//{
//	double menor;
//	int dimv = sizeof(vector) / sizeof(double);
//
//	if (vector(2) <= vector(1)) 
//		menor = vector(2);
//	else 
//		menor = vector(1);
//
//	for (int i=3; i<=dimv; i++) {
//		if (vector(i) <= menor)
//			menor = vector(i);
//	}
//
//	return menor;
//}
//
//double Util::maxfv(Array<double> &vector)
//{
//	double maior;
//	int dimv = sizeof(vector) / sizeof(double);
//
//	if (vector(2) > vector(1)) 
//		maior = vector(2);
//	else 
//		maior = vector(1);
//
//	for (int i=3; i<=dimv; i++) {
//		if (vector(i) >= maior)
//			maior = vector(i);
//	}
//
//	return maior;
//}
//
//double Util::mindbv(Array<double> &vector, int dim)
//{
//	double menor;
//
//	if (vector(2) <= vector(1)) {menor = vector(2);}
//	else {menor = vector(1);}
//
//	for (int i = 3; i <= dim; i++) {
//		if (vector(i) <= menor)
//			menor = vector(i);
//	}
//
//	return (menor);
//}
//
//double Util::maxdbv(Array<double> &vector, int dim)
//{
//
//			/*for( int i = 1; i <= 45; i++ )
//				printf("%e\n",vector( i ));
//			printf("\n");*/
//
//	double maior;
//
//	if (vector(2) > vector(1)) {maior = vector(2);}
//	else {maior = vector(1);}
//
//	for (int i = 3; i <= dim; i++) {
//		if (vector(i) >= maior)
//			maior = vector(i);
//	}
//
//	return (maior);
//}
//
//double Util::mindb (double n1, double n2) {
//Devolve o menor de dois números do tipo double
//	=============================================
//	n1: primeiro número
//	n2: segundo número */
//
//	double menor;
//
//	if (n1 == n2) 
//		menor = n1;
//	else if (n1 > n2) 
//		menor = n2;
//	else 
//		menor = n1;
//
//	return (menor);
//}
//
//int Util::char2int32(char *cbuffer)
//{
//	int int32val = (int)(cbuffer[0]*pow((double)2,24) + cbuffer[1]*pow((double)2,16) + 
//					cbuffer[2]*pow((double)2,8) + cbuffer[3]*pow((double)2,0));
//
//	return int32val;
//}
//
//int Util::char2int16(char *cbuffer)
//{
//	int int16val = (int)(cbuffer[0]*pow((double)2,8) + cbuffer[1]*pow((double)2,0));
//
//	return int16val;
//}
//
//void Util::impedance(Array<double> &dis)
//{
//	FILE *f = NULL;
//
//	// define qual frequencia
//	int k = 3;
//
//	f = fopen("tensor de impedancia 3.txt","w");
//	if(!f)
//	{
//		printf("Erro ao abrir 'tensor de impedancia.txt'\n");
//		return;
//	}
//
//	fprintf(f,"frequencia numero %d.\n",k);
//	fprintf(f,"zxxr\tzxxi\tzxyr\tzxyi\tzyxr\tzyxi\tzyyr\tzyyi\n");
//
//	for(int i = 1; i <= 243; i++)
//	{
//		fprintf(f,"%lf\t",dis(k, i, 1));
//		fprintf(f,"%lf\t",dis(k, i, 2));
//		fprintf(f,"%lf\t",dis(k, i, 3));
//		fprintf(f,"%lf\t",dis(k, i, 4));
//		fprintf(f,"%lf\t",dis(k, i, 5));
//		fprintf(f,"%lf\t",dis(k, i, 6));
//		fprintf(f,"%lf\t",dis(k, i, 7));
//		fprintf(f,"%lf\n",dis(k, i, 8));
//	}
//
//	fprintf(f,"\0");
//	fclose(f);
//
//	// define qual frequencia
//	k = 5;
//
//	f = fopen("tensor de impedancia 5.txt","w");
//	if(!f)
//	{
//		printf("Erro ao abrir 'tensor de impedancia.txt'\n");
//		return;
//	}
//
//	fprintf(f,"frequencia numero %d.\n",k);
//	fprintf(f,"zxxr\tzxxi\tzxyr\tzxyi\tzyxr\tzyxi\tzyyr\tzyyi\n");
//
//	for(int i = 1; i <= 243; i++)
//	{
//		fprintf(f,"%lf\t",dis(k, i, 1));
//		fprintf(f,"%lf\t",dis(k, i, 2));
//		fprintf(f,"%lf\t",dis(k, i, 3));
//		fprintf(f,"%lf\t",dis(k, i, 4));
//		fprintf(f,"%lf\t",dis(k, i, 5));
//		fprintf(f,"%lf\t",dis(k, i, 6));
//		fprintf(f,"%lf\t",dis(k, i, 7));
//		fprintf(f,"%lf\n",dis(k, i, 8));
//	}
//
//	fprintf(f,"\0");
//	fclose(f);
//
//	// define qual frequencia
//	k = 7;
//
//	f = fopen("tensor de impedancia 7.txt","w");
//	if(!f)
//	{
//		printf("Erro ao abrir 'tensor de impedancia.txt'\n");
//		return;
//	}
//
//	fprintf(f,"frequencia numero %d.\n",k);
//	fprintf(f,"zxxr\tzxxi\tzxyr\tzxyi\tzyxr\tzyxi\tzyyr\tzyyi\n");
//
//	for(int i = 1; i <= 243; i++)
//	{
//		fprintf(f,"%lf\t",dis(k, i, 1));
//		fprintf(f,"%lf\t",dis(k, i, 2));
//		fprintf(f,"%lf\t",dis(k, i, 3));
//		fprintf(f,"%lf\t",dis(k, i, 4));
//		fprintf(f,"%lf\t",dis(k, i, 5));
//		fprintf(f,"%lf\t",dis(k, i, 6));
//		fprintf(f,"%lf\t",dis(k, i, 7));
//		fprintf(f,"%lf\n",dis(k, i, 8));
//	}
//
//	fprintf(f,"\0");
//	fclose(f);
//
//	return;
//}
//
//void Util::copy_ampf_and_phsf(Header *files, int cont, Array<double> &ampf, Array<double> &phsf, Array<double> &ampf_new, Array<double> &phsf_new)
//{
//	int index = files[ cont ].nochnl + files[ 4 ].nochnl;
//	for( int i = 0; i < files[ cont ].nfq; i++ ) {
//		for( int j = 0; j < index; j++ ) {
//			ampf_new( i*index + j ) = ampf( j + 1 , i + 1 );
//			phsf_new( i*index + j ) = phsf( j + 1 , i + 1 );			
//		}
//	}
//}
//
//Complex Util::conjugate( Complex a ) {
//	return std::conj( a );
//}
//
//void Util::qr_decomposition_via_GramSchimidt( Array<Complex> &A, Array<Complex> &Q, Array<Complex> &R ) {
//
//	int m = A.GetDescriptor().GetDim( 0 ); // amount of rows
//	int n = A.GetDescriptor().GetDim( 1 ); // amount of columns
//
//	double s;
//	Complex t;
//	for( int k = 0; k < n; k++ ) {
//		s = 0;
//		for( int j = 0; j < m; j++ )
//			s += pow(A( j, k ).real(),2) + pow(A( j, k ).imag(),2);
//			//s += A( j, k ).magnitude2();
//		R( k, k ) = Complex( sqrt( s ), 0 );
//		for( int j = 0; j < m; j++ )
//			Q( j, k ) = A( j, k )/R( k, k );
//		for( int i = k + 1; i < n; i++ ) {
//			t = Complex( 0, 0 );
//			for( int j = 1; j < m; j++ )
//				t += A( j, i )*Q( j, k );
//			R( k, i ) = t;
//			for( int j = 0; j < m; j++ )
//				A( j, i ) -= R( k, i )*Q( j, k );
//		}
//	}
//
//	return;
//}
//
//void Util::invert_2x2_complex_matrix( Array<Complex> &matrix ) {
//
//	Array<Complex> aux1(ArrayUtil::GetArrayDescriptor( 2, 2, 2 ) );
//
//	Complex aux2 = Complex( 1, 0 )/( matrix( 0 , 0 )*matrix( 1 , 1 ) - matrix( 0 , 1 )*matrix( 1 , 0 ) );
//	
//	for( int i = 0; i < 2; i++ ) {
//		for( int j = 0; j < 2; j++ )
//			aux1( i, j ) = matrix( i, j );
//	}
//
//	matrix( 0, 0 ) = aux2*aux1( 1 , 1 );
//	matrix( 0, 1 ) = -aux2*aux1( 0 , 1 );
//	matrix( 1, 0 ) = -aux2*aux1( 1 , 0 );
//	matrix( 1, 1 ) = aux2*aux1( 0 , 0 );
//
//	return;
//}
//
//int Util::get_icp( Header *files, int cont ) {
//	return files[cont].tppb*files[cont].tnblk;
//}
//
//int Util::get_number_of_events( Header *files, int cont, int frequencyNumber ) {
//	return files[ cont ].ifq( frequencyNumber + 1 ) - files[ cont ].ifq( frequencyNumber );
//}
//
//int Util::getInicialPositionOfEvent( Header *files, int cont, int frequencyNumber ) {
//	return files[cont].ifq(frequencyNumber) - files[cont].ifq(1);
//}
//
//void Util::print_complex_array( Array<Complex> &A ) {
//
//	switch( A.GetDescriptor().GetNDim() ) {
//	
//	case 1:
//		for( int i = 0; i < A.GetDescriptor().GetDim( 0 ); i++ )
//			printf("%e i%e\t", A( i ).real(), A( i ).imag() );
//		printf("\n");
//	
//	case 2:
//		for( int i = 0; i < A.GetDescriptor().GetDim( 0 ); i++ ) {
//			for( int j = 0; j < A.GetDescriptor().GetDim( 1 ); j++ )
//				printf("%e i%e\t", A( i, j ).real(), A( i, j ).imag() );
//		printf("\n");
//		}
//	
//	case 3:
//		for( int i = 0; i < A.GetDescriptor().GetDim( 0 ); i++ ) {
//			for( int j = 0; j < A.GetDescriptor().GetDim( 1 ); j++ ) {
//				for( int k = 0; k < A.GetDescriptor().GetDim( 2 ); k++ )
//					printf("%e i%e\t", A( i, j, k ).real(), A( i, j, k ).imag() );
//				printf("\n");
//			}
//			printf("\n");
//		}
//	}	
//	printf("\n");	
//
//	return;
//}
//
//void Util::print_double_array( Array<double> &A ) {
//
//	switch( A.GetDescriptor().GetNDim() ) {
//	
//	case 1:
//		for( int i = 0; i < A.GetDescriptor().GetDim( 0 ); i++ )
//			printf("%e\t", A( i ) );
//		printf("\n");
//	
//	case 2:
//		for( int i = 0; i < A.GetDescriptor().GetDim( 0 ); i++ ) {
//			for( int j = 0; j < A.GetDescriptor().GetDim( 1 ); j++ )
//				printf("%e\t", A( i, j ) );
//		printf("\n");
//		}
//	
//	case 3:
//		for( int i = 0; i < A.GetDescriptor().GetDim( 0 ); i++ ) {
//			for( int j = 0; j < A.GetDescriptor().GetDim( 1 ); j++ ) {
//				for( int k = 0; k < A.GetDescriptor().GetDim( 2 ); k++ )
//					printf("%e\t", A( i, j, k ) );
//				printf("\n");
//			}
//			printf("\n");
//		}
//	}	
//	printf("\n");	
//
//	return;
//
//	return;
//}