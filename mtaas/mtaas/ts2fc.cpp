#include "ts2fc.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cmath>

// constructor
Extract_FCs_FixedWindowLength::Extract_FCs_FixedWindowLength( bool acquisitionContinuous, int N, int winLngth ) :
		factorOfEachDecimationLevel(4),
		overlap(0.25),
		minNumDataPerFreq(20),
		pLRatio(0.05),
		firstsLinesToRepeatFromFcDistribution(4),
		fcDraftMadeForWinLengthOf(4096)
{
	isAcquisitionContinuous = acquisitionContinuous;
	nPointsTS = N;
	windowLength = winLngth;
	this->set_parameters();
}

void SetUpRemoteReferenceConcomitance::find_concomitance_for_stations( std::vector<StationBase> &station )
{
}

void Extract_FCs_FixedWindowLength::get_all_FCs( std::vector<StationBase> &station )
{
	size_t auxN = 0;
	bool auxAcqstnCntns = true;
	for( int istn = 0; istn < station.size(); istn++ ) {
		for( int its = 0; its < station[istn].ts.size(); its++ ) {

			auxN = station[istn].ts[its].ch[0].timeSeries[0].getDim(0);
			auxAcqstnCntns = station[istn].ts[its].isAcquisitionContinuous;
			Extract_FCs_FixedWindowLength fc( auxAcqstnCntns, auxN );

			fc.set_corrections( station[istn].ts[its] );

			for( int imatch = 0; imatch < station[istn].ts[its].amountOfPossibleCombinationsForRR; imatch++ ) {
			}
		}
	}
}

void Extract_FCs_FixedWindowLength::set_corrections( StationFile &ts )
{
	// auxiliary variables
	size_t npts2;
	double pi2n, t, g, freq, period, w, RE, IM;
	ComplexDouble temp;
	bool auxID;

	npts2 = floor(this->windowLength/2);
	pi2n = 2*PI/(this->windowLength/2);

	// fill in dr with sampling periods for all decimation levels
	std::vector<double> dr( this->nDecimationLevel );
	for( int i = 0;  i < dr.size(); i++ )
		dr[i] = pow( this->factorOfEachDecimationLevel, i )/ts.samplingFrequency;

	// correct for first difference if appropriate
    // correct unist of fc's
	for( int ich = 0; ich < ts.ch.size(); ich++ ) {
		for( int id = 0; id < this->nDecimationLevel; id++ ) {
			t = sqrt(dr[id]);
			for( int i = 0; i < npts2; i++ )
				ts.ch[ich].correction[0]( id, i ) = ts.ch[ich].countConversion*t;
		}
	}
	
	// corrections for low pass decimation filters 
	// corrects for effect of low pass filters on all decimation levels
	for( int ich = 0; ich < ts.ch.size(); ich++ ) {
		for( int id = 1; id < this->nDecimationLevel; id++ ) {
			for( int jd = id; jd < this->nDecimationLevel; jd++ ) {
				pi2n = 2*PI*dr[id - 1]/dr[jd]/this->windowLength;
				for( int i = 0; i < npts2; i++ ) {
					g = this->fcDistribution[0]( id, 1 );
					for( int j = 1; j < this->fcDistribution[0].getDim(1); j++ ) {
						t = pi2n*i*( j - 1 );
						g += 2*this->fcDistribution[0]( id, j )*cos( t );
					}
					ts.ch[ich].correction[0]( id, i ) = ts.ch[ich].correction[0]( id, i )/g;
				}
			}
		}
	}

	// corrections for analogue instrument filters/system response
	for( int ich = 0; ich < ts.ch.size(); ich++ ) {
		for( int id = 0; id < this->nDecimationLevel; id++ ) {
			for( int i = 0; i < npts2; i++ ) {
				freq = (i + 1)/dr[id]/this->windowLength;
				period = 1/freq;
				temp = ComplexDouble( 0, 0 );
				for( int j = 0; j < ts.ch[ich].correction[0].getDim(0) - 1; j++ ) {
					if( freq > ts.ch[ich].systemResponse[0]( j, 0 ).real() && freq <= ts.ch[ich].systemResponse[0]( j + 1, 0 ).real() ) {
						w = (ts.ch[ich].systemResponse[0]( j + 1, 0 ).real() - freq)/(ts.ch[ich].systemResponse[0]( j + 1, 0 ).real() - ts.ch[ich].systemResponse[0]( j, 0 ).real());
						RE = ts.ch[ich].systemResponse[0]( j, 1 ).real()*w + ts.ch[ich].systemResponse[0]( j + 1, 1 ).real()*( 1 - w );
						IM = ts.ch[ich].systemResponse[0]( j, 1 ).imag()*w + ts.ch[ich].systemResponse[0]( j + 1, 1 ).imag()*( 1 - w );
						auxID = true;
						break;
					}
				}
				if(auxID) {
					temp = ComplexDouble( RE, IM )*ts.ch[ich].gain;
					ts.ch[ich].correction[0]( id, i ) = ts.ch[ich].correction[0]( id, i )/temp;
				}
				else
					ts.ch[ich].correction[0]( id, i ) = 0;
			}
		}
	}
}

void Extract_FCs_FixedWindowLength::set_parameters()
{
	matCUDA::Array<int> auxFcDistribution = this->get_draft_of_fc_distribution();

	if( this->isAcquisitionContinuous == false )
		this->nDecimationLevel = 1;
	else
		//this->nDecimationLevel = determine_decimation_level( auxFcDistribution );
		this->nDecimationLevel = 1;

	matCUDA::Array<int> fcDistribution = this->get_fc_distribution( auxFcDistribution );

	this->nArCoeff = static_cast<int>(this->pLRatio*(double)this->windowLength/(1 - this->pLRatio) + 0.5);
}

matCUDA::Array<int> Extract_FCs_FixedWindowLength::get_draft_of_fc_distribution()
{
	Array<int> draftOfFcDistribution(13,2);

	draftOfFcDistribution(0,0) = 800; draftOfFcDistribution(0,1) = 1024;
	draftOfFcDistribution(1,0) = 544; draftOfFcDistribution(1,1) = 800;
	draftOfFcDistribution(2,0) = 384; draftOfFcDistribution(2,1) = 544;
	draftOfFcDistribution(3,0) = 256; draftOfFcDistribution(3,1) = 384;
	draftOfFcDistribution(4,0) = 192; draftOfFcDistribution(4,1) = 256;
	draftOfFcDistribution(5,0) = 128; draftOfFcDistribution(5,1) = 192;
	draftOfFcDistribution(6,0) = 96; draftOfFcDistribution(6,1) = 128;
	draftOfFcDistribution(7,0) = 64; draftOfFcDistribution(7,1) = 96;
	draftOfFcDistribution(8,0) = 32; draftOfFcDistribution(8,1) = 64;
	draftOfFcDistribution(9,0) = 24; draftOfFcDistribution(9,1) = 32;
	draftOfFcDistribution(10,0) = 16; draftOfFcDistribution(10,1) = 24;
	draftOfFcDistribution(11,0) = 8; draftOfFcDistribution(11,1) = 16;
	draftOfFcDistribution(12,0) = 4; draftOfFcDistribution(12,1) = 8;

	return draftOfFcDistribution;
}

matCUDA::Array<int> Extract_FCs_FixedWindowLength::get_fc_distribution( matCUDA::Array<int> auxFcDistribution )
{
	// TODO - improve
	size_t nElements = this->firstsLinesToRepeatFromFcDistribution*this->nDecimationLevel + auxFcDistribution.getDim(0) - this->firstsLinesToRepeatFromFcDistribution;
	size_t irow= -1;
	double changeFactor = (double)this->windowLength/(double)this->fcDraftMadeForWinLengthOf;
	std::cout << changeFactor << std::endl;

	Array<int> fcDistribution(nElements,3);
	for( int nDeci = 0; nDeci < this->nDecimationLevel; nDeci++ ) {
		for( int iLines = 0; iLines < this->firstsLinesToRepeatFromFcDistribution; iLines++ ) {
			irow++;
			fcDistribution(irow,0) = nDeci + 1;
			fcDistribution(irow,1) = static_cast<int>( auxFcDistribution(iLines,0)*changeFactor + 0.5 );
			fcDistribution(irow,2) = static_cast<int>( auxFcDistribution(iLines,1)*changeFactor + 0.5 );
		}
	}

	for( int iLines = this->firstsLinesToRepeatFromFcDistribution; iLines < auxFcDistribution.getDim(0); iLines++ ) {
		irow++;
		fcDistribution(irow,0) = this->nDecimationLevel;
		fcDistribution(irow,1) = static_cast<int>( auxFcDistribution(iLines,0)*changeFactor + 0.5 );
		fcDistribution(irow,2) = static_cast<int>( auxFcDistribution(iLines,1)*changeFactor + 0.5 );
	}

	return fcDistribution;
}

size_t Extract_FCs_FixedWindowLength::determine_decimation_level( matCUDA::Array<int> auxFcDistribution )
{
	size_t nDecimationLevel = 0;
	size_t amountOfFcsInLongerPeriod = auxFcDistribution( auxFcDistribution.getDim(0), 0 ) - auxFcDistribution( auxFcDistribution.getDim(0), 1 ) + 1;

	size_t nSeg = (size_t)floor(this->nPointsTS/(this->windowLength*(1 - this->overlap)));

	while(nSeg*amountOfFcsInLongerPeriod > this->minNumDataPerFreq) {
		nDecimationLevel++;
		nSeg = (size_t)floor(this->nPointsTS/pow(this->factorOfEachDecimationLevel,nDecimationLevel)/(this->windowLength*(1 - this->overlap)));
	}

	return nDecimationLevel;
}

void Extract_FCs_VariableWindowLength::get_all_FCs( std::vector<StationBase> &station )
{
}