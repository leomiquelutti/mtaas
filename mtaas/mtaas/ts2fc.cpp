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

			//// checkouts
			//if( station[istn].ts[its].mtu.mtuTsBand == phoenixTsBand_t(5) )
			//	for( int ich = 0; ich < station[istn].ts[its].ch.size(); ich++ )
			//		station[istn].ts[its].ch[ich].correction->print();

			for( int icomb = 0; icomb < station[istn].ts[its].combination.size(); icomb++ )
				fc.extract_fcs_for_each_combination( station, istn, its, icomb );
		}
	}
}

void Extract_FCs_FixedWindowLength::extract_fcs_for_each_combination( std::vector<StationBase> &station, const size_t idxStn, const size_t idxTs, const size_t idxCombination )
{
	size_t nSegments;
	size_t istn = idxStn, its = idxTs, icomb = idxCombination;

	Combination auxComb = station[istn].ts[its].combination[icomb];
	
	matCUDA::Array<double> auxReal( this->windowLength );
	matCUDA::Array<ComplexDouble> auxCplx( this->windowLength );

	for( int ideci = 0; ideci < this->nDecimationLevel; ideci++ ) {

		nSegments = floor( station[istn].ts[its].ch[0].timeSeries[0].getDim(0)/( this->windowLength*( 1 - this->overlap ) ) );

		for( int iseg = 0; iseg < nSegments; iseg++ ) {

			//// single site
			//if( station[istn].ts[its].combination[icomb].numberOfConcomitantTs == 1 ) {
			//}
			//// remote references
			//else
			//{
			//}

			// transform segments of TS in FC
			for( int iblock = 0; iblock < auxComb.idxBgn[0].getDim(0); iblock++ ) {
					
				size_t istart = auxComb.idxBgn[0](iblock,0), iend = istart + this->windowLength - 1;

				while( iend < auxComb.idxEnd[0](iblock,0) ) {

					// TEST IT NOW!!!!
					for( int ich = 0; ich < station[istn].ts[its].ch.size(); ich++ ) {
						auxReal = station[istn].ts[its].ch[ich].timeSeries[0].submatrix( istart, iend, 0, 0 );
						auxReal = auxReal.detrend();
						auxCplx = fft( &auxReal );
					}

					istart += floor( this->windowLength*( 1 - this->overlap ) );
					iend = istart + this->windowLength - 1;
				}
			}
						
			// transform segments of RR in FC

			// get in

			// get inRR

			// get out
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
		ts.ch[ich].correction = new Array<ComplexDouble>( this->nDecimationLevel, npts2 );
		for( int id = 0; id < this->nDecimationLevel; id++ ) {
			t = sqrt(dr[id]);
			for( int i = 0; i < npts2; i++ )
				ts.ch[ich].correction[0]( id, i ) = ts.ch[ich].countConversion*t;
		}
	}

	//// checkouts
	//if( ts.mtu.mtuTsBand == phoenixTsBand_t(5) ) {
	//	for( int ich = 0; ich < ts.ch.size(); ich++ ) {
	//		cout << ts.ch[ich].name << endl;
	//		cout << ts.ch[ich].correction[0]( 0, 0) << endl;
	//		//for( int i = 0; i < ts.ch[ich].correction[0].getDim(1) - 1; i++ ) {
	//		//	cout << ts.ch[ich].correction[0]( 0, i) << endl;
	//		//}
	//		//cout << endl;
	//	}
	//}
	
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

	//// checkouts
	//if( ts.mtu.mtuTsBand == phoenixTsBand_t(5) ) {
	//	for( int ich = 0; ich < ts.ch.size(); ich++ ) {
	//		cout << ts.ch[ich].name << endl;
	//		cout << ts.ch[ich].correction[0]( 0, 0) << endl;
	//		//for( int i = 0; i < ts.ch[ich].correction[0].getDim(1) - 1; i++ ) {
	//		//	cout << ts.ch[ich].correction[0]( 0, i) << endl;
	//		//}
	//		//cout << endl;
	//	}
	//}

	// corrections for analogue instrument filters/system response
	for( int ich = 0; ich < ts.ch.size(); ich++ ) {
		for( int id = 0; id < this->nDecimationLevel; id++ ) {
			for( int i = 0; i < npts2; i++ ) {
				freq = (i + 1)/dr[id]/this->windowLength;
				period = 1/freq;
				temp = ComplexDouble( 0, 0 );
				for( int j = 0; j < ts.ch[ich].systemResponseFreqs[0].getDim(0) - 1; j++ ) {
					if( freq > ts.ch[ich].systemResponseFreqs[0]( j + 1 ) && freq <= ts.ch[ich].systemResponseFreqs[0]( j ) ) {
						
						w = (ts.ch[ich].systemResponseFreqs[0]( j ) - freq)/(ts.ch[ich].systemResponseFreqs[0]( j ) - ts.ch[ich].systemResponseFreqs[0]( j + 1 ));
						RE = ts.ch[ich].systemResponse[0]( j + 1, 1 ).real()*w + ts.ch[ich].systemResponse[0]( j, 1 ).real()*( 1 - w );
						IM = ts.ch[ich].systemResponse[0]( j + 1, 1 ).imag()*w + ts.ch[ich].systemResponse[0]( j, 1 ).imag()*( 1 - w );
						
						auxID = true;
						break;
					}
				}
				if(auxID) {
					temp = ComplexDouble( RE, IM )*ts.ch[ich].gain;
					ts.ch[ich].correction[0]( id, i ) = ts.ch[ich].correction[0]( id, i )/temp;
					auxID = false;
				}
				else
					ts.ch[ich].correction[0]( id, i ) = 0;
			}
		}
	}

	//// checkouts
	//if( ts.mtu.mtuTsBand == phoenixTsBand_t(5) ) {
	//	for( int ich = 0; ich < ts.ch.size(); ich++ ) {
	//		//for( int i = 0; i < ts.ch[ich].correction[0].getDim(1) - 1; i++ )
	//			cout << ts.ch[ich].correction[0](0,1279) << " ";
	//	}
	//		cout << endl;
	//}
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