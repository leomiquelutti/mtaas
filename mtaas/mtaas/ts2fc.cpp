#include "ts2fc.h"
//#include <stdlib.h>
//#include <stdio.h>
//#include <math.h>
//#include <cmath>

// constructor
Extract_FCs_FixedWindowLength::Extract_FCs_FixedWindowLength( std::vector<StationBase> &station, int winLngth ) :
		factorOfEachDecimationLevel(4),
		overlap(0.25),
		minNumDataPerFreq(20),
		pLRatio(0.05),
		firstsLinesToRepeatFromFcDistribution(4),
		fcDraftMadeForWinLengthOf(4096),
		dpssNW(1),
		dpssDegree(0)
{
	windowLength = winLngth;
	this->set_parameters();
	this->get_all_FCs( station );
}

void Extract_FCs_FixedWindowLength::get_all_FCs( std::vector<StationBase> &station )
{
	// allocate memory for ch[ich].fc
	this->allocate_memory_for_ch_fc( station );

	// allocate memory for in, inRR and out for each FrequencyResponse
	this->allocate_memory_for_FrequencyResponses( station );

	//// checkouts
	//for( int istn = 0; istn < station.size(); istn++ )
	//	for( int its = 0; its < station[istn].ts.size(); its++ )
	//		station[istn].ts[its].combination[0].fcDistribution[0].print();

	// set corrections to be used
	this->set_corrections( station );

	// DO IT - retrieve FCs
	for( int istn = 0; istn < station.size(); istn++ )
		for( int its = 0; its < station[istn].ts.size(); its++ )
			for( int icomb = 0; icomb < station[istn].ts[its].combination.size(); icomb++ ) {
				// cout << istn << " " << its << " " << icomb << endl;
				this->extract_fcs_for_each_combination( station, istn, its, icomb );
				

				//// checkouts
				//for( int ifreq = 0; ifreq < station[istn].ts[its].combination[icomb].fr.size(); ifreq++ ) 
				//for( int i = 500; i < 505; i++ ) {	
				//	// in
				//	cout << station[istn].ts[its].combination[icomb].fr[ifreq].in[0]( i, 0 ) << " " << station[istn].ts[its].combination[icomb].fr[ifreq].in[0]( i, 1 ) << " ";
				//	// inRR
				//	
				//	cout << station[istn].ts[its].combination[icomb].fr[ifreq].inRR[0]( i, 0 ) << " " << station[istn].ts[its].combination[icomb].fr[ifreq].inRR[0]( i, 1 ) << " ";
				//	// out
				//	
				//	cout << station[istn].ts[its].combination[icomb].fr[ifreq].out[0]( i, 0 ) << " " << station[istn].ts[its].combination[icomb].fr[ifreq].out[0]( i, 1 ) << " " << station[istn].ts[its].combination[icomb].fr[ifreq].out[0]( i, 2 ) << endl;
				//	cout << endl << endl;
				//}
			}

	// free memory
	this->deallocate_memory_for_ch_fc( station );
}

//void Extract_FCs_FixedWindowLength::get_all_FCs( std::vector<StationBase> &station )
//{
//	Extract_FCs_FixedWindowLength fc;
//
//	// allocate memory for ch[ich].fc
//	fc.allocate_memory_for_ch_fc( station );
//
//	// allocate memory for in, inRR and out for each FrequencyResponse
//	fc.allocate_memory_for_FrequencyResponses( station );
//
//	//// checkouts
//	//for( int istn = 0; istn < station.size(); istn++ )
//	//	for( int its = 0; its < station[istn].ts.size(); its++ )
//	//		station[istn].ts[its].combination[0].fcDistribution[0].print();
//
//	// set corrections to be used
//	fc.set_corrections( station );
//
//	// DO IT - retrieve FCs
//	for( int istn = 0; istn < station.size(); istn++ )
//		for( int its = 0; its < station[istn].ts.size(); its++ )
//			for( int icomb = 0; icomb < station[istn].ts[its].combination.size(); icomb++ ) {
//				// cout << istn << " " << its << " " << icomb << endl;
//				fc.extract_fcs_for_each_combination( station, istn, its, icomb );
//				
//
//				//// checkouts
//				//for( int ifreq = 0; ifreq < station[istn].ts[its].combination[icomb].fr.size(); ifreq++ ) 
//				//for( int i = 500; i < 505; i++ ) {	
//				//	// in
//				//	cout << station[istn].ts[its].combination[icomb].fr[ifreq].in[0]( i, 0 ) << " " << station[istn].ts[its].combination[icomb].fr[ifreq].in[0]( i, 1 ) << " ";
//				//	// inRR
//				//	
//				//	cout << station[istn].ts[its].combination[icomb].fr[ifreq].inRR[0]( i, 0 ) << " " << station[istn].ts[its].combination[icomb].fr[ifreq].inRR[0]( i, 1 ) << " ";
//				//	// out
//				//	
//				//	cout << station[istn].ts[its].combination[icomb].fr[ifreq].out[0]( i, 0 ) << " " << station[istn].ts[its].combination[icomb].fr[ifreq].out[0]( i, 1 ) << " " << station[istn].ts[its].combination[icomb].fr[ifreq].out[0]( i, 2 ) << endl;
//				//	cout << endl << endl;
//				//}
//			}
//
//	// free memory
//	fc.deallocate_memory_for_ch_fc( station );
//}

//void Extract_FCs_FixedWindowLength::get_in( std::vector<StationBase> &station, const size_t iStnIn, const size_t iTsIn, const size_t iComb, const size_t idxSeg  )
//{
//	size_t idx = 0;
//	size_t lowerLimit = 0, upperLimit = 0;
//
//	for( int ifreq = 0; ifreq < station[iStnIn].ts[iTsIn].combination[iComb].fr.size(); ifreq++ ) {
//		
//		idx = idxSeg*station[iStnIn].ts[iTsIn].combination[iComb].fr[ifreq].nFCsPerSegment;
//		lowerLimit = station[iStnIn].ts[iTsIn].combination[iComb].fr[ifreq].lowerLimitBandWidth;
//		upperLimit = station[iStnIn].ts[iTsIn].combination[iComb].fr[ifreq].upperLimitBandWidth;
//		
//		for( int ifc = lowerLimit; ifc < upperLimit; ifc++ ) {
//			for( int ich = 0; ich < N_CHANNEL_INPUTS; ich++ )
//				station[iStnIn].ts[iTsIn].combination[iComb].fr[ifreq].in[0]( idx, ich ) = station[iStnIn].ts[iTsIn].ch[ich].fc[0]( ifc );
//			idx++;
//		}
//	}
//}
//
//void Extract_FCs_FixedWindowLength::get_out( std::vector<StationBase> &station, const size_t iStnIn, const size_t iTsIn, const size_t iComb, const size_t idxSeg  )
//{
//	size_t idx = 0;
//	size_t lowerLimit = 0, upperLimit = 0;
//
//	for( int ifreq = 0; ifreq < station[iStnIn].ts[iTsIn].combination[iComb].fr.size(); ifreq++ ) {
//		
//		idx = idxSeg*station[iStnIn].ts[iTsIn].combination[iComb].fr[ifreq].nFCsPerSegment;
//		lowerLimit = station[iStnIn].ts[iTsIn].combination[iComb].fr[ifreq].lowerLimitBandWidth;
//		upperLimit = station[iStnIn].ts[iTsIn].combination[iComb].fr[ifreq].upperLimitBandWidth;
//		
//		for( int ifc = lowerLimit; ifc < upperLimit; ifc++ ) {
//			for( int ich = N_CHANNEL_INPUTS; ich < station[iStnIn].ts[iTsIn].nChannels; ich++ )
//				station[iStnIn].ts[iTsIn].combination[iComb].fr[ifreq].out[0]( idx, ich - N_CHANNEL_INPUTS ) = station[iStnIn].ts[iTsIn].ch[ich].fc[0]( ifc );
//			idx++;
//		}
//	}
//}
//
//void Extract_FCs_FixedWindowLength::get_in_out_inRR( std::vector<StationBase> &station, const size_t iStnIn, const size_t iTsIn, const size_t iComb, const size_t iStnOut, const size_t iTsOut, const size_t idxSeg  )
//void Extract_FCs_FixedWindowLength::get_in_out_inRR( std::vector<StationBase> &station, const size_t iStnIn, const size_t iTsIn, const size_t iComb, const size_t idxSeg  )
//{
//	size_t idx = 0;
//	size_t lowerLimit = 0, upperLimit = 0;
//
//	for( int ifreq = 0; ifreq < station[iStnIn].ts[iTsIn].combination[iComb].fr.size(); ifreq++ ) {
//		
//		idx = idxSeg*station[iStnIn].ts[iTsIn].combination[iComb].fr[ifreq].nFCsPerSegment;
//		lowerLimit = station[iStnIn].ts[iTsIn].combination[iComb].fr[ifreq].lowerLimitBandWidth;
//		upperLimit = station[iStnIn].ts[iTsIn].combination[iComb].fr[ifreq].upperLimitBandWidth;
//		
//		for( int ifc = lowerLimit; ifc < upperLimit; ifc++ ) {
//			for( int ich = 0; ich < N_CHANNEL_INPUTS; ich++ ) {
//				station[iStnIn].ts[iTsIn].combination[iComb].fr[ifreq].in[0]( idx, ich ) = station[iStnIn].ts[iTsIn].ch[ich].fc[0]( ifc );
//				station[iStnIn].ts[iTsIn].combination[iComb].fr[ifreq].inRR[0]( idx, ich ) = station[iStnOut].ts[iTsOut].ch[ich].fc[0]( ifc );
//			}
//			for( int ich = N_CHANNEL_INPUTS; ich < station[iStnIn].ts[iTsIn].nChannels; ich++ )
//				station[iStnIn].ts[iTsIn].combination[iComb].fr[ifreq].out[0]( idx, ich - N_CHANNEL_INPUTS ) = station[iStnIn].ts[iTsIn].ch[ich].fc[0]( ifc );
//				
//			idx++;
//		}
//	}
//}

void Extract_FCs_FixedWindowLength::get_in_out_inRR( std::vector<StationBase> &station, const size_t istn, const size_t its, const size_t icomb, const size_t iseg  )
{
	size_t idx = 0;
	size_t lowerLimit = 0, upperLimit = 0;
	size_t istn2 = 0, its2 = 0, idxInRR = 0;

	//// checkouts
	//for( int i = 500; i < 505; i++ ) {
	//
	//	// in
	//	for( int ich = 0; ich < N_CHANNEL_INPUTS; ich++ )
	//		cout << station[istn].ts[its].ch[ich].fc[0]( i ) << " ";
	//
	//	// out
	//	for( int ich = N_CHANNEL_INPUTS; ich < station[istn].ts[its].ch.size(); ich++ )
	//		cout << station[istn].ts[its].ch[ich].fc[0]( i ) << " ";		
	//	
	//	// inRR
	//	//cout << station[istn].ts[its].combination[icomb].nConcomitantTs << endl;
	//	for( int iaux = 0; iaux < station[istn].ts[its].combination[icomb].nConcomitantTs; iaux++ ) {
	//
	//		istn2 = station[istn].ts[its].combination[icomb].idxStn[iaux];
	//		its2 = station[istn].ts[its].combination[icomb].idxTs[iaux];
	//		
	//		for( int ich = 0; ich < N_CHANNEL_INPUTS; ich++ )
	//			cout << station[istn2].ts[its2].ch[ich].fc[0]( i ) << " ";
	//	}
	//	cout << endl << endl;
	//}

	for( int ifreq = 0; ifreq < station[istn].ts[its].combination[icomb].fr.size(); ifreq++ ) {
		
		idx = iseg*station[istn].ts[its].combination[icomb].fr[ifreq].nFCsPerSegment;
		lowerLimit = station[istn].ts[its].combination[icomb].fr[ifreq].lowerLimitBandWidth;
		upperLimit = station[istn].ts[its].combination[icomb].fr[ifreq].upperLimitBandWidth;
			
		// in
		for( int ich = 0; ich < N_CHANNEL_INPUTS; ich++ ) {

			//for( int i = 500; i < 505; i++ )
			//	cout << station[istn].ts[its].ch[ich].fc[0]( i ) << endl;
			//cout << endl;
			//station[istn].ts[its].ch[ich].fc[0].print();
			idx = iseg*station[istn].ts[its].combination[icomb].fr[ifreq].nFCsPerSegment;
			for( int ifc = lowerLimit; ifc < upperLimit; ifc++ )
				station[istn].ts[its].combination[icomb].fr[ifreq].in[0]( idx++, ich ) = station[istn].ts[its].ch[ich].fc[0]( ifc );
		}
			
		//out
		for( int ich = N_CHANNEL_INPUTS; ich < station[istn].ts[its].ch.size(); ich++ ) {

			//for( int i = 500; i < 505; i++ )
			//	cout << station[istn].ts[its].ch[ich].fc[0]( i ) << endl;
			//cout << endl;
			//station[istn].ts[its].ch[ich].fc[0].print();
			idx = iseg*station[istn].ts[its].combination[icomb].fr[ifreq].nFCsPerSegment;
			for( int ifc = lowerLimit; ifc < upperLimit; ifc++ )
				station[istn].ts[its].combination[icomb].fr[ifreq].out[0]( idx++, ich - N_CHANNEL_INPUTS ) = station[istn].ts[its].ch[ich].fc[0]( ifc );
		}
			
		//// inRR
		//idxInRR = 0;
		//for( int iaux = 0; iaux < station[istn].ts[its].combination[icomb].nConcomitantTs + 1; iaux++ ) {
		//
		//	if( iaux == 0 ) {
		//		istn2 = istn;
		//		its2 = its;
		//	}
		//	else {
		//		istn2 = station[istn].ts[its].combination[icomb].idxStn[iaux - 1];
		//		its2 = station[istn].ts[its].combination[icomb].idxTs[iaux - 1];
		//	}
		//	
		//	for( int ich = 0; ich < N_CHANNEL_INPUTS; ich++ ) {
		//		idx = iseg*station[istn].ts[its].combination[icomb].fr[ifreq].nFCsPerSegment;
		//		cout << idx << " " << idxInRR << " " << istn << " " << its << " " << istn2 << " " << its2 << endl;
		//		for( int ifc = lowerLimit; ifc < upperLimit; ifc++ )
		//			station[istn].ts[its].combination[icomb].fr[ifreq].inRR[0]( idx++, ich + idxInRR ) = station[istn2].ts[its2].ch[ich].fc[0]( ifc );
		//	}
		//
		//	idxInRR += N_CHANNEL_INPUTS;
		//}		

		// inRR
		idxInRR = 0;
		for( int iaux = 0; iaux < station[istn].ts[its].combination[icomb].nConcomitantTs; iaux++ ) {

			istn2 = station[istn].ts[its].combination[icomb].idxStn[iaux];
			its2 = station[istn].ts[its].combination[icomb].idxTs[iaux];
			
			for( int ich = 0; ich < N_CHANNEL_INPUTS; ich++ ) {

				//for( int i = 500; i < 505; i++ )
				//	cout << station[istn2].ts[its2].ch[ich].fc[0]( i ) << endl;
				//cout << endl;
				//station[istn2].ts[its2].ch[ich].fc[0].print();
				idx = iseg*station[istn].ts[its].combination[icomb].fr[ifreq].nFCsPerSegment;
				for( int ifc = lowerLimit; ifc < upperLimit; ifc++ )
					station[istn].ts[its].combination[icomb].fr[ifreq].inRR[0]( idx++, ich + idxInRR ) = station[istn2].ts[its2].ch[ich].fc[0]( ifc );
			}

			idxInRR += N_CHANNEL_INPUTS;
		}		

		idx++;
		
		//for( int ifc = lowerLimit; ifc < upperLimit; ifc++ ) {
		//	
		//	// in
		//	for( int ich = 0; ich < N_CHANNEL_INPUTS; ich++ )
		//		station[istn].ts[its].combination[icomb].fr[ifreq].in[0]( idx, ich ) = station[istn].ts[its].ch[ich].fc[0]( ifc );
		//	
		//	// out
		//	for( int ich = N_CHANNEL_INPUTS; ich < station[istn].ts[its].ch.size(); ich++ )
		//		station[istn].ts[its].combination[icomb].fr[ifreq].out[0]( idx, ich - N_CHANNEL_INPUTS ) = station[istn].ts[its].ch[ich].fc[0]( ifc );
		//	
		//	// inRR
		//	idxInRR = 0;
		//	for( int iaux = 0; iaux = station[istn].ts[its].combination[icomb].nConcomitantTs + 1; iaux++ ) {
		//
		//		if( iaux == 0 ) {
		//			istn2 = istn;
		//			its2 = its;
		//		}
		//		else {
		//			istn2 = station[istn].ts[its].combination[icomb].idxStn[iaux - 1];
		//			its2 = station[istn].ts[its].combination[icomb].idxTs[iaux - 1];
		//			upperLimit = N_CHANNEL_INPUTS;
		//		}
		//
		//		for( int ich = 0; ich < N_CHANNEL_INPUTS; ich++ )
		//			station[istn].ts[its].combination[icomb].fr[ifreq].inRR[0]( idx, ich + idxInRR ) = station[istn2].ts[its2].ch[ich].fc[0]( ifc );
		//		idxInRR += N_CHANNEL_INPUTS;
		//	}		
		//
		//	idx++;
		//}
	}
}

void Extract_FCs_FixedWindowLength::extract_fcs_for_each_combination( std::vector<StationBase> &station, const size_t idxStn, const size_t idxTs, const size_t idxCombination )
{
	size_t nSegments;
	size_t istn = idxStn, its = idxTs, icomb = idxCombination; // index just to help
	size_t istart = 0, iend = istart + this->windowLength - 1; // index for reading concomitant TSs
	size_t winPace = floor( this->windowLength*( 1 - this->overlap ) );

	Combination auxComb = station[istn].ts[its].combination[icomb];
	
	matCUDA::Array<double> auxReal( this->windowLength );		
					
	matCUDA::Array<ComplexDouble> auxCmplx( floor(this->windowLength/2) + 1 );

	// index for writing in, inRR and out variables
	size_t idx = 0;

	// sweep through the "blocks" of data - one block is a piecewise continuous set of data
	for( int iblock = 0; iblock < auxComb.idxBgn[0].getDim(0); iblock++ ) {

		// sweep through the decimation levels of each block
		for( int ideci = 0; ideci < auxComb.deciLevelEachBlock[iblock]; ideci++ ) {

			nSegments = floor( ( auxComb.idxEnd[0](iblock,0) - auxComb.idxBgn[0](iblock,0) )/( this->windowLength*( 1 - this->overlap ) ) );
			//cout << nSegments << endl;

			// sweep through the segments of each decimation level of each block
			for( int iseg = 0; iseg < nSegments; iseg++ ) {

				//cout << iseg << " " << nSegments << endl << endl;
				for( int iInvolvedTs = 0; iInvolvedTs < auxComb.nConcomitantTs + 1; iInvolvedTs++ ) {
					
					istart = auxComb.idxBgn[0](iblock,iInvolvedTs) + iseg*winPace;
					iend = istart + this->windowLength - 1;

					size_t iAuxStn = 0, iAuxTs = 0, upperLimit = 0;

					// define which will be treated now
					if( iInvolvedTs == 0 ) {
						iAuxStn = istn;
						iAuxTs = its;
						upperLimit = station[iAuxStn].ts[iAuxTs].ch.size();
					}
					else {
						iAuxStn = auxComb.idxStn[iInvolvedTs - 1];
						iAuxTs = auxComb.idxTs[iInvolvedTs - 1];
						upperLimit = N_CHANNEL_INPUTS;
					}

					// get FCs for stations of interest
					for( int ich = 0; ich < upperLimit; ich++ ) {
						
						size_t imin = 500, imax = 505;
						auxReal = station[iAuxStn].ts[iAuxTs].ch[ich].timeSeries[0].submatrix( istart, iend, 0, 0 );
						
						//for( int i = imin; i < imax; i++ )
						//	cout << auxReal( i ) << endl;
						//cout<< endl;

						auxReal = auxReal.detrend();
						
						//for( int i = imin; i < imax; i++ )
						//	cout << auxReal( i ) << endl;
						//cout<< endl;

						auxReal = auxReal.elementWiseMultiply( this->auxDpss );
						
						//for( int i = imin; i < imax; i++ )
						//	cout << auxReal( i ) << endl;
						//cout<< endl;

						auxCmplx = fft( &auxReal );
						
						//for( int i = imin; i < imax; i++ )
						//	cout << auxCmplx( i ) << endl;
						//cout<< endl;

						this->correct_fcs( &auxCmplx, &station[iAuxStn].ts[iAuxTs].ch[ich].correction[0], ideci );
						
						//for( int i = imin; i < imax; i++ )
						//	cout << auxCmplx( i ) << endl;
						//cout<< endl;

						//station[iAuxStn].ts[iAuxTs].ch[ich].fc[0].print();

						//cout << iAuxStn << " " << iAuxTs << " " << ich << " " << endl;
						station[iAuxStn].ts[iAuxTs].ch[ich].fc[0] = auxCmplx.submatrix( 0, this->nFreq - 1, 0, 0 );

						//station[iAuxStn].ts[iAuxTs].ch[ich].fc[0].print();
						
						//for( int i = imin; i < imax; i++ )
						//	cout << station[iAuxStn].ts[iAuxTs].ch[ich].fc[0]( i ) << " " << auxCmplx(i) << endl;
						//cout<< endl;
					}
				}
				//cout << iseg << endl << endl;
				this->get_in_out_inRR( station, istn, its, icomb, idx++ );
			}
	//station[0].ts[0].ch[0].fc[0].print();
		}
	}
}

void Extract_FCs_FixedWindowLength::allocate_memory_for_FrequencyResponses( std::vector<StationBase> &station )
{
	size_t nSeg = 0, nDeci = 1, lngth = 0, nDeciMax = 0, nFreq = 0, nFCsPerSegment = 0, nElements = 0;
	double auxSmplFreq = 0, auxMeasFreq = 0;

	for( int istn = 0; istn < station.size(); istn++ ) {
		for( int its = 0; its < station[istn].ts.size(); its++ ) {
			for( int icomb = 0; icomb < station[istn].ts[its].combination.size(); icomb++ ) {
				
				std::vector<int> idxFreqsToRemove, idxFreqsToAdd;

				for( int iblock = 0; iblock < station[istn].ts[its].combination[icomb].nBlocks; iblock++ ) {

					// amount of sequencial points in block
					lngth = station[istn].ts[its].combination[icomb].idxEnd[0](iblock,0) - station[istn].ts[its].combination[icomb].idxBgn[0](iblock,0);
					station[istn].ts[its].combination[icomb].lengthTsEachBlock.push_back( lngth );

					// maximum decimation level allowed for the block
					nDeci = this->determine_decimation_level( lngth );
					nDeci = 1; // COMMENT THIS LATER WHEN PROGRAM IS WORKING WITH DECIMATION
					station[istn].ts[its].combination[icomb].deciLevelEachBlock.push_back( nDeci );
				}

				// calculate fcDist of the Combination
				nDeciMax = *std::max_element( station[istn].ts[its].combination[icomb].deciLevelEachBlock.begin(), station[istn].ts[its].combination[icomb].deciLevelEachBlock.end() );
				matCUDA::Array<int> auxFc = this->get_fc_distribution( nDeciMax );

				// calculate sampling frequencies of all decimation levels for this Combination
				station[istn].ts[its].combination[icomb].samplFreqEachDeciLevel.resize( nDeciMax );
				for( int ideci = 0; ideci < nDeciMax; ideci++ )
					station[istn].ts[its].combination[icomb].samplFreqEachDeciLevel[ideci] = station[istn].ts[its].samplingFrequency/pow( this->factorOfEachDecimationLevel, ideci );

				// calculate measured frequencies of the Combination
				nFreq = auxFc.getDim(0);
				station[istn].ts[its].combination[icomb].fr.resize(nFreq);
				station[istn].ts[its].combination[icomb].measuredFrequencies.resize(nFreq);
				for( int ifreq = 0; ifreq < nFreq; ifreq++ ) {
					auxSmplFreq = station[istn].ts[its].combination[icomb].samplFreqEachDeciLevel[auxFc(ifreq,0) - 1];
					auxMeasFreq = auxSmplFreq/2/this->windowLength*( auxFc(ifreq,1) + auxFc(ifreq,2) );
					station[istn].ts[its].combination[icomb].measuredFrequencies[ifreq] = auxMeasFreq;
					
					// fill in FrequencyResponse parameters
					station[istn].ts[its].combination[icomb].fr[ifreq].frequency = auxMeasFreq;
					station[istn].ts[its].combination[icomb].fr[ifreq].deciLevel = auxFc( ifreq, 0 );
					station[istn].ts[its].combination[icomb].fr[ifreq].lowerLimitBandWidth = auxFc( ifreq, 1 );
					station[istn].ts[its].combination[icomb].fr[ifreq].upperLimitBandWidth = auxFc( ifreq, 2 );
					station[istn].ts[its].combination[icomb].fr[ifreq].nFCsPerSegment = auxFc( ifreq, 2 ) - auxFc( ifreq, 1 ) + 1;

					if( auxMeasFreq > station[istn].ts[its].maxFreqInHertz || auxMeasFreq < station[istn].ts[its].minFreqInHertz )
						idxFreqsToRemove.push_back( ifreq );
					else
						idxFreqsToAdd.push_back( ifreq );
				}

				// rebuild fcDist
				nFreq = idxFreqsToAdd.size();
				station[istn].ts[its].combination[icomb].fcDistribution = new matCUDA::Array<int> ( nFreq, auxFc.getDim(1) );
				for( int ifreq = 0; ifreq < nFreq; ifreq++ )
					for( int icol = 0; icol < auxFc.getDim(1); icol++ )
						station[istn].ts[its].combination[icomb].fcDistribution[0]( ifreq, icol ) = auxFc( idxFreqsToAdd[ifreq], icol );
				
				//// checkouts				
				//station[istn].ts[its].ch[0].systemResponseFreqs[0].print();
				//cout << endl;	
				//auxFc.print();
				//cout << endl;
				//station[istn].ts[its].combination[icomb].fcDistribution[0].print();
				//cout << endl << idxFreqsToRemove.size() << endl << endl;
				//for( int i = 0; i < idxFreqsToRemove.size(); i++ )
				//	cout << idxFreqsToRemove[i] << " ";
				//cout << endl << endl;

				// remove frequencies not allowed by system response 
				for( int ifreq = idxFreqsToRemove.size() - 1; ifreq >= 0; ifreq--) {
					station[istn].ts[its].combination[icomb].measuredFrequencies.erase( station[istn].ts[its].combination[icomb].measuredFrequencies.begin() + idxFreqsToRemove[ifreq] );
					station[istn].ts[its].combination[icomb].fr.erase( station[istn].ts[its].combination[icomb].fr.begin() + idxFreqsToRemove[ifreq] );
				}

				// set new maximum decimation level
				nDeciMax = station[istn].ts[its].combination[icomb].fcDistribution[0]( station[istn].ts[its].combination[icomb].fcDistribution[0].getDim(0) - 1, 0 );

				// rebuild samplFreqEachDeciLevel
				size_t auxMaxDeci = station[istn].ts[its].combination[icomb].samplFreqEachDeciLevel.size();
				for( int ideci = auxMaxDeci - 1; ideci >= nDeciMax; ideci-- )
					station[istn].ts[its].combination[icomb].samplFreqEachDeciLevel.erase( station[istn].ts[its].combination[icomb].samplFreqEachDeciLevel.begin() + ideci );

				// put info in vector station
				station[istn].ts[its].combination[icomb].nDeci = nDeciMax;
				station[istn].ts[its].combination[icomb].nFreq = nFreq;

				// allocate memory for "in", "inRR" and "out" of each FrequencyResponses
				for( int ifreq = 0; ifreq < nFreq; ifreq++ ) {
					for( int iblock = 0; iblock < station[istn].ts[its].combination[icomb].nBlocks; iblock++ ) {
						nDeci = station[istn].ts[its].combination[icomb].fr[ifreq].deciLevel;
						nSeg = (size_t)floor( lngth/pow(this->factorOfEachDecimationLevel,nDeci - 1)/( this->windowLength*( 1 - this->overlap ) ) );
						nFCsPerSegment = station[istn].ts[its].combination[icomb].fr[ifreq].nFCsPerSegment;
						station[istn].ts[its].combination[icomb].fr[ifreq].elementCounter += nSeg*nFCsPerSegment;
					}

					nElements = station[istn].ts[its].combination[icomb].fr[ifreq].elementCounter;
					//cout << nElements << " " << station[istn].ts[its].combination[icomb].fr[ifreq].frequency << " " << station[istn].ts[its].combination[icomb].fr[ifreq].deciLevel << endl << endl;
					station[istn].ts[its].combination[icomb].fr[ifreq].in = new matCUDA::Array<ComplexDouble>( nElements, N_CHANNEL_INPUTS );
					station[istn].ts[its].combination[icomb].fr[ifreq].out = new matCUDA::Array<ComplexDouble>( nElements, station[istn].ts[its].nChannels -  N_CHANNEL_INPUTS );
					station[istn].ts[its].combination[icomb].fr[ifreq].inRR = new matCUDA::Array<ComplexDouble>( nElements, N_CHANNEL_INPUTS*station[istn].ts[its].combination[icomb].nConcomitantTs );
				}
			}
			station[istn].ts[its].maxDecimationLevel = *std::max_element( station[istn].ts[its].combination[0].deciLevelEachBlock.begin(), station[istn].ts[its].combination[ station[istn].ts[its].combination.size() - 1 ].deciLevelEachBlock.end() );
		}
	}
}

void Extract_FCs_FixedWindowLength::allocate_memory_for_ch_fc( std::vector<StationBase> &station )
{
	//cout<< this->nFreq << endl;
	for( int istn = 0; istn < station.size(); istn++ ) 
		for( int its = 0; its < station[istn].ts.size(); its++ ) 
			for( int ich = 0; ich < station[istn].ts[its].ch.size(); ich++ )
				station[istn].ts[its].ch[ich].fc = new matCUDA::Array<ComplexDouble>(this->nFreq);
}

void Extract_FCs_FixedWindowLength::deallocate_memory_for_ch_fc( std::vector<StationBase> &station )
{
	for( int istn = 0; istn < station.size(); istn++ ) 
		for( int its = 0; its < station[istn].ts.size(); its++ ) 
			for( int ich = 0; ich < station[istn].ts[its].ch.size(); ich++ )
				delete station[istn].ts[its].ch[ich].fc;
}

void Extract_FCs_FixedWindowLength::correct_fcs( matCUDA::Array<ComplexDouble> *data, matCUDA::Array<ComplexDouble> *correction, size_t deciLevel )
{
	//for( int i = 1990; i < 2000; i++ )
	//	cout << correction[0]( i, deciLevel ) << endl;
	//cout<< endl;
	//cout<< this->nFreq << endl;
	//correction[0].print();
	for( int i = 0; i < this->nFreq; i++ )
		data[0]( i ) = data[0]( i + 1 )*correction[0]( deciLevel, i);
	//for( int i = 500; i < 505; i++ )
	//	cout << data[0]( i ) << endl;
	//cout<< endl;
}

void Extract_FCs_FixedWindowLength::set_corrections( std::vector<StationBase> &station )
{
	// auxiliary variables
	size_t npts2;
	double pi2n, t, g, freq, period, w, RE, IM;
	ComplexDouble temp;
	bool auxID;

	npts2 = floor(this->windowLength/2);
	pi2n = 2*PI/(this->windowLength/2);

	for( int istn = 0; istn < station.size(); istn++ ) {
		for( int its = 0; its < station[istn].ts.size(); its++ ) {

			// fill in dr with sampling periods for all decimation levels
			//cout << station[istn].ts[its].maxDecimationLevel << endl;
			std::vector<double> dr( station[istn].ts[its].maxDecimationLevel );
			for( int i = 0;  i < dr.size(); i++ )
				dr[i] = pow( this->factorOfEachDecimationLevel, i )/station[istn].ts[its].samplingFrequency;

			// correct for first difference if appropriate
			// correct unist of fc's
			for( int ich = 0; ich < station[istn].ts[its].ch.size(); ich++ ) {
				station[istn].ts[its].ch[ich].correction = new matCUDA::Array<ComplexDouble>( station[istn].ts[its].maxDecimationLevel, npts2 );
				for( int id = 0; id < station[istn].ts[its].maxDecimationLevel; id++ ) {
					t = sqrt(dr[id]);
					for( int i = 0; i < npts2; i++ )
						station[istn].ts[its].ch[ich].correction[0]( id, i ) = station[istn].ts[its].ch[ich].countConversion*t;
				}
			}

			//// checkouts
			//if( station[istn].ts[its].mtu.mtuTsBand == phoenixTsBand_t(5) ) {
			//	for( int ich = 0; ich < station[istn].ts[its].ch.size(); ich++ ) {
			//		cout << station[istn].ts[its].ch[ich].name << endl;
			//		cout << station[istn].ts[its].ch[ich].correction[0]( 0, 0) << endl;
			//		//for( int i = 0; i < station[istn].ts[its].ch[ich].correction[0].getDim(1) - 1; i++ ) {
			//		//	cout << station[istn].ts[its].ch[ich].correction[0]( 0, i) << endl;
			//		//}
			//		//cout << endl;
			//	}
			//}
	
			//// corrections for low pass decimation filters 
			//// corrects for effect of low pass filters on all decimation levels
			//this->fcDistribution[0].print();
			//for( int ich = 0; ich < station[istn].ts[its].ch.size(); ich++ ) {
			//	for( int id = 1; id < station[istn].ts[its].maxDecimationLevel; id++ ) {
			//		for( int jd = id; jd < station[istn].ts[its].maxDecimationLevel; jd++ ) {
			//			pi2n = 2*PI*dr[id - 1]/dr[jd]/this->windowLength;
			//			for( int i = 0; i < npts2; i++ ) {
			//				g = this->fcDistribution[0]( id, 1 );
			//				for( int j = 1; j < this->fcDistribution[0].getDim(1); j++ ) {
			//					t = pi2n*i*( j - 1 );
			//					g += 2*this->fcDistribution[0]( id, j )*cos( t );
			//				}
			//				station[istn].ts[its].ch[ich].correction[0]( id, i ) = station[istn].ts[its].ch[ich].correction[0]( id, i )/g;
			//			}
			//		}
			//	}
			//}

			//// checkouts
			//if( station[istn].ts[its].mtu.mtuTsBand == phoenixTsBand_t(5) ) {
			//	for( int ich = 0; ich < station[istn].ts[its].ch.size(); ich++ ) {
			//		cout << station[istn].ts[its].ch[ich].name << endl;
			//		cout << station[istn].ts[its].ch[ich].correction[0]( 0, 0) << endl;
			//		//for( int i = 0; i < station[istn].ts[its].ch[ich].correction[0].getDim(1) - 1; i++ ) {
			//		//	cout << station[istn].ts[its].ch[ich].correction[0]( 0, i) << endl;
			//		//}
			//		//cout << endl;
			//	}
			//}

			// corrections for analogue instrument filters/system response
			for( int ich = 0; ich < station[istn].ts[its].ch.size(); ich++ ) {
				for( int id = 0; id < station[istn].ts[its].maxDecimationLevel; id++ ) {
					for( int i = 0; i < npts2; i++ ) {
						freq = (i + 1)/dr[id]/this->windowLength;
						period = 1/freq;
						temp = ComplexDouble( 0, 0 );
						for( int j = 0; j < station[istn].ts[its].ch[ich].systemResponseFreqs[0].getDim(0) - 1; j++ ) {
							if( freq > station[istn].ts[its].ch[ich].systemResponseFreqs[0]( j + 1 ) && freq <= station[istn].ts[its].ch[ich].systemResponseFreqs[0]( j ) ) {
						
								w = (station[istn].ts[its].ch[ich].systemResponseFreqs[0]( j ) - freq)/(station[istn].ts[its].ch[ich].systemResponseFreqs[0]( j ) - station[istn].ts[its].ch[ich].systemResponseFreqs[0]( j + 1 ));
								RE = station[istn].ts[its].ch[ich].systemResponse[0]( j + 1, 1 ).real()*w + station[istn].ts[its].ch[ich].systemResponse[0]( j, 1 ).real()*( 1 - w );
								IM = station[istn].ts[its].ch[ich].systemResponse[0]( j + 1, 1 ).imag()*w + station[istn].ts[its].ch[ich].systemResponse[0]( j, 1 ).imag()*( 1 - w );
						
								auxID = true;
								break;
							}
						}
						if(auxID) {
							temp = ComplexDouble( RE, IM )*station[istn].ts[its].ch[ich].gain;
							station[istn].ts[its].ch[ich].correction[0]( id, i ) = station[istn].ts[its].ch[ich].correction[0]( id, i )/temp;
							auxID = false;
						}
						else
							station[istn].ts[its].ch[ich].correction[0]( id, i ) = 0;
					}
				}
			}
		}
	}

	//// checkouts
	//if( ts.mtu.mtuTsBand == phoenixTsBand_t(5) ) {
	//	for( int ich = 0; ich < ts.ch.size(); ich++ ) {
	//		//for( int i = 0; i < station[istn].ts[its].ch[ich].correction[0].getDim(1) - 1; i++ )
	//			cout << station[istn].ts[its].ch[ich].correction[0](0,1279) << " ";
	//	}
	//		cout << endl;
	//}
}

//void Extract_FCs_FixedWindowLength::set_corrections( StationFile &ts )
//{
//	// auxiliary variables
//	size_t npts2;
//	double pi2n, t, g, freq, period, w, RE, IM;
//	ComplexDouble temp;
//	bool auxID;
//
//	npts2 = floor(this->windowLength/2);
//	pi2n = 2*PI/(this->windowLength/2);
//
//	// fill in dr with sampling periods for all decimation levels
//	std::vector<double> dr( this->nDecimationLevel );
//	for( int i = 0;  i < dr.size(); i++ )
//		dr[i] = pow( this->factorOfEachDecimationLevel, i )/ts.samplingFrequency;
//
//	// correct for first difference if appropriate
//    // correct unist of fc's
//	for( int ich = 0; ich < ts.ch.size(); ich++ ) {
//		ts.ch[ich].correction = new Array<ComplexDouble>( this->nDecimationLevel, npts2 );
//		for( int id = 0; id < this->nDecimationLevel; id++ ) {
//			t = sqrt(dr[id]);
//			for( int i = 0; i < npts2; i++ )
//				ts.ch[ich].correction[0]( id, i ) = ts.ch[ich].countConversion*t;
//		}
//	}
//
//	//// checkouts
//	//if( ts.mtu.mtuTsBand == phoenixTsBand_t(5) ) {
//	//	for( int ich = 0; ich < ts.ch.size(); ich++ ) {
//	//		cout << ts.ch[ich].name << endl;
//	//		cout << ts.ch[ich].correction[0]( 0, 0) << endl;
//	//		//for( int i = 0; i < ts.ch[ich].correction[0].getDim(1) - 1; i++ ) {
//	//		//	cout << ts.ch[ich].correction[0]( 0, i) << endl;
//	//		//}
//	//		//cout << endl;
//	//	}
//	//}
//	
//	// corrections for low pass decimation filters 
//	// corrects for effect of low pass filters on all decimation levels
//	for( int ich = 0; ich < ts.ch.size(); ich++ ) {
//		for( int id = 1; id < this->nDecimationLevel; id++ ) {
//			for( int jd = id; jd < this->nDecimationLevel; jd++ ) {
//				pi2n = 2*PI*dr[id - 1]/dr[jd]/this->windowLength;
//				for( int i = 0; i < npts2; i++ ) {
//					g = this->fcDistribution[0]( id, 1 );
//					for( int j = 1; j < this->fcDistribution[0].getDim(1); j++ ) {
//						t = pi2n*i*( j - 1 );
//						g += 2*this->fcDistribution[0]( id, j )*cos( t );
//					}
//					ts.ch[ich].correction[0]( id, i ) = ts.ch[ich].correction[0]( id, i )/g;
//				}
//			}
//		}
//	}
//
//	//// checkouts
//	//if( ts.mtu.mtuTsBand == phoenixTsBand_t(5) ) {
//	//	for( int ich = 0; ich < ts.ch.size(); ich++ ) {
//	//		cout << ts.ch[ich].name << endl;
//	//		cout << ts.ch[ich].correction[0]( 0, 0) << endl;
//	//		//for( int i = 0; i < ts.ch[ich].correction[0].getDim(1) - 1; i++ ) {
//	//		//	cout << ts.ch[ich].correction[0]( 0, i) << endl;
//	//		//}
//	//		//cout << endl;
//	//	}
//	//}
//
//	// corrections for analogue instrument filters/system response
//	for( int ich = 0; ich < ts.ch.size(); ich++ ) {
//		for( int id = 0; id < this->nDecimationLevel; id++ ) {
//			for( int i = 0; i < npts2; i++ ) {
//				freq = (i + 1)/dr[id]/this->windowLength;
//				period = 1/freq;
//				temp = ComplexDouble( 0, 0 );
//				for( int j = 0; j < ts.ch[ich].systemResponseFreqs[0].getDim(0) - 1; j++ ) {
//					if( freq > ts.ch[ich].systemResponseFreqs[0]( j + 1 ) && freq <= ts.ch[ich].systemResponseFreqs[0]( j ) ) {
//						
//						w = (ts.ch[ich].systemResponseFreqs[0]( j ) - freq)/(ts.ch[ich].systemResponseFreqs[0]( j ) - ts.ch[ich].systemResponseFreqs[0]( j + 1 ));
//						RE = ts.ch[ich].systemResponse[0]( j + 1, 1 ).real()*w + ts.ch[ich].systemResponse[0]( j, 1 ).real()*( 1 - w );
//						IM = ts.ch[ich].systemResponse[0]( j + 1, 1 ).imag()*w + ts.ch[ich].systemResponse[0]( j, 1 ).imag()*( 1 - w );
//						
//						auxID = true;
//						break;
//					}
//				}
//				if(auxID) {
//					temp = ComplexDouble( RE, IM )*ts.ch[ich].gain;
//					ts.ch[ich].correction[0]( id, i ) = ts.ch[ich].correction[0]( id, i )/temp;
//					auxID = false;
//				}
//				else
//					ts.ch[ich].correction[0]( id, i ) = 0;
//			}
//		}
//	}
//
//	//// checkouts
//	//if( ts.mtu.mtuTsBand == phoenixTsBand_t(5) ) {
//	//	for( int ich = 0; ich < ts.ch.size(); ich++ ) {
//	//		//for( int i = 0; i < ts.ch[ich].correction[0].getDim(1) - 1; i++ )
//	//			cout << ts.ch[ich].correction[0](0,1279) << " ";
//	//	}
//	//		cout << endl;
//	//}
//}

void Extract_FCs_FixedWindowLength::set_parameters()
{
	this->get_draft_of_fc_distribution();

	this->nArCoeff = static_cast<int>(this->pLRatio*(double)this->windowLength/(1 - this->pLRatio) + 0.5);

	this->auxDpss = new matCUDA::Array<double>( matCUDA::read_file_vector<double>( std::string("dpss_4096_1_0.txt") ) );
	//this->auxDpss = &matCUDA::read_file_vector<double>( std::string("dpss_4096_1_0.txt") );
	
	//if( this->isAcquisitionContinuous == false )
	//	this->nDecimationLevel = 1;
	//else
	//	//this->nDecimationLevel = determine_decimation_level( this->auxFcDistribution );
	//	this->nDecimationLevel = 1;
	//
	//this->nDecimationLevel = determine_decimation_level( this->auxFcDistribution );
	//this->nDecimationLevel = 1;
	//
	//matCUDA::Array<int> fcDistribution = this->get_fc_distribution( this->auxFcDistribution );
}

void Extract_FCs_FixedWindowLength::get_draft_of_fc_distribution()
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

	this->nFreq = draftOfFcDistribution(0,1);
	this->auxFcDistribution = new matCUDA::Array<int>( draftOfFcDistribution.getDim(0), draftOfFcDistribution.getDim(1) );
	*(this->auxFcDistribution) = draftOfFcDistribution;
}

matCUDA::Array<int> Extract_FCs_FixedWindowLength::get_fc_distribution( size_t maxDecimationLevel )
{
	// TODO - improve
	size_t nElements = this->firstsLinesToRepeatFromFcDistribution*maxDecimationLevel + this->auxFcDistribution[0].getDim(0) - this->firstsLinesToRepeatFromFcDistribution;
	size_t irow= -1;
	double changeFactor = (double)this->windowLength/(double)this->fcDraftMadeForWinLengthOf;

	Array<int> fcDistribution(nElements,3);
	for( int nDeci = 0; nDeci < maxDecimationLevel; nDeci++ ) {
		for( int iLines = 0; iLines < this->firstsLinesToRepeatFromFcDistribution; iLines++ ) {
			irow++;
			fcDistribution(irow,0) = nDeci + 1;
			fcDistribution(irow,1) = static_cast<int>( this->auxFcDistribution[0](iLines,0)*changeFactor + 0.5 );
			fcDistribution(irow,2) = static_cast<int>( this->auxFcDistribution[0](iLines,1)*changeFactor + 0.5 );
		}
	}

	for( int iLines = this->firstsLinesToRepeatFromFcDistribution; iLines < this->auxFcDistribution[0].getDim(0); iLines++ ) {
		irow++;
		fcDistribution(irow,0) = maxDecimationLevel;
		fcDistribution(irow,1) = static_cast<int>( this->auxFcDistribution[0](iLines,0)*changeFactor + 0.5 );
		fcDistribution(irow,2) = static_cast<int>( this->auxFcDistribution[0](iLines,1)*changeFactor + 0.5 );
	}

	fcDistribution.print();
	//fcDistribution.write2file("teste.txt");

	return fcDistribution;
}

size_t Extract_FCs_FixedWindowLength::determine_decimation_level( size_t nPoints )
{
	size_t nDecimationLevel = 0;
	size_t amountOfFcsInLongerPeriod = this->auxFcDistribution[0]( this->auxFcDistribution[0].getDim(0), 0 ) - this->auxFcDistribution[0]( this->auxFcDistribution[0].getDim(0), 1 ) + 1;

	this->nPointsTS = nPoints;

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