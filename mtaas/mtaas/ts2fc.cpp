#include "ts2fc.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cmath>

// constructor
Extract_FCs_FixedWindowLength::Extract_FCs_FixedWindowLength( int winLngth ) :
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
}

void Extract_FCs_FixedWindowLength::get_all_FCs( std::vector<StationBase> &station )
{
	//size_t auxN = 0;
	//bool auxAcqstnCntns = true;
		
	Extract_FCs_FixedWindowLength fc;

	// allocate memory for ch[ich].fc
	fc.allocate_memory_for_ch_fc( station );

	// allocate memory for in, inRR and out for each FrequencyResponse
	fc.allocate_memory_for_FrequencyResponses( station );

	//// checkouts
	//for( int istn = 0; istn < station.size(); istn++ )
	//	for( int its = 0; its < station[istn].ts.size(); its++ )
	//		station[istn].ts[its].combination[0].fcDistribution[0].print();


	// set corrections to be used
	fc.set_corrections( station );

	// DO IT - retrieve FCs
	for( int istn = 0; istn < station.size(); istn++ ) {
		for( int its = 0; its < station[istn].ts.size(); its++ ) {

			for( int icomb = 0; icomb < station[istn].ts[its].combination.size(); icomb++ )
				fc.extract_fcs_for_each_combination( station, istn, its, icomb );
		}
	}

	// free memory
	fc.deallocate_memory_for_ch_fc( station );
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
				
				// checkouts				
				station[istn].ts[its].ch[0].systemResponseFreqs[0].print();
				cout << endl;	
				auxFc.print();
				cout << endl;
				station[istn].ts[its].combination[icomb].fcDistribution[0].print();
				cout << endl << idxFreqsToRemove.size() << endl << endl;
				for( int i = 0; i < idxFreqsToRemove.size(); i++ )
					cout << idxFreqsToRemove[i] << " ";
				cout << endl << endl;

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
						nFCsPerSegment = station[istn].ts[its].combination[icomb].fr[ifreq].upperLimitBandWidth - station[istn].ts[its].combination[icomb].fr[ifreq].lowerLimitBandWidth + 1;
						station[istn].ts[its].combination[icomb].fr[ifreq].elementCounter += nSeg*nFCsPerSegment;
					}

					nElements = station[istn].ts[its].combination[icomb].fr[ifreq].elementCounter;
					cout << nElements << endl << endl;
					station[istn].ts[its].combination[icomb].fr[ifreq].in = new matCUDA::Array<ComplexDouble>( nElements, N_CHANNEL_INPUTS );
					station[istn].ts[its].combination[icomb].fr[ifreq].out = new matCUDA::Array<ComplexDouble>( nElements, station[istn].ts[its].nChannels -  N_CHANNEL_INPUTS );
					station[istn].ts[its].combination[icomb].fr[ifreq].inRR = new matCUDA::Array<ComplexDouble>( nElements, N_CHANNEL_INPUTS*station[istn].ts[its].combination[icomb].numberOfConcomitantTs );
				}
			}
			station[istn].ts[its].maxDecimationLevel = *std::max_element( station[istn].ts[its].combination[0].deciLevelEachBlock.begin(), station[istn].ts[its].combination[ station[istn].ts[its].combination.size() - 1 ].deciLevelEachBlock.end() );
		}
	}
}

void Extract_FCs_FixedWindowLength::allocate_memory_for_ch_fc( std::vector<StationBase> &station )
{
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
				delete( station[istn].ts[its].ch[ich].fc );
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

	// sweep through the "blocks" of data - one block is a piecewise continuous set of data
	for( int iblock = 0; iblock < auxComb.idxBgn[0].getDim(0); iblock++ ) {

		// sweep through the decimation levels of each block
		for( int ideci = 0; ideci < auxComb.deciLevelEachBlock[iblock]; ideci++ ) {

			nSegments = floor( ( auxComb.idxEnd[0](iblock,0) - auxComb.idxBgn[0](iblock,0) )/( this->windowLength*( 1 - this->overlap ) ) );
			cout << nSegments << endl;

			// sweep through the segments of each decimation level of each block
			for( int iseg = 0; iseg < nSegments; iseg++ ) {

				for( int iInvolvedTs = 0; iInvolvedTs < auxComb.numberOfConcomitantTs + 1; iInvolvedTs++ ) {
					
					istart = auxComb.idxBgn[0](iblock,iInvolvedTs);
					iend = istart + this->windowLength - 1;

					size_t iAuxStn = 0, iAuxTs = 0, upperLimit = 0;
					// transform segments of "this" TS and its concomitants combinations in FC
					while( iend < auxComb.idxEnd[0](iblock,iInvolvedTs) ) {

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

						// get FCs for station of interest (this one)
						for( int ich = 0; upperLimit; ich++ ) {
							auxReal = station[iAuxStn].ts[iAuxTs].ch[ich].timeSeries[0].submatrix( istart, iend, 0, 0 );
							auxReal = auxReal.detrend();
							auxReal = auxReal.elementWiseMultiply( this->auxDpss );
							auxCmplx = fft( &auxReal );
							this->correct_fcs( &auxCmplx, &station[iAuxStn].ts[iAuxTs].ch[ich].correction[0] );
							station[iAuxStn].ts[iAuxTs].ch[ich].fc = &auxCmplx;
						}

						//// TEST IT NOW!!!!
						//// get FCs for station of interest (this one)
						//if( iInvolvedTs == 0 ) {
						//	for( int ich = 0; ich < station[istn].ts[its].ch.size(); ich++ ) {
						//		auxReal = station[istn].ts[its].ch[ich].timeSeries[0].submatrix( istart, iend, 0, 0 );
						//		auxReal = auxReal.detrend();
						//		auxReal = auxReal.elementWiseMultiply( this->auxDpss );
						//		auxCmplx = fft( &auxReal );
						//		this->correct_fcs( &auxCmplx, &station[istn].ts[its].ch[ich].correction[0] );
						//		station[istn].ts[its].ch[ich].fc = &auxCmplx;
						//	}
						//}
						//// get FCs for remote stations
						//else {
						//	size_t auxCombStn = auxComb.idxStn[iInvolvedTs - 1];
						//	size_t auxCombTs = auxComb.idxTs[iInvolvedTs - 1];
						//	for( int ich = 0; ich < N_CHANNEL_INPUTS; ich++ ) { // TODO - IMPROVE this number"2"
						//		auxReal = station[auxCombStn].ts[auxCombTs].ch[ich].timeSeries[0].submatrix( istart, iend, 0, 0 );
						//		auxReal = auxReal.detrend();
						//		auxReal = auxReal.elementWiseMultiply( this->auxDpss );
						//		auxCmplx = fft( &auxReal );
						//		this->correct_fcs( &auxCmplx, &station[auxCombStn].ts[auxCombTs].ch[ich].correction[0] );
						//		station[auxCombStn].ts[auxCombTs].ch[ich].fc = &auxCmplx;
						//	}
						//}

						istart += winPace;
						iend = istart + this->windowLength - 1;
					}
				}

				// get in

				// get inRR

				// get out
			}
		}
	}
}

void Extract_FCs_FixedWindowLength::correct_fcs( matCUDA::Array<ComplexDouble> *data, matCUDA::Array<ComplexDouble> *correction )
{
	for( int i = 0; i < this->nFreq; i++ )
		data[0]( i ) = data[0]( i + 1 )*correction[0]( i );
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
			cout << station[istn].ts[its].maxDecimationLevel << endl;
			std::vector<double> dr( station[istn].ts[its].maxDecimationLevel );
			for( int i = 0;  i < dr.size(); i++ )
				dr[i] = pow( this->factorOfEachDecimationLevel, i )/station[istn].ts[its].samplingFrequency;

			// correct for first difference if appropriate
			// correct unist of fc's
			for( int ich = 0; ich < station[istn].ts[its].ch.size(); ich++ ) {
				station[istn].ts[its].ch[ich].correction = new Array<ComplexDouble>( station[istn].ts[its].maxDecimationLevel, npts2 );
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

	this->auxDpss = new matCUDA::Array<double>( matCUDA::dpss<double>( this->windowLength, 1, 0 ) );
	
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