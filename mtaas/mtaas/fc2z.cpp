#include "fc2z.h"

void Evaluate_Z_LS::get_Z( std::vector<StationBase> &station )
{
	for( int istn = 0; istn < station.size(); istn++ ) {
		for( int its = 0; its < station[istn].ts.size(); its++ ) {
			for( int icomb = 0; icomb < station[istn].ts[its].combination.size(); icomb++ ) {
				for( int ifreq = 0; ifreq < station[istn].ts[its].combination[icomb].fr.size(); ifreq++ ) {
					FrequencyResponses auxFr;
					auxFr = station[istn].ts[its].combination[icomb].fr[ifreq];
					//auxFr.transferTensor = new matCUDA::Array<ComplexDouble>( N_CHANNEL_INPUTS, station[istn].ts[its].ch.size() - N_CHANNEL_INPUTS );
					matCUDA::Array<ComplexDouble> auxA = (auxFr.inRR[0].hermitian())*(auxFr.in[0]);
					matCUDA::Array<ComplexDouble> auxY = (auxFr.inRR[0].hermitian())*(auxFr.out[0]);
					matCUDA::Array<ComplexDouble> auxX = auxY.ls( &auxA );
					auxX.print();
					//matCUDA::Array<ComplexDouble> auxHRrHerm = auxFr.inRR[0].hermitian();
					//matCUDA::Array<ComplexDouble> auxHHerm = auxFr.in[0].hermitian();
					///*matCUDA::Array<ComplexDouble>*/ auxX = (auxHRrHerm*auxHHerm).invert()*auxHRrHerm*auxFr.out[0];
					//auxX.print();
					auxFr.transferTensor = new matCUDA::Array<ComplexDouble>( auxX );
					auxFr.param = new Parameters( (auxFr.transferTensor[0]), station[istn].ts[its].nChannels, auxFr.frequency );
					station[istn].ts[its].combination[icomb].fr[ifreq] = auxFr;
					auxX = 0;
					//auxHRrHerm = 0;
					//auxHHerm = 0;
					//delete &auxHRrHerm;
					//delete &auxHHerm;
					//delete( &auxX );
				}
				//station[istn].ts[its].combination[icomb].write_parameters_to_file();
			}
		}
	}
}