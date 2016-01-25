#include "parameters.h"

Parameters::Parameters( matCUDA::Array<ComplexDouble> transferTensor, size_t nChannels, double frequency ) :
	impedanceTensor( 2, 2 ), 
	apparentResistivity( 2, 2 ),
	phase( 2, 2 ),
	tipper( 1, 2 )
{
	switch (nChannels)
	{
	case 4:
		this->initialize_Z_only( transferTensor );
	case 5:
		this->initialize_Z_and_tipper( transferTensor );
	}

	this->frequency = frequency;

	this->evaluate_rhoapp();
	this->evaluate_phase();

	//this->apparentResistivity.print();
	//this->phase.print();
}

void Parameters::evaluate_rhoapp()
{
	// apparent resistivity
	matCUDA::Array<ComplexDouble> z2 = this->impedanceTensor.abs2();

	for( int irow = 0; irow < this->impedanceTensor.getDim(0); irow++ )
		for( int icol = 0; icol < this->impedanceTensor.getDim(1); icol++ )
			this->apparentResistivity( irow, icol ) = z2( irow, icol ).real()/2/PI/this->frequency;
}


void Parameters::evaluate_phase()
{
	// phase
	matCUDA::Array<double> ratioImRe( this->impedanceTensor.getDim(0), this->impedanceTensor.getDim(1) );

	for( int irow = 0; irow < this->impedanceTensor.getDim(0); irow++ )
		for( int icol = 0; icol < this->impedanceTensor.getDim(1); icol++ )
			ratioImRe( irow, icol ) = (this->impedanceTensor( irow, icol )).imag()/(this->impedanceTensor( irow, icol )).real();

	this->phase = ratioImRe.atand();
}
void Parameters::initialize_Z_only( matCUDA::Array<ComplexDouble> transferTensor )
{
	// impedance tensor
	this->impedanceTensor = transferTensor.transpose();
}

void Parameters::initialize_Z_and_tipper( matCUDA::Array<ComplexDouble> transferTensor )
{
	// impedance tensor
	this->impedanceTensor( 0, 0 ) = transferTensor( 0, 1 ); //Zxx
	this->impedanceTensor( 0, 1 ) = transferTensor( 1, 1 ); //Zxy
	this->impedanceTensor( 1, 0 ) = transferTensor( 0, 2 ); //Zyx
	this->impedanceTensor( 1, 1 ) = transferTensor( 1, 2 ); //Zyy

	// tipper
	this->tipper( 0, 0 ) = transferTensor( 0, 0 ); // Tx
	this->tipper( 0, 1 ) = transferTensor( 1, 0 ); // Ty
}