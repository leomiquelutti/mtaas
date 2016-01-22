#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "matCUDA.h"

// class that contains MT and other useful parameters
class Parameters
{
public:

	Parameters( matCUDA::Array<ComplexDouble> transferTensor, size_t nChannels, double frequency );
	~Parameters();

private:

	matCUDA::Array<double>			apparentResistivity;
	double							frequency;
	matCUDA::Array<ComplexDouble>	impedanceTensor;
	matCUDA::Array<double>			phase;
	matCUDA::Array<ComplexDouble>	tipper;

	void evaluate_rhoapp();
	void evaluate_phase();

	void initialize_Z_only( matCUDA::Array<ComplexDouble> transferTensor );
	void initialize_Z_and_tipper( matCUDA::Array<ComplexDouble> transferTensor );
};

#endif