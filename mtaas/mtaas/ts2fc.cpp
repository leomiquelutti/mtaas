#include "ts2fc.h"

void SetUpRemoteReferenceConcomitance::find_concomitance_for_stations( std::vector<StationBase> &station )
{
}

void Extract_FCs_FixedWindowLength::get_FCs( std::vector<StationBase> &station )
{
	/*for( int i = 0; i < station.size(); i++ );*/
	Extract_FCs_FixedWindowLength fc( 10 );
}

void Extract_FCs_FixedWindowLength::set_parameters()
{
}

void Extract_FCs_VariableWindowLength::get_FCs( std::vector<StationBase> &station )
{
}