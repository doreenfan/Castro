
#include <Maestro_Data.H>
#include <Castro.H>

using namespace amrex;

void Castro::MAESTRO_init ()
{
    BL_PROFILE("Castro::MAESTRO_init()");

    if (AMREX_SPACEDIM > 2) {
	amrex::Abort("ERROR: Restart from MAESTROeX only works for DIM=2.");
    }

    MultiFab& S_new = get_new_data(State_Type);
    
    amrex::Print() << "Initialize from MAESTROeX plotfile " << maestrodata::maestro_plotfile << "\n";

    MaestroData maestro_data;
    
    maestro_data.read_params();

    maestro_data.setup();

    // maestro_data.regrid(level, S_new, geom);
    
    // maestro_data.init(S_new);

    // DEBUG 
    maestro_data.test(level, S_new, geom);
    
}
