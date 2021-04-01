#include <AMReX_VisMF.H>
#include <Castro.H>

using namespace amrex;

namespace {
    const std::string level_prefix{"Level_"};
}

void Castro::MAESTRO_init ()
{
    BL_PROFILE("Castro::MAESTRO_init()");

    amrex::Print() << "Initialize from MAESTROeX plotfile " << MAESTRO_plotfile << "\n";

    VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

    std::string line, word;
    int step;

    // Header
    {
        std::string File(MAESTRO_plotfile + "/Header");
        Vector<char> fileCharPtr;
        ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
        std::string fileCharPtrString(fileCharPtr.dataPtr());
        std::istringstream is(fileCharPtrString, std::istringstream::in);

        // read in title line
        std::getline(is, line);

        // read in time step number
        is >> start_step;
        GotoNextLine(is);
        ++start_step;

        // read in finest_level
        is >> finest_level;
        GotoNextLine(is);

        // read in step
        is >> step;
        GotoNextLine(is);

        // read in dt
        is >> dt;
        GotoNextLine(is);

        // read in time
        is >> t_old;
        GotoNextLine(is);
        t_new = t_old + dt;

        // read in rel_eps
        is >> rel_eps;
        GotoNextLine(is);

	// define MAESTROeX multifabs
	Vector<MultiFab> sold(finest_level + 1);
	Vector<MultiFab> uold(finest_level + 1);
	
        for (int lev = 0; lev <= finest_level; ++lev) {
            // read in level 'lev' BoxArray from Header
            BoxArray ba;
            ba.readFrom(is);
            GotoNextLine(is);

            // create a distribution mapping
            DistributionMapping dm{ba, ParallelDescriptor::NProcs()};

            // set BoxArray grids and DistributionMapping dmap in AMReX_AmrMesh.H class
            SetBoxArray(lev, ba);
            SetDistributionMap(lev, dm);

            // build MultiFab data
            sold[lev].define(ba, dm, Nscal, ng_s);
            uold[lev].define(ba, dm, AMREX_SPACEDIM, ng_s);
            
            // build FluxRegister data
            if (lev > 0 && reflux_type == 2) {
                flux_reg_s[lev] = std::make_unique<FluxRegister>(
                    ba, dm, refRatio(lev - 1), lev, Nscal);
            }
        }
    }

    // read in the MultiFab data - put it in the "old" MultiFabs
    for (int lev = 0; lev <= finest_level; ++lev) {
        VisMF::Read(sold[lev], amrex::MultiFabFileFullPrefix(lev, MAESTRO_plotfile,
                                                             level_prefix, "snew"));
        VisMF::Read(uold[lev], amrex::MultiFabFileFullPrefix(lev, MAESTRO_plotfile,
                                                             level_prefix, "unew"));
    }

    
    // model file (cell-centered)
    {
        std::string File(MAESTRO_modelfile);
        Vector<char> fileCharPtr;
        ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
        std::string fileCharPtrString(fileCharPtr.dataPtr());
        std::istringstream is(fileCharPtrString, std::istringstream::in);

        // read in cell-centered base state
        for (int i = 0; i < MAESTRO_npts_model; ++i) {
            std::getline(is, line);
            std::istringstream lis(line);
            lis >> word;
            rho0_old.array()(i) = std::stod(word);
            lis >> word;
            p0_old.array()(i) = std::stod(word);
            lis >> word;
            gamma1bar_old.array()(i) = std::stod(word);
            lis >> word;
            rhoh0_old.array()(i) = std::stod(word);
            lis >> word;
            beta0_old.array()(i) = std::stod(word);
            lis >> word;
            psi.array()(i) = std::stod(word);
            lis >> word;
            tempbar.array()(i) = std::stod(word);
            lis >> word;
            etarho_cc.array()(i) = std::stod(word);
            lis >> word;
            tempbar_init.array()(i) = std::stod(word);
        }
    }

}
