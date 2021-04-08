#ifndef MAESTRO_DATA_H
#define MAESTRO_DATA_H

#include <AMReX_AmrCore.H>
#include <AMReX_PlotFileUtil.H>

// This vector can be accessed on the GPU
using RealVector = amrex::Gpu::ManagedVector<amrex::Real>;

///
/// Maestro data
/// 
namespace maestrodata {
    extern std::string maestro_plotfile;
    extern std::string maestro_modelfile;
    extern AMREX_GPU_MANAGED int maestro_npts_model;
    extern std::string maestro_first_species;
    extern AMREX_GPU_MANAGED int maestro_nspec;
    extern AMREX_GPU_MANAGED amrex::Real maestro_cutoff_density;
    extern AMREX_GPU_MANAGED int maestro_init_type;
    extern AMREX_GPU_MANAGED bool maestro_spherical;
}

///
/// @class MaestroData
///
class MaestroData {
   public:
    /* public member functions */

    /// constructor
    MaestroData();

    /// destructor
    ~MaestroData();

    /// read maestro-related parameters from parameter file
    void read_params();
    
    /// read in Maestro data from plot file
    void setup();

    /// initialize Castro data using Maestro data
    void init();

    ///
    /// write Castro plotfile for test case
    /// 
    void test();

    
    /* private data members */

    // AMReX plotfile
    amrex::PlotFileData* pltfile;
    
    // geometry data
    int finest_level = 0;
    amrex::Vector<amrex::BoxArray> grid;
    amrex::Vector<amrex::DistributionMapping> dmap;
    amrex::Vector<amrex::Geometry> geom;

    // define Maestro multifabs
    amrex::Vector<amrex::MultiFab> state_mf;  // includes rho, X, rhoh
    amrex::Vector<amrex::MultiFab> p0_mf;
    amrex::Vector<amrex::MultiFab> temp_mf;
    amrex::Vector<amrex::MultiFab> u_mf;
    amrex::Vector<amrex::MultiFab> w0_mf;

    // define Maestro radial base states
    RealVector r_cc_loc;
    RealVector rho0;
    RealVector rhoh0;
    RealVector p0;
    RealVector gamma1bar;
    RealVector tempbar;    
};

#endif
