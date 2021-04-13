#include <regex>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <AMReX_VisMF.H>
#include <AMReX_ParmParse.H>
#include <Maestro_Data.H>
#include <Castro.H>

using namespace amrex;

/// Maestro parameter data
std::string maestrodata::maestro_plotfile;
std::string maestrodata::maestro_modelfile;
AMREX_GPU_MANAGED int maestrodata::maestro_npts_model;
std::string maestrodata::maestro_first_species;
AMREX_GPU_MANAGED int maestrodata::maestro_nspec;
AMREX_GPU_MANAGED amrex::Real maestrodata::maestro_cutoff_density;
AMREX_GPU_MANAGED int maestrodata::maestro_init_type;
AMREX_GPU_MANAGED bool maestrodata::maestro_spherical;

/// constructor
MaestroData::MaestroData() = default;

/// destructor
MaestroData::~MaestroData() = default;

//
// Read in parameters from input file
//
void MaestroData::read_params()
{        
    // read in maestro-related input parameters
    {
    ParmParse pp("castro");
#include <maestro_queries.H>
    }
    
    // ERRORS: missing necessary parameters
    if (maestrodata::maestro_plotfile == "fillme") {
        Abort("ERROR: MAESTRO_plotfile was not specified!\n");
    }

    if (maestrodata::maestro_first_species == "fillme") {
        Abort("ERROR: MAESTRO_first_species was not specified!\n");
    }

    // WARNINGS: set parameters to default values
    if (maestrodata::maestro_modelfile == "fillme") {
        Print() << "WARNING: MAESTRO_modelfile was not specified!\n"
		<< "Setting modelfile to default: "
		<< maestrodata::maestro_plotfile << "/BaseCC_0 \n";
	maestrodata::maestro_modelfile =
	    maestrodata::maestro_plotfile + "/BaseCC_0";
    }
    
    if (maestrodata::maestro_nspec == 0) {
        Print() << "WARNING: MAESTRO_nspec was not specified!\n"
		<< "Setting nspec to Castro default: " << NumSpec << "\n";
	maestrodata::maestro_nspec = NumSpec;
    }   
}

//
// Initialize Maestro data
// 
void MaestroData::setup()
{
    BL_PROFILE("MaestroData::setup()");
    
    // set input parameters
    pltfile = new amrex::PlotFileData(maestrodata::maestro_plotfile);
    finest_level = pltfile->finestLevel();

    // setup geometry
    const auto probLo = pltfile->probLo();
    const auto probHi = pltfile->probHi();

    grid.resize(finest_level + 1);
    dmap.resize(finest_level + 1);
    geom.resize(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        grid[lev] = pltfile->boxArray(lev);
        dmap[lev] = pltfile->DistributionMap(lev);

        Box domain = pltfile->probDomain(lev);
        RealBox real_box(probLo, probHi);
        geom[lev].define(domain, &real_box);
    }

    // read input MultiFabs
    p0_mf.resize(finest_level + 1);
    temp_mf.resize(finest_level + 1);

    for (int lev = 0; lev <= finest_level; ++lev) {
	p0_mf[lev].define(grid[lev], dmap[lev], 1, 0);
        temp_mf[lev].define(grid[lev], dmap[lev], 1, 0);
    }
    
    // store p0 and add pi if necessary
    for (int lev = 0; lev <= finest_level; ++lev) {
	if (maestrodata::maestro_init_type == 1) {
	    p0_mf[lev] = pltfile->get(lev, "p0");
	} else {
	    p0_mf[lev] = pltfile->get(lev, "p0pluspi");
	}
	temp_mf[lev] = pltfile->get(lev, "tfromp");
    }

    // Note: state stores (rho, rhoh, X_k)
    state_mf.resize(finest_level + 1);  // includes rho, X, rhoh
    vel_mf.resize(finest_level + 1);
    w0_mf.resize(finest_level + 1);
    
    for (int lev = 0; lev <= finest_level; ++lev) {
	state_mf[lev].define(grid[lev], dmap[lev], 2 + NumSpec, 0);
        vel_mf[lev].define(grid[lev], dmap[lev], AMREX_SPACEDIM, 0);
        w0_mf[lev].define(grid[lev], dmap[lev], AMREX_SPACEDIM, 0);
    }

    // full states
    // get list of variables in maestro plotfile
    auto maestro_var_names = pltfile->varNames();
    
    for (int lev = 0; lev <= finest_level; ++lev) {
	// set default values to zero
	state_mf[lev].setVal(0.);
	
	MultiFab::Copy(state_mf[lev], pltfile->get(lev, "rho"), 0, 0, 1, 0);
	MultiFab::Copy(state_mf[lev], pltfile->get(lev, "rhoh"), 0, 1, 1, 0);

	// get all species
	for (int i = 0; i < NumSpec; ++i) {
	    std::string spec_string = "X(";
	    spec_string += short_spec_names_cxx[i];
	    spec_string += ")";

	    // check that the species is in maestro plotfile
	    auto r = std::find(std::begin(maestro_var_names), std::end(maestro_var_names), spec_string);
	    if (r == std::end(maestro_var_names)) {
		Print() << "WARNING: Could not find " << spec_string << " in MAESTROeX plotfile!\n";
	    } else {
		MultiFab::Copy(state_mf[lev], pltfile->get(lev, spec_string), 0, 2+i, 1, 0);
	    }
	}
    }
    
    // velocities
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        std::string x = "vel";
        std::string w = "w0";
        x += (120 + i);
        w += (120 + i);

        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(vel_mf[lev], pltfile->get(lev, x), 0, i, 1, 0);
            MultiFab::Copy(w0_mf[lev], pltfile->get(lev, w), 0, i, 1, 0);
        }
    }

    
    // model file (cell-centered)
    std::string line, word;
    int npts_model = maestrodata::maestro_npts_model;
    r_model.resize(npts_model);
    rho0_model.resize(npts_model);
    rhoh0_model.resize(npts_model);
    p0_model.resize(npts_model);
    gamma1bar_model.resize(npts_model);
    tempbar_model.resize(npts_model);
    
    {
        std::string File(maestrodata::maestro_modelfile);
        Vector<char> fileCharPtr;
        ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
        std::string fileCharPtrString(fileCharPtr.dataPtr());
        std::istringstream is(fileCharPtrString, std::istringstream::in);

	// read in header
	std::getline(is, line);
	
        // read in cell-centered base state
        for (int i = 0; i < npts_model; ++i) {
            std::getline(is, line);
            std::istringstream lis(line);
            lis >> word;
	    r_model[i] = std::stod(word);
	    lis >> word;
            rho0_model[i] = std::stod(word);
            lis >> word;
            rhoh0_model[i] = std::stod(word);
            lis >> word;
            p0_model[i] = std::stod(word);
            lis >> word;
            gamma1bar_model[i] = std::stod(word);
            lis >> word;
            tempbar_model[i] = std::stod(word);
        }
    }
    
}

///
/// Regrid maestro data onto Castro grid
/// Initializes mstate 
///
void MaestroData::regrid(MultiFab& s_in)
{
    BL_PROFILE("MaestroData::regrid()");

    // Castro grid
    const amrex::Geometry& dgeom = DefaultGeometry();
    const amrex::BoxArray& ba = s_in.boxArray();
    const amrex::DistributionMapping& dm = s_in.DistributionMap();

    // define new Maestro data grid
    mstate.define(ba, dm, s_in.nComp(), s_in.nGrow());
    mstate.setVal(0.);

    // put Maestro data onto new grid
    
}

//
// Initialize Castro data using Maestro data
//
void MaestroData::init(MultiFab& s_in)
{
    BL_PROFILE("MaestroData::init()");

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(mstate, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
	const Box& bx = mfi.tilebox();

	initdata(bx, mstate.array(mfi));
    }

    MultiFab::Copy(s_in, mstate, 0, 0, mstate.nComp(), mstate.nGrow());
    
}

void MaestroData::initdata(const Box& bx,
			   Array4<Real> const& state)
{
    const Real minpres = p0_model[maestrodata::maestro_npts_model-1];
    
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    {
	// local variables
	Real ekin, pressure, entropy;

	eos_t eos_state;

	// set pressure 
	// state(i,j,k,UEDEN) = state(i,j,k,UEDEN);

	if (maestrodata::maestro_init_type == 1 || maestrodata::maestro_init_type == 2) {

	    // load pressure from our temporary storage field
	    pressure = amrex::max(state(i,j,k,UEDEN),minpres);

	    // compute e and T
	    eos_state.p   = pressure;
	    eos_state.rho = state(i,j,k,URHO);
	    eos_state.T   = state(i,j,k,UTEMP);
	    for (int n = 0; n < NumSpec; ++n) {
		eos_state.xn[n]  = state(i,j,k,UFS+n);
	    }

	    eos(eos_input_rp, eos_state);

	    // set tempbar
	    // state(i,j,k,UTEMP) = state(i,j,k,UTEMP);
	    state(i,j,k,UEINT) = eos_state.e;

	    // compute kinetic energy
	    ekin = 0.5*state(i,j,k,URHO)*(state(i,j,k,UMX)*state(i,j,k,UMX)
					  + state(i,j,k,UMY)*state(i,j,k,UMY)
#if AMREX_SPACEDIM == 3
					  + state(i,j,k,UMZ)*state(i,j,k,UMZ)
#endif
					  );

	    // convert velocity to momentum
	    state(i,j,k,UMX) = state(i,j,k,UMX)*state(i,j,k,URHO);
	    state(i,j,k,UMY) = state(i,j,k,UMY)*state(i,j,k,URHO);
#if AMREX_SPACEDIM == 3
	    state(i,j,k,UMZ) = state(i,j,k,UMZ)*state(i,j,k,URHO);
#endif

	    // compute rho*e
	    state(i,j,k,UEINT) = state(i,j,k,URHO) * state(i,j,k,UEINT);

	    // compute rho*E = rho*e + ke
	    state(i,j,k,UEDEN) = state(i,j,k,UEINT) + ekin;

	    // convert X to rhoX
	    for (int n = 0; n < NumSpec; ++n) {
		state(i,j,k,UFS+n) = state(i,j,k,URHO) * state(i,j,k,UFS+n);
	    }

	} else if (maestrodata::maestro_init_type == 3) {

	    // load pressure from our temporary storage field
	    pressure = amrex::max(state(i,j,k,UEDEN),minpres);

	    // compute rho and e
	    eos_state.T  = state(i,j,k,UTEMP);
	    eos_state.p  = pressure;
	    for (int n = 0; n < NumSpec; ++n) {
		eos_state.xn[n] = state(i,j,k,UFS+n);
	    }

	    eos(eos_input_tp, eos_state);

	    state(i,j,k,URHO)  = eos_state.rho;
	    state(i,j,k,UEINT) = eos_state.e;

	    // compute kinetic energy
	    ekin = 0.5*state(i,j,k,URHO)*(state(i,j,k,UMX)*state(i,j,k,UMX)
					  + state(i,j,k,UMY)*state(i,j,k,UMY)
#if AMREX_SPACEDIM == 3
					  + state(i,j,k,UMZ)*state(i,j,k,UMZ)
#endif
					  );

	    // convert velocity to momentum
	    state(i,j,k,UMX) = state(i,j,k,UMX)*state(i,j,k,URHO);
	    state(i,j,k,UMY) = state(i,j,k,UMY)*state(i,j,k,URHO);
#if AMREX_SPACEDIM == 3
	    state(i,j,k,UMZ) = state(i,j,k,UMZ)*state(i,j,k,URHO);
#endif

	    // compute rho*e
	    state(i,j,k,UEINT) = state(i,j,k,URHO) * state(i,j,k,UEINT);

	    // compute rho*E = rho*e + ke
	    state(i,j,k,UEDEN) = state(i,j,k,UEINT) + ekin;

	    // convert X to rhoX
	    for (int n = 0; n < NumSpec; ++n) {
		state(i,j,k,UFS+n) = state(i,j,k,URHO) * state(i,j,k,UFS+n);
	    }

	} else if (maestrodata::maestro_init_type == 4) {

	    // load pressure from our temporary storage field
	    pressure = amrex::max(state(i,j,k,UEDEN),minpres);

	    // load entropy from our temporary storage field
	    entropy = state(i,j,k,UEINT);

	    // compute kinetic energy
	    ekin = 0.5*state(i,j,k,URHO)*(state(i,j,k,UMX)*state(i,j,k,UMX)
					  + state(i,j,k,UMY)*state(i,j,k,UMY)
#if AMREX_SPACEDIM == 3
					  + state(i,j,k,UMZ)*state(i,j,k,UMZ)
#endif
					  );

	    // compute rho, T, and e
	    eos_state.p  = pressure;
	    eos_state.s  = entropy;
	    for (int n = 0; n < NumSpec; ++n) {
		eos_state.xn[n] = state(i,j,k,UFS+n);
	    }
	
	    eos(eos_input_ps, eos_state);
	
	    state(i,j,k,URHO)  = eos_state.rho;
	    state(i,j,k,UEINT) = eos_state.e;
	    state(i,j,k,UTEMP) = eos_state.T;

	    // convert velocity to momentum
	    state(i,j,k,UMX) = state(i,j,k,UMX)*state(i,j,k,URHO);
	    state(i,j,k,UMY) = state(i,j,k,UMY)*state(i,j,k,URHO);
#if AMREX_SPACEDIM == 3
	    state(i,j,k,UMZ) = state(i,j,k,UMZ)*state(i,j,k,URHO);
#endif

	    // compute rho*e
	    state(i,j,k,UEINT) = state(i,j,k,URHO) * state(i,j,k,UEINT);

	    // compute rho*E = rho*e + ke
	    state(i,j,k,UEDEN) = state(i,j,k,UEINT) + ekin;

	    // convert X to rhoX
	    for (int n = 0; n < NumSpec; ++n) {
		state(i,j,k,UFS+n) = state(i,j,k,URHO) * state(i,j,k,UFS+n);
	    }
	
	}
    });
}

//
// Test for debugging: read in Maestro data and output on Castro grid
//
void MaestroData::test(MultiFab& s_in) {
    
    // Write out data (level 0) on initial maestro grid
    VisMF::Write(state_mf[0], "maestro_state0");

    regrid(s_in);

    init(s_in);
    
    // Write out final castro state
    VisMF::Write(s_in, "castro_Snew");
}

