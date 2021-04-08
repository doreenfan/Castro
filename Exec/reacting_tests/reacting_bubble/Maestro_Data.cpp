#include <regex>

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

    // setup geometry
    const auto domainBoxFine = pltfile->probDomain(finest_level);
    const auto dxFine = pltfile->cellSize(finest_level);

    // read input MultiFabs
    Vector<MultiFab> p0_mf(finest_level + 1);
    Vector<MultiFab> temp_mf(finest_level + 1);
    
    for (int lev = 0; lev <= finest_level; ++lev) {
	if (maestrodata::maestro_init_type == 1) {
	    p0_mf[lev] = pltfile->get(lev, "p0");
	} else {
	    p0_mf[lev] = pltfile->get(lev, "p0pluspi");
	}
	temp_mf[lev] = pltfile->get(lev, "tfromp");
    }

    // Note: state stores (rho, rhoh, X_k)
    Vector<MultiFab> state_mf(finest_level + 1);  // includes rho, X, rhoh
    Vector<MultiFab> vel_mf(finest_level + 1);
    Vector<MultiFab> w0_mf(finest_level + 1);
    
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

    
    ///
    /// Regrid maestro data onto Castro grid
    /// Initializes multifabs: state, p0, temp, vel, w0 
    ///
    regrid(state_mf, p0_mf, temp_mf, vel_mf, w0_mf);
    
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

void MaestroData::regrid(amrex::Vector<amrex::MultiFab>& state_mf,
			 amrex::Vector<amrex::MultiFab>& p0_mf,
			 amrex::Vector<amrex::MultiFab>& temp_mf,
			 amrex::Vector<amrex::MultiFab>& vel_mf,
			 amrex::Vector<amrex::MultiFab>& w0_mf)
{
    BL_PROFILE("MaestroData::regrid()");

}

//
// Initialize Castro data using Maestro data
//
void MaestroData::init(MultiFab& state)
{
    BL_PROFILE("MaestroData::init()");
    
    
}


//
// Test case: read in Maestro data and output on Castro grid
//
void MaestroData::test(MultiFab& state) {
    
    
}

