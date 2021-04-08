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
void MaestroData::read_params() {
        
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
		<< "Setting nspec to default: " << NumSpec << "\n";
	maestrodata::maestro_nspec = NumSpec;
    }
    
}

//
// Initialize Maestro data
// 
void MaestroData::setup() {
    
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
    p0_mf.resize(finest_level + 1);
    temp_mf.resize(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
	if (maestrodata::maestro_init_type == 1) {
	    p0_mf[lev] = pltfile->get(lev, "p0");
	} else {
	    p0_mf[lev] = pltfile->get(lev, "p0pluspi");
	}
	temp_mf[lev] = pltfile->get(lev, "tfromp");
    }

    // Note: state stores (rho, rhoh, X_k)
    state_mf.resize(finest_level + 1);
    u_mf.resize(finest_level + 1);
    w0_mf.resize(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
	state_mf[lev].define(grid[lev], dmap[lev], 2 + NumSpec, 0);
        u_mf[lev].define(grid[lev], dmap[lev], AMREX_SPACEDIM, 0);
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
	    if (r < std::end(maestro_var_names)) {
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
            MultiFab::Copy(u_mf[lev], pltfile->get(lev, x), 0, i, 1, 0);
            MultiFab::Copy(w0_mf[lev], pltfile->get(lev, w), 0, i, 1, 0);
        }
    }


    // model file (cell-centered)
    std::string line, word;
    int npts_model = maestrodata::maestro_npts_model;
    r_cc_loc.resize(npts_model);
    rho0.resize(npts_model);
    rhoh0.resize(npts_model);
    p0.resize(npts_model);
    gamma1bar.resize(npts_model);
    tempbar.resize(npts_model);
    
    {
        std::string File(maestrodata::maestro_modelfile);
        Vector<char> fileCharPtr;
        ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
        std::string fileCharPtrString(fileCharPtr.dataPtr());
        std::istringstream is(fileCharPtrString, std::istringstream::in);

        // read in cell-centered base state
        for (int i = 0; i < npts_model; ++i) {
            std::getline(is, line);
            std::istringstream lis(line);
            lis >> word;
	    r_cc_loc[i] = std::stod(word);
	    lis >> word;
            rho0[i] = std::stod(word);
            lis >> word;
            rhoh0[i] = std::stod(word);
            lis >> word;
            p0[i] = std::stod(word);
            lis >> word;
            gamma1bar[i] = std::stod(word);
            lis >> word;
            tempbar[i] = std::stod(word);
        }
    }
    
}

//
// Initialize Castro data using Maestro data
//
void MaestroData::init() {
    
    
}


//
// Test case: read in Maestro data and output on Castro grid
//
void MaestroData::test() {
    
    
}

