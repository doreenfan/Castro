///
/// Returns true if flag corresponding to source type ``src`` is set.
///
/// @param src      integer, index corresponding to source type
///
    bool source_flag(int src);


///
/// Returns whether any sources are actually applied.
///
    bool apply_sources();

    static int get_output_at_completion();


///
/// Construct source terms at old time
///
/// @param source   MultiFab to save sources to
/// @param state_old   Old state
/// @param state_new   New state
/// @param time     the current simulation time
/// @param dt       the timestep to advance (e.g., go from time to
///                    time + dt)
/// @param apply_to_state   whether to apply the source to state_new
///
    void do_old_sources(
#ifdef MHD
                        amrex::MultiFab& Bx,
                        amrex::MultiFab& By,
                        amrex::MultiFab& Bz,
#endif
                        amrex::MultiFab& source,
                        amrex::MultiFab& state_old, amrex::MultiFab& state_new,
                        amrex::Real time, amrex::Real dt, bool apply_to_state);


///
/// Construct the old-time source
///
/// @param src      integer corresponding to source type to construct
/// @param source   MultiFab to save sources to
/// @param state    State data
/// @param time     the current simulation time
/// @param dt       the timestep to advance (e.g., go from time to
///                    time + dt)
///
    void construct_old_source(int src, amrex::MultiFab& source,
                              amrex::MultiFab& state,
                              amrex::Real time, amrex::Real dt);

///
/// Construct new time sources
///
/// @param source   MultiFab to save sources to
/// @param state_old    Old state
/// @param state_new    New state
/// @param time     the current simulation time
/// @param dt       the timestep to advance (e.g., go from time to
///                    time + dt)
/// @param apply_to_state   whether to apply the source to state_new
///
    void do_new_sources(
#ifdef MHD
                        amrex::MultiFab& Bx,
                        amrex::MultiFab& By,
                        amrex::MultiFab& Bz,
#endif
                        amrex::MultiFab& source,
                        amrex::MultiFab& state_old, amrex::MultiFab& state_new,
                        amrex::Real time, amrex::Real dt, bool apply_to_state);


///
/// Construct the new-time sources.
///
/// @param src          Integer corresponding to type of source to construct
/// @param source       MultiFab to save source to
/// @param state_old    Old state
/// @param state_new    New state
/// @param time         the current simulation time
/// @param dt       t   he timestep to advance (e.g., go from time to
///                    time + dt)
///
    void construct_new_source(int src, amrex::MultiFab& source,
                              amrex::MultiFab& state_old, amrex::MultiFab& state_new,
                              amrex::Real time, amrex::Real dt);


///
/// Evaluate diagnostics quantities describing the effect of an
/// update on the state. The optional parameter local determines
/// whether we want to do this calculation globally over all processes
/// or locally just on this processor. The latter is useful if you
/// are evaluating the contribution from multiple source changes at once
/// and want to amortize the cost of the parallel reductions.
/// Note that the resultant output is volume-weighted.
///
/// @param update   MultiFab containing source
/// @param dt       timestep
/// @param local    boolean, do we want to do this calculation
///                 globally over all processes
///                 or locally just on this processor?
///
    amrex::Vector<amrex::Real> evaluate_source_change(const amrex::MultiFab& update, amrex::Real dt,
                                                      bool local = false);


///
/// Print the change due to a given source term update.
/// We assume here that the input array is lined up with
/// the ``NUM_STATE`` components of ``State_Type`` because we are
/// interested in printing changes to energy, mass, etc.
///
/// @param update       Real Vector of changes
///
    void print_source_change(amrex::Vector<amrex::Real> update);


///
/// For the old-time or new-time sources update, evaluate the change in the state
/// for all source terms, then print the results.
///
/// @param dt       timestep
/// @param is_new   boolean, do we get the new or old state data?
///
    void print_all_source_changes(amrex::Real dt, bool is_new);


///
/// Evaluate the change to the state due to sources and then print it.
///
/// @param source       update to the state
/// @param dt           timestep (will multiply the update)
/// @param source_name  string to associate with this source
///
    void evaluate_and_print_source_change(const amrex::MultiFab& source, amrex::Real dt,
                                          std::string source_name);



///
/// Center source terms in time, storing in ``S_new``
///
/// @param S_new    MultiFab to store centered sources in
/// @param src_old  Old source term
/// @param src_new  New source term
/// @param dt       timestep
///
    void time_center_source_terms (amrex::MultiFab& S_new,
                                   amrex::MultiFab& src_old,
                                   amrex::MultiFab& src_new,
                                   amrex::Real dt = 1.0);

#ifdef SPONGE

///
/// Construct sponge source terms at old timestep
///
/// @param source   MultiFab to save source terms to
/// @param state    Old state
/// @param time     current time
/// @param dt       timestep
///
    void construct_old_sponge_source(amrex::MultiFab& source, amrex::MultiFab& state,
                                     amrex::Real time, amrex::Real dt);


///
/// Construct sponge source terms at new timestep
///
/// @param source       MultiFab to save source terms to
/// @param state_old    Old state
/// @param state_new    New state
/// @param time         current time
/// @param dt           timestep
///
    void construct_new_sponge_source(amrex::MultiFab& source, amrex::MultiFab& state_old,
                                     amrex::MultiFab& state_new, amrex::Real time, amrex::Real dt);


///
/// Compute the sponge source given a state
///
/// @param bx           Box to operate over
/// @param state        input state
/// @param source       source incremented with the sponge
/// @param dt           timestep
///
    void apply_sponge(const amrex::Box& bx,
                      amrex::Array4<amrex::Real const> const state,
                      amrex::Array4<amrex::Real> const source,
                      amrex::Real dt);

#endif


///
/// Create external source terms at old timestep
///
/// @param source   MultiFab to save source terms to
/// @param state    Old state
/// @param time     current time
/// @param dt       timestep
///
    void construct_old_ext_source(amrex::MultiFab& source, amrex::MultiFab& state,
                                  amrex::Real time, amrex::Real dt);


///
/// Create external source terms at new timestep
///
/// @param source       MultiFab to save source terms to
/// @param state_old    Old state
/// @param state_new    New state
/// @param time         current time
/// @param dt           timestep
///
    void construct_new_ext_source(amrex::MultiFab& source,
                                  amrex::MultiFab& state_old, amrex::MultiFab& state_new,
                                  amrex::Real time, amrex::Real dt);


///
/// Fill ``ext_src`` with external sources
///
/// @param time     current time
/// @param dt       timestep
/// @param S_old    Old state
/// @param S_new    New state
/// @param ext_src MultiFab to save sources to
///
    void fill_ext_source(const amrex::Real time, const amrex::Real dt,
                         const amrex::MultiFab& S_old, const amrex::MultiFab& S_new,
                         amrex::MultiFab& ext_src);


///
/// Construct thermal sources at old timestep
///
/// @param source   MultiFab to save sources to
/// @param state    Old state
/// @param time     current time
/// @param dt       timestep
///
    void construct_old_thermo_source(amrex::MultiFab& source, amrex::MultiFab& state,
                                     amrex::Real time, amrex::Real dt);


///
/// Construct thermal sources at new timestep
///
/// @param source       MultiFab to save sources to
/// @param state_old    Old state
/// @param state_new    New state
/// @param time         current time
/// @param dt           timestep
///
    void construct_new_thermo_source(amrex::MultiFab& source, amrex::MultiFab& state_old,
                                     amrex::MultiFab& state_new, amrex::Real time, amrex::Real dt);


///
/// Fill ``ext_src`` with thermal sources
///
/// @param S_old    Old state
/// @param S_new    New state
/// @param ext_src  MultiFab to fill with sources
///
    void fill_thermo_source(amrex::MultiFab& S_new, amrex::MultiFab& ext_src);


///
/// Construct axisymmetric geometry source at old timestep
///
/// @param source   MultiFab to save sources to
/// @param state    Old state
/// @param time     current time
/// @param dt       timestep
///
    void construct_old_geom_source(amrex::MultiFab& source, amrex::MultiFab& state,
                                   amrex::Real time, amrex::Real dt);


///
/// Construct axisymmetric geometry source at new timestep
///
/// @param source       MultiFab to save sources to
/// @param state_old    Old state
/// @param state_new    New state
/// @param time         current time
/// @param dt           timestep
///
    void construct_new_geom_source(amrex::MultiFab& source, amrex::MultiFab& state_old,
                                   amrex::MultiFab& state_new, amrex::Real time, amrex::Real dt);


///
/// Fill ``ext_src`` with axisymmetric geometry sources
///
/// @param time     current time
/// @param dt       timestep
/// @param S        state
/// @param ext_src  MultiFab to fill with sources
///
    void fill_geom_source(amrex::Real time, amrex::Real dt,
                          amrex::MultiFab& S, amrex::MultiFab& ext_src);


#ifdef REACTIONS
///
/// Obtain the effective source term due to reactions on the primitive variables
///
/// @param source   MultiFab to save source term to
/// @param dt       timestep
///
    void get_react_source_prim(amrex::MultiFab& source, amrex::Real time, amrex::Real dt);
#endif

