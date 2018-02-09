
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include "myfunc.H"
#include "myfunc_F.H"

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    
    main_main();
    
    amrex::Finalize();
    return 0;
}

void main_main ()
{
    //bool runChapter4 = true;
    bool runChapter4 = false;

    // What time is it now?  We'll use this to compute total run time.
    Real strt_time = ParallelDescriptor::second();

    // AMREX_SPACEDIM: number of dimensions
    int n_cell, max_grid_size, nsteps, plot_int;
    Vector<int> is_periodic(AMREX_SPACEDIM,1);  // periodic in all direction by default

    // inputs parameters
    {
        // ParmParse is way of reading inputs from the inputs file
        ParmParse pp;

        // We need to get n_cell from the inputs file - this is the number of cells on each side of 
        //   a square (or cubic) domain.
        pp.get("n_cell",n_cell);

        // The domain is broken into boxes of size max_grid_size
        pp.get("max_grid_size",max_grid_size);

        // Default plot_int to -1, allow us to set it to something else in the inputs file
        //  If plot_int < 0 then no plot files will be writtenq
        plot_int = -1;
        pp.query("plot_int",plot_int);

        // Default nsteps to 0, allow us to set it to something else in the inputs file
        nsteps = 10;
        pp.query("nsteps",nsteps);

        pp.queryarr("is_periodic", is_periodic);

        Array<Real> xr {-1.0, 1.0};
        if (!pp.queryarr("xrange", xr))
        {
          amrex::Print() << "Cannot find xrange in inputs, so the default of {-1.0, 1.0} will be used\n";
        }
    }

    // make BoxArray and Geometry
    BoxArray ba;
    Geometry geom;
    {
        IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
        IntVect dom_hi(AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1));
        Box domain(dom_lo, dom_hi);

        // Initialize the boxarray "ba" from the single box "domain"
        ba.define(domain);

        // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
        ba.maxSize(max_grid_size);

        // This defines the physical box, [-1,1] in each direction.
        RealBox real_box({AMREX_D_DECL(-1.0,-1.0,-1.0)},
                         {AMREX_D_DECL( 1.0, 1.0, 1.0)});

        // This defines a Geometry object
        geom.define(domain,&real_box,CoordSys::cartesian,is_periodic.data());
    }

    // Nghost = number of ghost cells for each array 
    int Nghost = 1;
    
    // Ncomp = number of components for each array
    int Ncomp  = 1;
  
    // How Boxes are distrubuted among MPI processes
    Print() << "create the dm from ba >>> \n";
    DistributionMapping dm(ba);
    Print() << dm;
    Print() << "<<< done creating dm \n";

    //===================================================================================
    // Attempt to loop over all the Boxes in the BoxArray
    int N = ba.size();
    int m1, m2;
    for (int i = 0; i < N; i++)
    {
        Box bx = ba[i];
        IntVect boxSize = bx.size();
        m1 = boxSize[0];
        m2 = boxSize[1];
        Print() << "Box #" << i << " is " << m1 << "-by-" << m2 << " (" << m1*m2  << ")" << "\n";
    }
    //===================================================================================

    // we allocate two phi multifabs; one will store the old state, the other the new.
    MultiFab phi_old(ba, dm, Ncomp, Nghost);
    MultiFab phi_new(ba, dm, Ncomp, Nghost);

    // Initialize phi_new by calling a Fortran routine.
    // MFIter = MultiFab Iterator
    //int countIter = 0;
    for ( MFIter mfi(phi_new); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.validbox();

        //countIter++;
        //amrex::Print() << "\n\n NEW ITERATION" << countIter << endl;

        //amrex::Print() << bx << endl;
        //amrex::Print() << mfi << endl;
        //amrex::Print() << geom.CellSize() << endl;
        //amrex::Print() << geom.ProbLo() << endl;
        //amrex::Print() << geom.ProbHi() << endl;

        init_phi(BL_TO_FORTRAN_BOX(bx),
                 BL_TO_FORTRAN_ANYD(phi_new[mfi]),
                 geom.CellSize(), geom.ProbLo(), geom.ProbHi());
    }

    // compute the time step
    const Real* dx = geom.CellSize();
    Real dt = 0.9*dx[0]*dx[0] / (2.0*AMREX_SPACEDIM);

    // time = starting time in the simulation
    Real time = 0.0;

    // Write a plotfile of the initial data if plot_int > 0 (plot_int was defined in the inputs file)
    if (plot_int > 0)
    {
        int n = 0;
        const std::string& pltfile = amrex::Concatenate("plt",n,5);
        WriteSingleLevelPlotfile(pltfile, phi_new, {"phi"}, geom, time, 0);
    }

    // build the flux multifabs
    std::array<MultiFab, AMREX_SPACEDIM> flux;
    for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
    {
        // flux(dir) has one component, zero ghost cells, and is nodal in direction dir
        BoxArray edge_ba = ba;
        edge_ba.surroundingNodes(dir);
        flux[dir].define(edge_ba, dm, 1, 0);
    }

    // TIME LOOP
    for (int n = 1; n <= nsteps; ++n)
    {
        MultiFab::Copy(phi_old, phi_new, 0, 0, 1, 0);

        // new_phi = old_phi + dt * (something)
        advance(phi_old, phi_new, flux, dt, geom); 
        time = time + dt;
        
        // Tell the I/O Processor to write out which step we're doing
        Print() << "Advanced step " << n << "\n";

        // Write a plotfile of the current data (plot_int was defined in the inputs file)
        if (plot_int > 0 && n%plot_int == 0)
        {
            const std::string& pltfile = amrex::Concatenate("plt",n,5);
            WriteSingleLevelPlotfile(pltfile, phi_new, {"phi"}, geom, time, n);
        }

        //n = nsteps + 1;
    }

    // Call the timer again and compute the maximum difference between the start time and stop time
    //   over all processors
    Real stop_time = ParallelDescriptor::second() - strt_time;
    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    ParallelDescriptor::ReduceRealMax(stop_time,IOProc);

    // Tell the I/O Processor to write out the "run time"
    amrex::Print() << "Run time = " << stop_time << std::endl;



    if (runChapter4) 
    {
    // Tutorial Ch4 : The Basics
    // ========================================================================
    
    // Box
    IntVect lo(AMREX_D_DECL(64, 64, 64));
    IntVect hi(AMREX_D_DECL(127, 127, 127));
    IndexType typ({AMREX_D_DECL(1, 1, 1)});
    Box cc(lo, hi);
    Box nd(lo, hi+1, typ);
    amrex::Print() << "\na cell-centered Box : " << cc << "\n";
    amrex::Print() << "an all nodal Box : " << nd << "\n\n";

    // Box conversions
    Box b0(lo, hi);  // index type: (cell, cell, cell)

    Box b1 = surroundingNodes(b0);  // new box with index type: (node, node, node)
    amrex::Print() << b1 << "\n";
    amrex::Print() << b0 << "\n";

    Box b2 = enclosedCells(b1);  // a new box with index type: (cell, cell, cell)
    if (b2 == b0)
    {
      amrex::Print() << "both b0 and b2 are identical\n";
    }

    Box b3 = convert(b0, {AMREX_D_DECL(0, 1, 0)});  
    // convert b0 tp type: (cell, node, cell)
    amrex::Print() << b3 << "\n";

    amrex::Print() << b3.surroundingNodes() << "\n";
    amrex::Print() << b3.enclosedCells() << "\n\n";

    // Box refinement and coarsening
    Box ccbox({AMREX_D_DECL(16, 16, 16)}, {AMREX_D_DECL(31, 31, 31)});
    ccbox.refine(2);
    amrex::Print() << ccbox << "\n";
    amrex::Print() << ccbox.coarsen(2) << "\n\n";

    // RealBox and Geometry
    int ncell = 64;
    // define a Box with ncell # of cells in each spatial dimension
    IntVect lower(AMREX_D_DECL(0, 0, 0));
    IntVect upper(AMREX_D_DECL(ncell-1, ncell-1, ncell-1));
    Box myDomain(lower, upper);
    // define the RealBox with the proper physical limits of the geometry
    RealBox myRealBox({AMREX_D_DECL(-1.0, -1.0, -1.0)}, {AMREX_D_DECL(1.0, 1.0, 1.0)});
    // choose the Cartesian coordinate system (0)
    int coord = 0;
    // establish periodic BC in each dimension
    std::array<int, AMREX_SPACEDIM> is_periodic2 {AMREX_D_DECL(1, 1, 1)};
    // define a new geometry object using myDomain
    Geometry myGeom(myDomain, &myRealBox, coord, is_periodic2.data());
    // the Geometry object can return various informaiton about the physical domain 
    // and the indexing domain
    const Real* problo = myGeom.ProbLo();   // lower corner of pyhical domain
    Real yhi = myGeom.ProbHi(1);            // y-direction of upper corner
    const Real* mesh = myGeom.CellSize();     // cell size for each dimension
    const Box& domain2 = myGeom.Domain();   // Index domain

    amrex::Print() << *problo << "\n";
    amrex::Print() << yhi << "\n";
    amrex::Print() << *mesh << "\n";
    amrex::Print() << domain2 << "\n\n";
    
    // BoxArray: class used to store collection of Box objects on single AMR level
    Box b({AMREX_D_DECL(0, 0, 0)}, {AMREX_D_DECL(127, 127, 127)});
    BoxArray ba2(b);
    amrex::Print() << "BoxArray size is " << ba2.size() << "\n";
    ba2.maxSize(64);   // chop into boxes of 64^3 cells
    amrex::Print() << "BoxArray size is now " << ba2.size() << "\n";
    amrex::Print() << "BoxArray contents:\n" << ba2 << "\n";

    // Section 4.12 about BaseFab, etc
    Box bx({AMREX_D_DECL(-4, 8, 32)}, {AMREX_D_DECL(32, 64, 48)});
    int numComponents = 4;
    BaseFab<Real> fabTest(bx, numComponents);

    // Using arithmetic functions in BaseFab
    int ncomp = 2;
    FArrayBox fab1(bx, ncomp);
    FArrayBox fab2(bx, ncomp);
    fab1.setVal(1.0);   // fill fab1 with 1.0
    fab1.mult(10.0, 0); // multiply component #0 by 10.0
    fab2.setVal(2.0);   // fill fab2 with 2.0
    Real a = 3.0;
    fab2.saxpy(a, fab1);    // defines fab2 <- a * fab1 + fab2
    amrex::Print() << "\ncount in fab1 = " << fab1.nComp() << "\n";
    amrex::Print() << "count in fab2 = " << fab2.nComp() << "\n";

    // practice with MFIter without tiling...
    int ngrow = 1;
    DistributionMapping dm2(ba2);
    MultiFab mf2(ba2, dm2, ncomp, ngrow);
    for (MFIter mfi(mf2); mfi.isValid(); ++mfi)  // loop over grids
    {        
        const Box& box = mfi.validbox();
        FArrayBox& fab = mf2[mfi];
        Real* a = fab.dataPtr();
        const Box& abox = fab.box();
        // call a user function on the active data
    }


    }


}

