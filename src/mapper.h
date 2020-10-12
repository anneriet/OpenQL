/**
 * @file   mapper.h
 * @date   06/2018 - now
 * @author Hans van Someren
 * @brief  openql virtual to real qubit mapping and routing
 */

#ifndef QL_MAPPER_H
#define QL_MAPPER_H

#include <random>
#include <chrono>
#include <ctime>
#include <ratio>
#include "utils.h"
#include "platform.h"
#include "kernel.h"
#include "resource_manager.h"
#include "gate.h"
#include "scheduler.h"
//#include "metrics.h"

// Note on the use of constructors and Init functions for classes of the mapper
// -----------------------------------------------------------------------------
// Almost all classes of the mapper have one or more members that require initialization
// using a value that was passed on to the Mapper.Init function as a parameter (i.e. platform, cycle_time).
// Dealing with those initializations in the nested constructors was cumbersome.
// Hence, the constructors create just skeleton objects which need explicit initialization before use.
// Such initialization is provided by a class local Init function for a virgin object,
// or by copying an existing object to it.
// The constructors are trivial by this and can be synthesized by default.
//
// Construction of skeleton objects requires the used classes to provide such (non-parameterized) constructors.


// =========================================================================================
// the standard assert didn't stop the compiler on its failure (causing long debug sessions)
// so added the next to get around this and enforce assert working
// it can be replaced by the standard assert if it works
static void assert_fail(const char *f, int l, const char *s)
{
    FATAL("assert " << s << " failed in file " << f << " at line " << l);
}
#define MapperAssert(condition)   { if (!(condition)) { assert_fail(__FILE__, __LINE__, #condition); } }


// =========================================================================================
// Grid: definition and access functions to the grid of qubits that supports the real qubits.
// Maintain several maps to ease navigating in the grid; these are constant after initialization.

// Grid
//
// Config file definitions:
//  nq:                 hardware_settings.qubit_number
//  ncores:             hardware_settings.number_of_cores
//  topology.conn;      gc_specified/gc_full: topology.connectivity: how connectivity between qubits is specified
//  topology.form;      gf_xy/gf_irregular: topology.form: how relation between neighbors is specified
//  topology.x_size/y_size: x/y space, defines underlying grid (only gf_xy)
//  topology.qubits:    mapping of qubit to x/y coordinates (defines x[i]/y[i] for each qubit i) (only gf_xy)
//  topology.edges:     mapping of edge (physical connection between 2 qubits) to its src and dst qubits (defines nbs)
//
// Grid public members (apart from nq):
//  form:               how relation between neighbors is specified
//  Distance(qi,qj):    distance in physical connection hops from real qubit qi to real qubit qj;
//                      - computing it relies on nbs (and Floyd-Warshall) (gf_xy and gf_irregular)
//  nbs[qi]:            list of neighbor real qubits of real qubit qi
//                      - nbs can be derived from topology.edges (gf_xy and gf_irregular)
//  Normalize(qi, neighborlist):    rotate neighborlist such that largest angle diff around qi is behind last element
//                      relies on nbs, and x[i]/y[i] (gf_xy only)
//
// For an irregular grid form, only nq and edges (so nbs) need to be specified; distance is computed from nbs:
// - there is no underlying rectangular grid, so there are no defined x and y coordinates of qubits;
//   this means that Normalize as needed by mappathselect==borders cannot work
// Below, we support regular (xy) grids which need not be fully assigned; this requires edges (so nbs) to be defined,
//   from which distance is computed; also we have x/y coordinates per qubit specified in the configuration file
//   An underlying grid with x/y coordinates comes in use for:
//   - crossbars
//   - cclight qwg assignment (done in another manner now)
//   - when mappathselectopt==borders
//
// Not implemented:
// forms gf_cross and gf_plus: given x_size and y_size, the relations are implicitly defined by the internal
//      diagonal (gf_cross) or horizontal/vertical (gf_plus) connections between grid points (qubits);
//      with gf_cross only half the grid is occupied by a qubit; the grid point (0,0) doesn't have a qubit, (1,0) and (0,1) do;
//      topology.qubits and topology.edges need not be present in the configuration file;
//      Distance in both forms would be defined by a formula, not a function
typedef
enum GridConnectivity
{
    gc_specified,   // "specified": edges are specified in "edges" section
    gc_full         // "full": qubits are fully connected by edges
} gridconn_t;

typedef
enum GridForms
{
    gf_xy,          // nodes have explicit neighbor definitions, qubits have explicit x/y coordinates
    gf_irregular    // nodes have explicit neighbor definitions, qubits don't have x/y coordinates
} gridform_t;

class Grid
{
public:
    const ql::quantum_platform* platformp;    // current platform: topology
    size_t nq;                          // number of qubits in the platform
    size_t ncores;                      // number of cores in the platform
                                        // Grid configuration, all constant after initialization
    gridform_t form;                    // form of grid
    gridconn_t conn;                    // connectivity of grid
    int nx;                             // length of x dimension (x coordinates count 0..nx-1)
    int ny;                             // length of y dimension (y coordinates count 0..ny-1)

    typedef std::list<size_t> neighbors_t;  // neighbors is a list of qubits
    std::map<size_t,neighbors_t> nbs;   // nbs[i] is list of neighbor qubits of qubit i
    std::map<size_t,int> x;             // x[i] is x coordinate of qubit i
    std::map<size_t,int> y;             // y[i] is y coordinate of qubit i
    std::vector<std::vector<size_t>>  dist; // dist[i][j] is computed distance between qubits i and j;

// Grid initializer
// initialize mapper internal grid maps from configuration
// this remains constant over multiple kernels on the same platform
void Init(const ql::quantum_platform* p)
{
    DOUT("Grid::Init");
    platformp = p;
    nq = platformp->qubit_number;
    DOUT("... number of real qbits=" << nq);

    std::string formstr;
    if (platformp->topology.count("form") <= 0)
    {
        formstr = "xy";
    }
    else
    {
        formstr = platformp->topology["form"].get<std::string>();
    }
    if (formstr == "xy") { form = gf_xy; }
    if (formstr == "irregular") { form = gf_irregular; }

    if (form == gf_irregular)
    {
        // irregular can do without topology.x_size, topology.y_size, and topology.qubits
        nx = 0;
        ny = 0;
    }
    else
    {
        // gf_xy have an x/y space; coordinates are explicitly specified
        nx = platformp->topology["x_size"];
        ny = platformp->topology["y_size"];
    }
    DOUT("... formstr=" << formstr << "; form=" << form << "; nx=" << nx << "; ny=" << ny);

    InitCores();
    InitXY();
    InitNbs();
    SortNbs();
    ComputeDist();
    DPRINTGrid();
}

// core index from qubit index
// when multi-core assumes full and uniform core connectivity
size_t CoreOf(size_t qi)
{
    if (ncores == 1) return 0;
    MapperAssert(conn == gc_full);
    size_t nqpc = nq/ncores;
    return qi/nqpc;
}

// inter-core hop from qs to qt?
bool IsInterCoreHop(size_t qs, size_t qt)
{
    return CoreOf(qs) != CoreOf(qt);
}

// distance between two qubits
// formulae for convex (hole free) topologies with underlying grid and with bidirectional edges:
//      gf_cross:   std::max( std::abs( x[to_realqi] - x[from_realqi] ), std::abs( y[to_realqi] - y[from_realqi] ))
//      gf_plus:    std::abs( x[to_realqi] - x[from_realqi] ) + std::abs( y[to_realqi] - y[from_realqi] )
// when the neighbor relation is defined (topology.edges in config file), Floyd-Warshall is used, which currently is always
size_t Distance(size_t from_realqi, size_t to_realqi)
{
    return dist[from_realqi][to_realqi];
}

// coredistance between two qubits
// when multi-core assumes full and uniform core connectivity
size_t CoreDistance(size_t from_realqi, size_t to_realqi)
{
    if (CoreOf(from_realqi) == CoreOf(to_realqi)) return 0;
    return 1;
}

// minimum number of hops between two qubits is always >= distance(from, to)
// and inside one core (or without multi-core) the minimum number of hops == distance
//
// however, in multi-core with inter-core hops, an inter-core hop cannot execute a 2qgate
// so when the minimum number of hops are all inter-core hops (so distance(from,to) == coredistance(from,to))
// and no 2qgate has been placed yet, then at least one additional inter-core hop is needed for the 2qgate,
// the number of hops required being at least distance+1;
//
// we assume below that a valid path exists with distance+1 hops;
// this fails when not all qubits in a core support connections to all other cores;
// see the check in InitNbs
size_t MinHops(size_t from_realqi, size_t to_realqi)
{
    size_t d = Distance(from_realqi, to_realqi);
    size_t cd = CoreDistance(from_realqi, to_realqi);
    MapperAssert (cd <= d);
    if (cd == d)
    {
        return d+1;
    }
    else
    {
        return d;
    }
}


// return clockwise angle around (cx,cy) of (x,y) wrt vertical y axis with angle 0 at 12:00, 0<=angle<2*pi
double Angle(int cx, int cy, int x, int y)
{
    const double pi = 4*std::atan(1);
    double a = std::atan2((x-cx),(y-cy));
    if (a < 0) a += 2*pi;
    return a;
}

// rotate neighbors list such that largest angle difference between adjacent elements is behind back;
// this is needed when a given subset of variations from a node is wanted (mappathselect==borders);
// and this can only be computed when there is an underlying x/y grid (so not for form==gf_irregular)
void Normalize( size_t src, neighbors_t& nbl )
{
    if (form != gf_xy)
    {
        // there are no implicit/explicit x/y coordinates defined per qubit, so no sense of nearness
        std::string mappathselectopt = ql::options::get("mappathselect");
        MapperAssert ("borders" != mappathselectopt);
        return;
    }

    // std::cout << "Normalizing list from src=" << src << ": ";
    // for (auto dn : nbl) { std::cout << dn << " "; } std::cout << std::endl;

    const double pi = 4*std::atan(1);
    if (nbl.size() == 1)
    {
        // DOUT("... size was 1; unchanged");
        return;
    }

    // find maxinx index in neighbor list before which largest angle difference occurs
    int maxdiff = 0;                            // current maximum angle difference in loop search below
    neighbors_t::iterator maxinx = nbl.begin(); // before which max diff occurs

    // for all indices in and its next one inx compute angle difference and find largest of these
    for (neighbors_t::iterator in = nbl.begin(); in != nbl.end(); in++)
    {
        double a_in = Angle(x[src], y[src], x[*in], y[*in]);

        neighbors_t::iterator inx = std::next(in); if (inx == nbl.end()) inx = nbl.begin();
        double a_inx = Angle(x[src], y[src], x[*inx], y[*inx]);

        int diff = a_inx - a_in; if (diff < 0) diff += 2*pi;
        if (diff > maxdiff)
        {
            maxdiff = diff;
            maxinx = inx;
        }
    }

    // and now rotate neighbor list so that largest angle difference is behind last one
    neighbors_t   newnbl;
    for (neighbors_t::iterator in = maxinx; in != nbl.end(); in++)
    {
        newnbl.push_back(*in);
    }
    for (neighbors_t::iterator in = nbl.begin(); in != maxinx; in++)
    {
        newnbl.push_back(*in);
    }
    nbl = newnbl;

    // std::cout << "... rotated; result: ";
    // for (auto dn : nbl) { std::cout << dn << " "; } std::cout << std::endl;
}

// Floyd-Warshall dist[i][j] = shortest distances between all nq qubits i and j
void ComputeDist()
{
    // initialize all distances to maximum value, to neighbors to 1, to itself to 0
    dist.resize(nq); for (size_t i=0; i<nq; i++) dist[i].resize(nq, MAX_CYCLE);
    for (size_t i=0; i<nq; i++)
    {
        dist[i][i] = 0;
        for (size_t j: nbs[i])
        {
            dist[i][j] = 1;
        }
    }

    // find shorter distances by gradually including more qubits (k) in path
    for (size_t k=0; k<nq; k++)
    {
        for (size_t i=0; i<nq; i++)
        {
            for (size_t j=0; j<nq; j++)
            {
               if (dist[i][j] > dist[i][k] + dist[k][j])
               {
                   dist[i][j] = dist[i][k] + dist[k][j];
               }
            }
        }
    }
#ifdef debug
    for (size_t i=0; i<nq; i++)
    {
        for (size_t j=0; j<nq; j++)
        {
            if (form == gf_cross)
            {
                MapperAssert (dist[i][j] == (std::max( std::abs( x[i] - x[j] ), std::abs( y[i] - y[j] ))) );
            }
            else if (form == gf_plus)
            {
                MapperAssert (dist[i][j] == (std::abs( x[i] - x[j] ) + std::abs( y[i] - y[j] )) );
            }

        }
    }
#endif
}

void DPRINTGrid()
{
    if ( ql::utils::logger::LOG_LEVEL >= ql::utils::logger::log_level_t::LOG_DEBUG )
        PrintGrid();
}
void PrintGrid()
{
    if (form != gf_irregular)
    {
	    for (size_t i=0; i<nq; i++)
	    {
	        std::cout << "qubit[" << i << "]=(" << x[i] << "," << y[i] << ")";
	        std::cout << " has neighbors ";
	        for (auto & n : nbs[i])
	        {
	            std::cout << "qubit[" << n << "]=(" << x[n] << "," << y[n] << ") ";
	        }
	        std::cout << std::endl;
	    }
    }
    else
    {
	    for (size_t i=0; i<nq; i++)
	    {
	        std::cout << "qubit[" << i << "]";
	        std::cout << " has neighbors ";
	        for (auto & n : nbs[i])
	        {
	            std::cout << "qubit[" << n << "] ";
	        }
	        std::cout << std::endl;
	    }
    }
    for (size_t i=0; i<nq; i++)
    {
        std::cout << "qubit[" << i << "] distance(" << i << ",j)=";
        for (size_t j=0; j<nq; j++)
        {
            std::cout << Distance(i,j) << " ";
        }
        std::cout << std::endl;
    }
    for (size_t i=0; i<nq; i++)
    {
        std::cout << "qubit[" << i << "] minhops(" << i << ",j)=";
        for (size_t j=0; j<nq; j++)
        {
            std::cout << MinHops(i,j) << " ";
        }
        std::cout << std::endl;
    }
}

// init multi-core attributes
void InitCores()
{
    if (platformp->topology.count("number_of_cores") <= 0)
    {
        ncores = 1;
        DOUT("Number of cores (topology.number_of_cores) not defined");
    }
    else
    {
        ncores = platformp->topology["number_of_cores"];
    }
    DOUT("Numer of cores= " << ncores);
}

// init x, and y maps
void InitXY()
{
    if (form != gf_irregular)
    {
	    if (platformp->topology.count("qubits") == 0)
	    {
	        FATAL("Regular configuration doesn't specify qubits and their coordinates");
	    }
	    else
	    {
	        for (auto & aqbit : platformp->topology["qubits"] )
	        {
	            size_t qi = aqbit["id"];
	            int qx = aqbit["x"];
	            int qy = aqbit["y"];
	    
	            // sanity checks
	            if ( !(0<=qi && qi<nq) )
	            {
	                FATAL(" qbit in platform topology with id=" << qi << " is configured with id that is not in the range 0..nq-1 with nq=" << nq);
	            }
	            if (x.count(qi) > 0)
	            {
	                FATAL(" qbit in platform topology with id=" << qi << ": duplicate definition of x coordinate");
	            }
	            if (y.count(qi) > 0)
	            {
	                FATAL(" qbit in platform topology with id=" << qi << ": duplicate definition of y coordinate");
	            }
	            if ( !(0<=qx && qx<nx) )
	            {
	                FATAL(" qbit in platform topology with id=" << qi << " is configured with x that is not in the range 0..x_size-1 with x_size=" << nx);
	            }
	            if ( !(0<=qy && qy<ny) )
	            {
	                FATAL(" qbit in platform topology with id=" << qi << " is configured with y that is not in the range 0..y_size-1 with y_size=" << ny);
	            }
	
	            x[qi] = qx;
	            y[qi] = qy;
	        }
	    }
    }
}

// init nbs map
void InitNbs()
{
    if (platformp->topology.count("connectivity") <= 0)
    {
        DOUT("Configuration doesn't specify topology.connectivity: assuming connectivity is specified by edges section");
        conn = gc_specified;
    }
    else
    {
	    std::string connstr;
	    connstr = platformp->topology["connectivity"].get<std::string>();
	    if (connstr == "specified") { conn = gc_specified; }
        else if (connstr == "full") { conn = gc_full; }
        else FATAL("connectivity " << connstr << " not supported");
        DOUT("topology.connectivity=" << connstr );
    }
    if (conn == gc_specified)
    {
	    if (platformp->topology.count("edges") == 0)
	    {
            FATAL(" There aren't edges configured in the platform's topology");
	    }
	    for (auto & anedge : platformp->topology["edges"] )
	    {
            DOUT("connectivity is specified by edges section, reading ...");
	        size_t qs = anedge["src"];
	        size_t qd = anedge["dst"];
	
	        // sanity checks
	        if ( !(0<=qs && qs<nq) )
	        {
	            FATAL(" edge in platform topology has src=" << qs << " that is not in the range 0..nq-1 with nq=" << nq);
	        }
	        if ( !(0<=qd && qd<nq) )
	        {
	            FATAL(" edge in platform topology has dst=" << qd << " that is not in the range 0..nq-1 with nq=" << nq);
	        }
	        for (auto & n : nbs[qs])
	        {
	            if (n == qd)
	            {
	                FATAL(" redefinition of edge with src=" << qs << " and dst=" << qd);
	            }
	        }
	
	        nbs[qs].push_back(qd);
            DOUT("connectivity has been stored in nbs map");
	    }
    }
    if (conn == gc_full)
    {
        DOUT("connectivity is full");
        for (size_t qs=0; qs<nq; qs++)
        {
            for (size_t qd=0; qd<nq; qd++)
	        {
                if (qs != qd)
                {
                    DOUT("connecting qubit[" << qs << "] to qubit[" << qd << "]");
                    nbs[qs].push_back(qd);
                }
            }
	    }
    }
}

void SortNbs()
{
    if (form != gf_irregular)
    {
	    // sort neighbor list by angles
	    for (size_t qi=0; qi<nq; qi++)
	    {
	        // sort nbs[qi] to have increasing clockwise angles around qi, starting with angle 0 at 12:00
	        nbs[qi].sort(
	            [this,qi](const size_t& i, const size_t& j)
	            {
	                return Angle(x[qi], y[qi], x[i], y[i]) < Angle(x[qi], y[qi], x[j], y[j]);
	            }
	        );
	    }
    }
}

};  // end class Grid


// =========================================================================================
// Virt2Real: map of a virtual qubit index to its real qubit index
//
// Mapping maps each used virtual qubit to a real qubit index, but which one that is, may change.
// For a 2-qubit gate its operands should be nearest neighbor; when its virtual operand qubits
// are not mapping to nearest neighbors, that should be accomplished by moving/swapping
// the virtual qubits from their current real qubits to real qubits that are nearest neighbors:
// those moves/swaps are inserted just before that 2-qubit gate.
// Anyhow, the virtual operand qubits of gates must be mapped to the real ones, holding their state.
//
// The number of virtual qubits is less equal than the number of real qubits,
// so their indices use the same data type (size_t) and the same range type 0<=index<nq.
//
// Virt2Real maintains two maps:
// - a map (v2rMap[]) for each virtual qubit that is in use to its current real qubit index.
//      Virtual qubits are in use as soon as they have been encountered as operands in the program.
//      When a virtual qubit is not in use, it maps to UNDEFINED_QUBIT, the undefined real index.
//      The reverse map (GetVirt()) is implemented by a reverse look-up:
//      when there is no virtual qubit that maps to a particular real qubit,
//      the reverse map maps the real qubit index to UNDEFINED_QUBIT, the undefined virtual index.
//      At any time, the virtual to real and reverse maps are 1-1 for qubits that are in use.
// - a map for each real qubit whether there is state in it, and, if so, which (rs[]).
//      When a gate (except for swap/move) has been executed on a real qubit,
//      its state becomes valuable and must be preserved (rs_hasstate below).
//      But before that, it can be in a garbage state (rs_nostate below) or in a known state (rs_wasinited below).
//      The latter is used to replace a swap using a real qubit with such state by a move, which is cheaper.
// There is no support yet to make a virtual qubit not in use (which could be after a measure),
// nor to bring a real qubit in the rs_wasinited or rs_nostate state (perhaps after measure or prep).
//
// Some special situations are worth mentioning:
// - while a virtual qubit is being swapped/moved near to an other one,
//      along the trip real qubits may be used which have no virtual qubit mapping to them;
//      a move can then be used which assumes the 2nd real operand in the |0> (inited) state, and leaves
//      the 1st real operand in that state (while the 2nd has assumed the state of the former 1st).
//      the mapper implementation assumes that all real qubits in the rs_wasinited state are in that state.
// - on program start, no virtual qubit has a mapping yet to a real qubit;
//      mapping is initialized while virtual qubits are encountered as operands.
// - with multiple kernels, kernels assume the (unified) mapping from their predecessors and leave
//      the result mapping to their successors in the kernels' Control Flow Graph;
//      i.e. Virt2Real is what is passed between kernels as dynamic state;
//      statically, the grid, the maximum number of real qubits and the current platform stay unchanged.
// - while evaluating sets of swaps/moves as variations to continue mapping, Virt2Real is passed along
//      to represent the mapping state after such swaps/moves where done; when deciding on a particular
//      variation, the v2r mapping in the mainPast is made to replect the swaps/moves done.

typedef
enum realstate {
    rs_nostate,     // real qubit has no relevant state needing preservation, i.e. is garbage
    rs_wasinited,   // real qubit has initialized state suitable for replacing swap by move
    rs_hasstate     // real qubit has a unique state which must be preserved
} realstate_t;

#define UNDEFINED_QUBIT    MAX_CYCLE

class Virt2Real
{
private:

    size_t              nq;                 // size of the map; after initialization, will always be the same
    std::vector<size_t> v2rMap;             // v2rMap[virtual qubit index] -> real qubit index | UNDEFINED_QUBIT
    std::vector<realstate_t>rs;             // rs[real qubit index] -> {nostate|wasinited|hasstate}

public:

// map real qubit to the virtual qubit index that is mapped to it (i.e. backward map);
// when none, return UNDEFINED_QUBIT;
// a second vector next to v2rMap (i.e. an r2vMap) would speed this up;
size_t GetVirt(size_t r)
{
    MapperAssert(r != UNDEFINED_QUBIT);
    for (size_t v=0; v<nq; v++)
    {
        if (v2rMap[v] == r) return v;
    }
    return UNDEFINED_QUBIT;
}

realstate_t GetRs(size_t q)
{
    return rs[q];
}

void SetRs(size_t q, realstate_t rsvalue)
{
    rs[q] = rsvalue;
}

// expand to desired size
//
// mapping starts off undefined for all virtual qubits
// (unless option mapinitone2one is set, then virtual qubit i maps to real qubit i for all qubits)
//
// real qubits are assumed to have a garbage state
// (unless option mapassumezeroinitstate was set,
//  then all real qubits are assumed to have a state suitable for replacing swap by move)
//
// the rs initializations are done only once, for a whole program
void Init(size_t n)
{
    auto mapinitone2oneopt = ql::options::get("mapinitone2one");
    auto mapassumezeroinitstate = ql::options::get("mapassumezeroinitstate");

    nq = n;
    if ("yes" == mapinitone2oneopt)
    {
        DOUT("Virt2Real::Init(n=" << nq << "), initializing 1-1 mapping");
    }
    else
    {
        DOUT("Virt2Real::Init(n=" << nq << "), initializing on demand mapping");
    }
    if ("yes" == mapassumezeroinitstate)
    {
        DOUT("Virt2Real::Init(n=" << nq << "), assume all qubits in initialized state");
    }
    else
    {
        DOUT("Virt2Real::Init(n=" << nq << "), assume all qubits in garbage state");
    }
    v2rMap.resize(nq);
    rs.resize(nq);
    for (size_t i=0; i<nq; i++)
    {
        if ("yes" == mapinitone2oneopt)
        {
            v2rMap[i] = i;
        }
        else
        {
            v2rMap[i] = UNDEFINED_QUBIT;
        }
        if ("yes" == mapassumezeroinitstate)
        {
            rs[i] = rs_wasinited;
        }
        else
        {
            rs[i] = rs_nostate;
        }
    }
}

// map virtual qubit index to real qubit index
size_t& operator[] (size_t v)
{
    MapperAssert(v < nq);   // implies v != UNDEFINED_QUBIT
    return v2rMap[v];
}

// allocate a new real qubit for an unmapped virtual qubit v (i.e. v2rMap[v] == UNDEFINED_QUBIT);
// note that this may consult the grid or future gates to find a best real
// and thus should not be in Virt2Real but higher up
size_t AllocQubit(size_t v)
{
    // check all real indices for being in v2rMap
    // first one that isn't, is free and is returned
    for (size_t r=0; r<nq; r++)
    {
        size_t vt;
        for (vt=0; vt<nq; vt++)
        {
            if (v2rMap[vt] == r)
            {
                break;
            }
        }
        if (vt >= nq)
        {
            // real qubit r was not found in v2rMap
            // use it to map v
            v2rMap[v] = r;
            MapperAssert(rs[r]==rs_wasinited||rs[r]==rs_nostate);
            DOUT("AllocQubit(v=" << v << ") in r=" << r);
            return r;
        }
    }
    MapperAssert(0);    // number of virt qubits <= number of real qubits
    return UNDEFINED_QUBIT;
}

// r0 and r1 are real qubit indices;
// by execution of a swap(r0,r1), their states are exchanged at runtime;
// so when v0 was in r0 and v1 was in r1, then v0 is now in r1 and v1 is in r0;
// update v2r accordingly
void Swap(size_t r0, size_t r1)
{
    MapperAssert(r0 != r1);
    size_t v0 = GetVirt(r0);
    size_t v1 = GetVirt(r1);
    // DOUT("... swap between ("<< v0<<"<->"<<r0<<","<<v1<<"<->"<<r1<<") and ("<<v0<<"<->"<<r1<<","<<v1<<"<->"<<r0<<" )");
    // DPRINT("... before swap");
    MapperAssert(v0 != v1);         // also holds when vi == UNDEFINED_QUBIT

    if (v0 == UNDEFINED_QUBIT)
    {
        MapperAssert(rs[r0] != rs_hasstate);
    }
    else
    {
        MapperAssert(v0 < nq);
        v2rMap[v0] = r1;
    }

    if (v1 == UNDEFINED_QUBIT)
    {
        MapperAssert(rs[r1] != rs_hasstate);
    }
    else
    {
        MapperAssert(v1 < nq);
        v2rMap[v1] = r0;
    }

    realstate_t ts = rs[r0];
    rs[r0] = rs[r1];
    rs[r1] = ts;
    // DPRINT("... after swap");
}

void DPRINTReal(size_t r)
{
    if ( ql::utils::logger::LOG_LEVEL >= ql::utils::logger::log_level_t::LOG_DEBUG )
        PrintReal(r);
}
void PrintReal(size_t r)
{
    std::cout << " (r" << r;
    switch(rs[r])
    {
    case rs_nostate:
        std::cout << ":no";
        break;
    case rs_wasinited:
        std::cout << ":in";
        break;
    case rs_hasstate:
        std::cout << ":st";
        break;
    default:
        MapperAssert(0);
    }
    size_t v = GetVirt(r);
    if (v == UNDEFINED_QUBIT)
    {
        std::cout << "<-UN)";
    }
    else
    {
        std::cout << "<-v" << v << ")";
    }
}

void PrintVirt(size_t v)
{
    std::cout << " (v" << v;
    size_t r = v2rMap[v];
    if (r == UNDEFINED_QUBIT)
    {
        std::cout << "->UN)";
    }
    else
    {
        std::cout << "->r" << r;
        switch(rs[r])
        {
        case rs_nostate:
            std::cout << ":no)";
            break;
        case rs_wasinited:
            std::cout << ":in)";
            break;
        case rs_hasstate:
            std::cout << ":st)";
            break;
        default:
            MapperAssert(0);
        }
    }
}

void DPRINTReal(std::string s, size_t r0, size_t r1)
{
    if ( ql::utils::logger::LOG_LEVEL >= ql::utils::logger::log_level_t::LOG_DEBUG )
        PrintReal(s, r0, r1);
}
void PrintReal(std::string s, size_t r0, size_t r1)
{
    // DOUT("v2r.PrintReal ...");
    std::cout << s << ":";
//  std::cout << "... real2Virt(r<-v) " << s << ":";

    PrintReal(r0);
    PrintReal(r1);
    std::cout << std::endl;
}

void DPRINT(std::string s)
{
    if ( ql::utils::logger::LOG_LEVEL >= ql::utils::logger::log_level_t::LOG_DEBUG )
        Print(s);
}
void Print(std::string s)
{
    // DOUT("v2r.Print ...");
    std::cout << s << ":";
//  std::cout << "... virt2Real(r<-v) " << s << ":";
    for (size_t v=0; v<nq; v++)
    {
        PrintVirt(v);
    }
    std::cout << std::endl;

#ifdef needed
    std::cout << "... real2virt(r->v) " << s << ":";
    for (size_t r=0; r<nq; r++)
    {
        PrintReal(r);
    }
    std::cout << std::endl;
#endif
}

void Export(std::vector<size_t>& kv2rMap)
{
    kv2rMap = v2rMap;
}

void Export(std::vector<int>& krs)
{
    krs.resize(rs.size());
    for (unsigned int i=0; i < rs.size(); i++)
    {
        krs[i] = (int)rs[i];
    }
}

};  // end class Virt2Real





// =========================================================================================
// FreeCycle: map each real qubit to the first cycle that it is free for use
//
// in scheduling gates, qubit dependencies cause latencies
// for each real qubit, the first cycle that it is free to use is the cycle that the
// last gate that was scheduled in the qubit, has just finished (i.e. in the previous cycle);
// the map serves as a summary to ease scheduling next gates
//
// likewise, while mapping, swaps are scheduled just before a non-NN two-qubit gate,
// moreover, such swaps may involve real qubits on the path between the real operand qubits of the gate,
// which may be different from the real operand qubits;
// the evaluation of which path of swaps is best is, among other data, based
// on which path causes the latency of the whole circuit to be extended the least;
// this latency extension is measured from the data in the FreeCycle map;
// so a FreeCycle map is part of each path of swaps that is evaluated for a particular non-NN 2-qubit gate
// next to a FreeCycle map that is part of the output stream (the main past)
//
// since gate durations are in nano-seconds, and one cycle is some fixed number of nano-seconds,
// the duration is converted to a rounded-up number of cycles when computing the added latency
class FreeCycle
{
private:

    const ql::quantum_platform   *platformp;// platform description
    size_t                  nq;      // size of the map; after initialization, will always be the same
    size_t                  ct;      // multiplication factor from cycles to nano-seconds (unit of duration)
    std::vector<size_t>     fcv;     // fcv[real qubit index i]: qubit i is free from this cycle on
    ql::arch::resource_manager_t rm;  // actual resources occupied by scheduled gates


// access free cycle value of qubit i
size_t& operator[] (size_t i)
{
    return fcv[i];
}

public:

// explicit FreeCycle constructor
// needed for virgin construction
// default constructor was deleted because it cannot construct resource_manager_t without parameters
FreeCycle()
{
    DOUT("Constructing FreeCycle");
}

void Init(const ql::quantum_platform *p)
{
    DOUT("FreeCycle::Init()");
    ql::arch::resource_manager_t lrm(*p, ql::forward_scheduling);   // allocated here and copied below to rm because of platform parameter
    DOUT("... created FreeCycle Init local resource_manager");
    platformp = p;
    nq = platformp->qubit_number;
    ct = platformp->cycle_time;
    DOUT("... FreeCycle: nq=" << nq << ", ct=" << ct << "), initializing to all 0 cycles");
    fcv.clear();
    fcv.resize(nq, 1);   // this 1 implies that cycle of first gate will be 1 and not 0; OpenQL convention!?!?
    DOUT("... about to copy FreeCycle Init local resource_manager to FreeCycle member rm");
    rm = lrm;
    DOUT("... done copy FreeCycle Init local resource_manager to FreeCycle member rm");
}

// depth of the FreeCycle map
// equals the max of all entries minus the min of all entries
// not used yet; would be used to compute the max size of a top window on the past
size_t Depth()
{
    return Max() - Min();
}

// min of the FreeCycle map equals the min of all entries;
size_t Min()
{
    size_t  minFreeCycle = MAX_CYCLE;
    for (auto& v : fcv)
    {
        if (v < minFreeCycle)
        {
            minFreeCycle = v;
        }
    }
    return minFreeCycle;
}

// max of the FreeCycle map equals the max of all entries;
size_t Max()
{
    size_t  maxFreeCycle = 0;
    for (auto& v : fcv)
    {
        if (maxFreeCycle < v)
        {
            maxFreeCycle = v;
        }
    }
    return maxFreeCycle;
}

void DPRINT(std::string s)
{
    if ( ql::utils::logger::LOG_LEVEL >= ql::utils::logger::log_level_t::LOG_DEBUG )
        Print(s);
}
void Print(std::string s)
{
    size_t  minFreeCycle = Min();
    size_t  maxFreeCycle = Max();
    std::cout << "... FreeCycle" << s << ":";
    for (size_t i=0; i<nq; i++)
    {
        size_t v = fcv[i];
        std::cout << " [" << i << "]=";
        if (v == minFreeCycle)
        {
            std::cout << "_";
        }
        if (v == maxFreeCycle)
        {
            std::cout << "^";
        }
        std::cout << v;
    }
    std::cout << std::endl;
    // rm.Print("... in FreeCycle: ");
}

// return whether gate with first operand qubit r0 can be scheduled earlier than with operand qubit r1
bool IsFirstOperandEarlier(size_t r0, size_t r1)
{
    DOUT("... fcv[" << r0 << "]=" << fcv[r0]<< " fcv[" << r1 << "]=" << fcv[r1] << " IsFirstOperandEarlier=" << (fcv[r0] < fcv[r1]));
    return fcv[r0] < fcv[r1];
}

// will a swap(fr0,fr1) start earlier than a swap(sr0,sr1)?
// is really a short-cut ignoring config file and perhaps several other details
bool IsFirstSwapEarliest(size_t fr0, size_t fr1, size_t sr0, size_t sr1)
{
    std::string mapreverseswapopt = ql::options::get("mapreverseswap");
    if ("yes" == mapreverseswapopt)
    {
        if (fcv[fr0] < fcv[fr1])
        {
            size_t  tmp = fr1; fr1 = fr0; fr0 = tmp;
        }
        if (fcv[sr0] < fcv[sr1])
        {
            size_t  tmp = sr1; sr1 = sr0; sr0 = tmp;
        }
    }
    size_t  startCycleFirstSwap = std::max<size_t>(fcv[fr0]-1, fcv[fr1]);
    size_t  startCycleSecondSwap = std::max<size_t>(fcv[sr0]-1, fcv[sr1]);

    DOUT("... fcv[" << fr0 << "]=" << fcv[fr0] << " fcv[" << fr1 << "]=" << fcv[fr1] << " start=" << startCycleFirstSwap << " fcv[" << sr0 << "]=" << fcv[sr0] << " fcv[" << sr1 << "]=" << fcv[sr1] << " start=" << startCycleSecondSwap << " IsFirstSwapEarliest=" << (startCycleFirstSwap < startCycleSecondSwap));
    return startCycleFirstSwap < startCycleSecondSwap;
}

// when we would schedule gate g, what would be its start cycle? return it
// gate operands are real qubit indices
// is purely functional, doesn't affect state
size_t StartCycleNoRc(ql::gate *g)
{
    auto&       q = g->operands;
    size_t      operandCount = q.size();

    size_t      startCycle;
    if (operandCount == 1)
    {
        startCycle = fcv[q[0]];
    }
    else // if (operandCount == 2)
    {
        startCycle = std::max<size_t>(fcv[q[0]], fcv[q[1]]);
    }
    MapperAssert (startCycle < MAX_CYCLE);

    return startCycle;
}

// when we would schedule gate g, what would be its start cycle? return it
// gate operands are real qubit indices
// is purely functional, doesn't affect state
size_t StartCycle(ql::gate *g)
{
    size_t      startCycle = StartCycleNoRc(g);
    
    auto        mapopt = ql::options::get("mapper");
    if (mapopt == "baserc" || mapopt == "minextendrc")
    {
        size_t      baseStartCycle = startCycle;

        while (startCycle < MAX_CYCLE)
        {
            // DOUT("Startcycle for " << g->qasm() << ": available? at startCycle=" << startCycle);
            if (rm.available(startCycle, g, *platformp))
            {
                // DOUT(" ... [" << startCycle << "] resources available for " << g->qasm());
                break;
            }   
            else
            {
                // DOUT(" ... [" << startCycle << "] Busy resource for " << g->qasm());
                startCycle++;
            }
        }
        if (baseStartCycle != startCycle)
        {
            // DOUT(" ... from [" << baseStartCycle << "] to [" << startCycle-1 << "] busy resource(s) for " << g->qasm());
        }
    }
    MapperAssert (startCycle < MAX_CYCLE);

    return startCycle;
}

// schedule gate g in the FreeCycle map
// gate operands are real qubit indices
// the FreeCycle map is updated, not the resource map
// this is done, because AddNoRc is used to represent just gate dependences, avoiding a build of a dep graph
void AddNoRc(ql::gate *g, size_t startCycle)
{
    auto&       q = g->operands;
    size_t      operandCount = q.size();
    size_t      duration = (g->duration+ct-1)/ct;   // rounded-up unsigned integer division

    if (operandCount == 1)
    {
        fcv[q[0]] = startCycle + duration;
    }
    else // if (operandCount == 2)
    {
        fcv[q[0]] = startCycle + duration;
        fcv[q[1]] = fcv[q[0]];
    }
}

// schedule gate g in the FreeCycle and resource maps
// gate operands are real qubit indices
// both the FreeCycle map and the resource map are updated
// startcycle must be the result of an earlier StartCycle call (with rc!)
void Add(ql::gate *g, size_t startCycle)
{
    AddNoRc(g, startCycle);

    auto        mapopt = ql::options::get("mapper");
    if (mapopt == "baserc" || mapopt == "minextendrc")
    {
        rm.reserve(startCycle, g, *platformp);
    }
}

};  // end class FreeCycle


// =========================================================================================
// Past: state of the mapper while somewhere in the mapping process
//
// there is a Past attached to the output stream, that is a kind of window with a list of gates in it,
// to which gates are added after mapping; this is called the 'main' Past.
// while mapping, several alternatives are evaluated, each of which also has a Past attached,
// and each of which for most of the parts starts off as a copy of the 'main' Past;
// but it is in fact a temporary extension of this main Past
// 
// Past contains gates of which the schedule might influence a future path selected for mapping binary gates
// It maintains for each qubit from which cycle on it is free, so that swap insertion
// can exploit this to hide its overall circuit latency overhead by increasing ILP.
// Also it maintains the 1 to 1 (reversible) virtual to real qubit map: all gates in past
// and beyond are mapped and have real qubits as operands.
// While experimenting with path alternatives, a clone is made of the main past,
// to insert swaps and evaluate the latency effects; note that inserting swaps changes the mapping.
//
// On arrival of a quantum gate(s):
// - [isempty(waitinglg)]
// - if 2q nonNN clone mult. pasts, in each clone Add swap/move gates, Schedule, evaluate clones, select, Add swaps to mainPast
// - Add, Add, ...: add quantum gates to waitinglg, waiting to be scheduled in [!isempty(waitinglg)]
// - Schedule: schedules all quantum gates of waitinglg into lg [isempty(waitinglg) && !isempty(lg)]
// On arrival of a classical gate:
// - FlushAll: lg flushed to outlg [isempty(waitinglg) && isempty(lg) && !isempty(outlg)]
// - ByPass: classical gate added to outlg [isempty(waitinglg) && isempty(lg) && !isempty(outlg)]
// On no gates:
// - [isempty(waitinglg)]
// - FlushAll: lg flushed to outlg [isempty(waitinglg) && isempty(lg) && !isempty(outlg)]
// On end:
// - Out: outlg flushed to outCirc [isempty(waitinglg) && isempty(lg) && isempty(outlg)]
class Past
{
private:

    size_t                  nq;         // width of Past, Virt2Real, UseCount maps in number of real qubits
    size_t                  ct;         // cycle time, multiplier from cycles to nano-seconds
    const ql::quantum_platform    *platformp; // platform describing resources for scheduling
    ql::quantum_kernel      *kernelp;   // current kernel for creating gates

    Virt2Real               v2r;        // state: current Virt2Real map, imported/exported to kernel
    FreeCycle               fc;         // state: FreeCycle map (including resource_manager) of this Past
    typedef ql::gate *      gate_p;
    std::list<gate_p>       waitinglg;  // . . .  list of q gates in this Past, topological order, waiting to be scheduled in
                                        //        waitinglg only contains gates from Add and final Schedule call
                                        //        when evaluating alternatives, it is empty when Past is cloned; so no state
public:
    std::list<gate_p>       lg;         // state: list of q gates in this Past, scheduled by their (start) cycle values
                                        //        so this is the result list of this Past, to compare with other Alters
private:
    std::list<gate_p>       outlg;      // . . .  list of gates flushed out of this Past, not yet put in outCirc
                                        //        when evaluating alternatives, outlg stays constant; so no state
    std::map<gate_p,size_t> cycle;      // state: gate to cycle map, startCycle value of each past gatecycle[gp]
                                        //        cycle[gp] can be different for each gp for each past
                                        //        gp->cycle is not used by MapGates
                                        //        although updated by set_cycle called from MakeAvailable/TakeAvailable
    size_t                  nswapsadded;// number of swaps (including moves) added to this past
    size_t                  nmovesadded;// number of moves added to this past

public:

// explicit Past constructor
// needed for virgin construction
Past()
{
    DOUT("Constructing Past");
}

// past initializer
void Init(const ql::quantum_platform *p, ql::quantum_kernel *k)
{
    DOUT("Past::Init");
    platformp = p;
    kernelp = k;

    nq = platformp->qubit_number;
    ct = platformp->cycle_time;

    MapperAssert(kernelp->c.empty());   // kernelp->c will be used by new_gate to return newly created gates into
    v2r.Init(nq);               // v2r initializtion until v2r is imported from context
    fc.Init(platformp);         // fc starts off with all qubits free, is updated after schedule of each gate
    waitinglg.clear();          // no gates pending to be scheduled in; Add of gate to past entered here
    lg.clear();                 // no gates scheduled yet in this past; after schedule of gate, it gets here
    outlg.clear();              // no gates output yet by flushing from or bypassing this past
    nswapsadded = 0;            // no swaps or moves added yet to this past; AddSwap adds one here
    nmovesadded = 0;            // no moves added yet to this past; AddSwap may add one here
    cycle.clear();              // no gates have cycles assigned in this past; scheduling gate updates this
}

// import Past's v2r from v2r_value
void ImportV2r(Virt2Real& v2r_value)
{
    v2r = v2r_value;
}

// export Past's v2r into v2r_destination
void ExportV2r(Virt2Real& v2r_destination)
{
    v2r_destination = v2r;
}

void DFcPrint()
{
    if ( ql::utils::logger::LOG_LEVEL >= ql::utils::logger::log_level_t::LOG_DEBUG )
        fc.Print("");
}
void FcPrint()
{
    fc.Print("");
}

void Print(std::string s)
{
    std::cout << "... Past " << s << ":";
    v2r.Print("");
    fc.Print("");
    // DOUT("... list of gates in past");
    for ( auto & gp: lg)
    {
        DOUT("[" << cycle[gp] << "] " << gp->qasm());
    }
}

// all gates in past.waitinglg are scheduled here into past.lg
// note that these gates all are mapped and so have real operand qubit indices
// the FreeCycle map reflects for each qubit the first free cycle
// all new gates, now in waitinglist, get such a cycle assigned below, increased gradually, until definitive
void Schedule()
        // the copy includes the resource manager.
{
    // DOUT("Schedule ...");

    while (!waitinglg.empty())
    {
        size_t      startCycle = MAX_CYCLE;
        gate_p      gp;

        // find the gate with the minimum startCycle
        //
        // IMPORTANT: this assumes that the waitinglg gates list is in topological order,
        // which is ok because the pair of swap lists use distict qubits and
        // the gates of each are added to the back of the list in the order of execution.
        // Using tryfc.Add, the tryfc (try FreeCycle map) reflects the earliest startCycle per qubit,
        // and so dependences are respected, so we can find the gate that can start first ...
        // Note that tryfc includes the free cycle vector AND the resource map,
        // so using tryfc.StartCycle/tryfc.Add we get a realistic ASAP rc schedule.
        // We use a copy of fc and not fc itself, since the latter reflects the really scheduled gates
        // and that shouldn't be changed.
        //
        // This search is really a hack to avoid
        // the construction of a dependence graph and a set of schedulable gates
        FreeCycle   tryfc = fc;
        for (auto & trygp : waitinglg)
        {
            size_t tryStartCycle = tryfc.StartCycle(trygp);
            tryfc.Add(trygp, tryStartCycle);

            if (tryStartCycle < startCycle)
            {
                startCycle = tryStartCycle;
                gp = trygp;
            }
        }

        // add this gate to the maps, scheduling the gate (doing the cycle assignment)
        // DOUT("... add " << gp->qasm() << " startcycle=" << startCycle << " cycles=" << ((gp->duration+ct-1)/ct) );
        fc.Add(gp, startCycle);
        cycle[gp] = startCycle; // cycle[gp] is private to this past but gp->cycle is private to gp
        gp->cycle = startCycle; // so gp->cycle gets assigned for each alter' Past and finally definitively for mainPast
        // DOUT("... set " << gp->qasm() << " at cycle " << startCycle);
    
        // insert gate gp in lg, the list of gates, in cycle[gp] order, and inside this order, as late as possible
        //
        // reverse iterate because the insertion is near the end of the list
        // insert so that cycle values are in order afterwards and the new one is nearest to the end
        std::list<gate_p>::reverse_iterator rigp = lg.rbegin();
        for (; rigp != lg.rend(); rigp++)
        {
            if (cycle[*rigp] <= startCycle)
            {
                // rigp.base() because insert doesn't work with reverse iteration
                // rigp.base points after the element that rigp is pointing at
                // which is lucky because insert only inserts before the given element
                // the end effect is inserting after rigp
                lg.insert(rigp.base(), gp);
                break;
            }
        }
        // when list was empty or no element was found, just put it in front
        if (rigp == lg.rend())
        {
            lg.push_front(gp);
        }
    
        // having added it to the main list, remove it from the waiting list
        waitinglg.remove(gp);
    }

    // DPRINT("Schedule:");
}

// compute costs in cycle extension of optionally scheduling initcirc before the inevitable circ
int InsertionCost(ql::circuit& initcirc, ql::circuit& circ)
{
     // first fake-schedule initcirc followed by circ in a private freecyclemap
     size_t initmax;
     FreeCycle   tryfcinit = fc;
     for (auto & trygp : initcirc)
     {
         size_t tryStartCycle = tryfcinit.StartCycleNoRc(trygp);
         tryfcinit.AddNoRc(trygp, tryStartCycle);
     }
     for (auto & trygp : circ)
     {
         size_t tryStartCycle = tryfcinit.StartCycleNoRc(trygp);
         tryfcinit.AddNoRc(trygp, tryStartCycle);
     }
     initmax = tryfcinit.Max(); // this reflects the depth afterwards

     // then fake-schedule circ alone in a private freecyclemap
     size_t max;
     FreeCycle   tryfc = fc;
     for (auto & trygp : circ)
     {
         size_t tryStartCycle = tryfc.StartCycleNoRc(trygp);
         tryfc.AddNoRc(trygp, tryStartCycle);
     }
     max = tryfc.Max();         // this reflects the depth afterwards

     DOUT("... scheduling init+circ => depth " << initmax << ", scheduling circ => depth " << max << ", init insertion cost " << (initmax-max));
     MapperAssert(initmax >= max);
     // scheduling initcirc would be for free when initmax == max, so the cost is (initmax - max)
     return (initmax - max);
}

// add the mapped gate to the current past
// means adding it to the current past's waiting list, waiting for it to be scheduled later
void Add(gate_p gp)
{
    waitinglg.push_back(gp);
}

// create a new gate with given name and qubits
// return whether this was successful
// return the created gate(s) in circ (which is supposed to be empty on entry)
//
// since kernel.h only provides a gate interface as method of class quantum_kernel, that adds the gate to kernel.c,
// and we want the gate (or its decomposed sequence) here to be added to circ,
// the kludge is implemented to make sure that kernel.c (the current kernel's mapper input/output circuit)
// is available for this:
// in class Future, kernel.c is copied into the dependence graph or copied to a local circuit; and
// in Mapper::MapCircuit, a temporary local output circuit is used, which is written to kernel.c only at the very end
bool new_gate(ql::circuit& circ, std::string gname, std::vector<size_t> qubits, std::vector<size_t> cregs = {}, size_t duration=0, double angle = 0.0)

{
    bool    added;
    MapperAssert(circ.empty());
    MapperAssert(kernelp->c.empty());
    added = kernelp->gate_nonfatal(gname, qubits, cregs, duration, angle);   // creates gates in kernelp->c
    circ = kernelp->c;
    kernelp->c.clear();
    for (auto gp: circ)
    {
        DOUT("new_gate added: " << gp->qasm());
    }
    return added;
}

// return number of swaps added to this past
size_t NumberOfSwapsAdded()
{
    return nswapsadded;
}

// return number of moves added to this past
size_t NumberOfMovesAdded()
{
    return nmovesadded;
}

void new_gate_exception(std::string s)
{
    FATAL("gate is not supported by the target platform: '" << s << "'");
}

// will a swap(fr0,fr1) start earlier than a swap(sr0,sr1)?
// is really a short-cut ignoring config file and perhaps several other details
bool IsFirstSwapEarliest(size_t fr0, size_t fr1, size_t sr0, size_t sr1)
{
    return fc.IsFirstSwapEarliest(fr0, fr1, sr0, sr1);
}

// generate a move into circ with parameters r0 and r1 (which GenMove may reverse)
// whether this was successfully done can be seen from whether circ was extended
// please note that the reversal of operands may have been done also when GenMove was not successful
void GenMove(ql::circuit& circ, size_t& r0, size_t& r1)
{
    if (v2r.GetRs(r0)!=rs_hasstate)
    {
        MapperAssert(v2r.GetRs(r0)==rs_nostate||v2r.GetRs(r0)==rs_wasinited);
        // interchange r0 and r1, so that r1 (right-hand operand of move) will be the state-less one
        size_t  tmp = r1; r1 = r0; r0 = tmp;
        // DOUT("... reversed operands for move to become move(q" << r0 << ",q" << r1 << ") ...");
    }
    MapperAssert (v2r.GetRs(r0)==rs_hasstate);    // and r0 will be the one with state
    MapperAssert (v2r.GetRs(r1)!=rs_hasstate);    // and r1 will be the one without state (rs_nostate || rs_wasinited)

    // first (optimistically) create the move circuit and add it to circ
    bool created;
    auto mapperopt = ql::options::get("mapper");
    if ("maxfidelity" == mapperopt)
    {
        created = new_gate(circ, "move_prim", {r0,r1});    // gates implementing move returned in circ
    }
    else
    {
        created = new_gate(circ, "move_real", {r0,r1});    // gates implementing move returned in circ
    }
    if (!created)
    {
        created = new_gate(circ, "move", {r0,r1});
        if (!created) new_gate_exception("move or move_real");
    }

    if (v2r.GetRs(r1) == rs_nostate)
    {
        // r1 is not in inited state, generate in initcirc the circuit to do so
        // DOUT("... initializing non-inited " << r1 << " to |0> (inited) state preferably using move_init ...");
        ql::circuit initcirc;

        created = new_gate(initcirc, "move_init", {r1});
        if (!created)
        {
            created = new_gate(initcirc, "prepz", {r1});
            // if (created)
            // {
            //     created = new_gate(initcirc, "h", {r1});
            //     if (!created) new_gate_exception("h");
            // }
            if (!created) new_gate_exception("move_init or prepz");
        }

        // when difference in extending circuit after scheduling initcirc+circ or just circ
        // is less equal than threshold cycles (0 would mean scheduling initcirc was for free),
        // commit to it, otherwise abort
        int threshold;
        std::string mapusemovesopt = ql::options::get("mapusemoves");
        if ("yes" == mapusemovesopt)
        {
            threshold = 0;
        }
        else
        {
            threshold = atoi(mapusemovesopt.c_str());
        }
        if (InsertionCost(initcirc, circ) <= threshold)
        {
            // so we go for it!
            // circ contains move; it must get the initcirc before it ...
            // do this by appending circ's gates to initcirc, and then swapping circ and initcirc content
            DOUT("... initialization is for free, do it ...");
            for (auto& gp : circ)
            {
                initcirc.push_back(gp);
            }
            circ.swap(initcirc);
            v2r.SetRs(r1, rs_wasinited);
        }
        else
        {
            // undo damage done, will not do move but swap, i.e. nothing created thisfar
            DOUT("... initialization extends circuit, don't do it ...");
            circ.clear();       // circ being cleared also indicates creation wasn't successful
        }
        // initcirc getting out-of-scope here so gets destroyed
    }
}

// generate a single swap/move with real operands and add it to the current past's waiting list;
// note that the swap/move may be implemented by a series of gates (circuit circ below),
// and that a swap/move essentially is a commutative operation, interchanging the states of the two qubits;
//
// a move is implemented by 2 CNOTs, while a swap by 3 CNOTs, provided the target qubit is in |0> (inited) state;
// so, when one of the operands is the current location of an unused virtual qubit,
// use a move with that location as 2nd operand,
// after first having initialized the target qubit in |0> (inited) state when that has not been done already;
// but this initialization must not extend the depth so can only be done when cycles for it are for free
void AddSwap(size_t r0, size_t r1)
{
    bool created = false;

    DOUT("... extending with swap(q" << r0 << ",q" << r1 << ") ...");
    v2r.DPRINTReal("... adding swap/move", r0, r1);

    MapperAssert(v2r.GetRs(r0)==rs_wasinited||v2r.GetRs(r0)==rs_nostate||v2r.GetRs(r0)==rs_hasstate);
    MapperAssert(v2r.GetRs(r1)==rs_wasinited||v2r.GetRs(r1)==rs_nostate||v2r.GetRs(r1)==rs_hasstate);

    if (v2r.GetRs(r0)!=rs_hasstate && v2r.GetRs(r1)!=rs_hasstate)
    {
        DOUT("... no state in both operand of intended swap/move; don't add swap/move gates");
        v2r.Swap(r0,r1);
        return;
    }

    ql::circuit circ;   // current kernel copy, clear circuit
    std::string mapusemovesopt = ql::options::get("mapusemoves");
    if ("no" != mapusemovesopt && (v2r.GetRs(r0)!=rs_hasstate || v2r.GetRs(r1)!=rs_hasstate))
    {
        GenMove(circ, r0, r1);
        created = circ.size()!=0;
        if (created)
        {
            // generated move
            // move is in circ, optionally with initialization in front of it
            // also rs of its 2nd operand is 'rs_wasinited'
            // note that after swap/move, r0 will be in this state then
            nmovesadded++;                       // for reporting at the end
            DOUT("... move(q" << r0 << ",q" << r1 << ") ...");
        }
        else
        {
            DOUT("... move(q" << r0 << ",q" << r1 << ") cancelled, go for swap");
        }
    }
    if (!created)
    {
        // no move generated so do swap
        std::string mapreverseswapopt = ql::options::get("mapreverseswap");
        if ("yes" == mapreverseswapopt)
        {
            // swap(r0,r1) is about to be generated
            // it is functionally symmetrical,
            // but in the implementation r1 starts 1 cycle earlier than r0 (we should derive this from json file ...)
            // so swap(r0,r1) with interchanged operands might get scheduled 1 cycle earlier;
            // when fcv[r0] < fcv[r1], r0 is free for use 1 cycle earlier than r1, so a reversal will help
            if (fc.IsFirstOperandEarlier(r0, r1))
            {
                size_t  tmp = r1; r1 = r0; r0 = tmp;
                DOUT("... reversed swap to become swap(q" << r0 << ",q" << r1 << ") ...");
            }
        }
        auto mapperopt = ql::options::get("mapper");
        if ("maxfidelity" == mapperopt)
        {
            created = new_gate(circ, "swap_prim", {r0,r1});    // gates implementing swap returned in circ
        }
        else
        {
            created = new_gate(circ, "swap_real", {r0,r1});    // gates implementing swap returned in circ
        }
        if (!created)
        {
            created = new_gate(circ, "swap", {r0,r1});
            if (!created) new_gate_exception("swap or swap_real");
        }
        DOUT("... swap(q" << r0 << ",q" << r1 << ") ...");
    }
    nswapsadded++;                       // for reporting at the end
    for (auto &gp : circ)
    {
        Add(gp);
    }
    v2r.Swap(r0,r1);        // reflect in v2r that r0 and r1 interchanged state, i.e. update the map to reflect the swap
}

// add the mapped gate (with real qubit indices as operands) to the past
// by adding it to the waitinglist and scheduling it into the past
void AddAndSchedule(gate_p gp)
{
    Add(gp);
    Schedule();
}

// find real qubit index implementing virtual qubit index;
// if not yet mapped, allocate a new real qubit index and map to it
size_t MapQubit(size_t v)
{
    size_t  r = v2r[v];
    if (r == UNDEFINED_QUBIT)
    {
        r = v2r.AllocQubit(v);
    }
    return r;
}

// MakeReal gp
// assume gp points to a virtual gate with virtual qubit indices as operands;
// when a gate can be created with the same name but with "_real" appended, with the real qubits as operands, then create that gate
// otherwise keep the old gate; replace the virtual qubit operands by the real qubit indices
// since creating a new gate may result in a decomposition to several gates, the result is returned as a circuit vector
//
// So each gate in the circuit (optionally) passes through the following phases:
// 1. it is created:
//      when a decomposition in config file, decompose immediately, otherwise just create (k.gate)
//      so we expect gates like: x, cz, cnot to be specified in config file;
//      on the resulting (decomposed) gates, the routing is done including depth/cost estimation
// 2a.if needed for mapping, swap/move is created:
//      first try creating swap_real/move_real as above, otherwise just swap/real (AddSwap)
//      so we expect gates like: swap_real, move_real to be specified in config file,
//      swap_real/move_real unlike swap/real allow immediate decomposition;
//      when no swap_real/move_real are specified, just swap/move must be present
//      and swap/move are created usually without decomposition;
//      on the resulting (decomposed) gates, the routing is done including depth/cost estimation;
//      when the resulting gates end in _prim, see step 3
// 2b.the resulting gates of step 1: map operands/gate:
//      first try creating gate_real as above, otherwise just gate (MakeReal)
//      gate_real unlike gate allows immediate decomposition;
//      when the resulting gates end in _prim, see step 3
// 3. make primitive gates:
//      for each gate try recreating it with _prim appended to its name, otherwise keep it; this decomposes those with corresponding _prim entries
// 4. final schedule:
//      the resulting gates are subject to final scheduling (the original resource-constrained scheduler)

void stripname(std::string& name)
{
    DOUT("stripname(name=" << name << ")");
    size_t p = name.find(" ");
    if (p != std::string::npos)
    {
        name = name.substr(0,p);
    }
    DOUT("... after stripname name=" << name);
}

void MakeReal(ql::gate* gp, ql::circuit& circ)
{
    DOUT("MakeReal: " << gp->qasm());

    std::string gname = gp->name;
    stripname(gname);

    std::vector<size_t> real_qubits  = gp->operands;// starts off as copy of virtual qubits!
    for (auto& qi : real_qubits)
    {
        qi = MapQubit(qi);          // and now they are real
        auto mapprepinitsstateopt = ql::options::get("mapprepinitsstate");
        if (mapprepinitsstateopt == "yes" && (gname == "prepz" || gname == "Prepz"))
        {
            v2r.SetRs(qi, rs_wasinited);
        }
        else
        {
            v2r.SetRs(qi, rs_hasstate);
        }
    }

    auto mapperopt = ql::options::get("mapper");
    std::string real_gname = gname;
    if ("maxfidelity" == mapperopt)
    {
        DOUT("MakeReal: with mapper==maxfidelity generate _prim");
        real_gname.append("_prim");
    }
    else
    {
        real_gname.append("_real");
    }

    bool created = new_gate(circ, real_gname, real_qubits, gp->creg_operands, gp->duration, gp->angle);
    if (!created)
    {
        created = new_gate(circ, gname, real_qubits, gp->creg_operands, gp->duration, gp->angle);
        if (!created)
        {
            FATAL("MakeReal: failed creating gate " << real_gname << " or " << gname);
        }
    }
    DOUT("... MakeReal: new gate created for: " << real_gname << " or " << gname);
}

// as mapper after-burner
// make primitives of all gates that also have an entry with _prim appended to its name
// and decomposing it according to the .json file gate decomposition
void MakePrimitive(ql::gate* gp, ql::circuit& circ)
{
    std::string gname = gp->name;
    stripname(gname);
    std::string prim_gname = gname;
    prim_gname.append("_prim");
    bool created = new_gate(circ, prim_gname, gp->operands, gp->creg_operands, gp->duration, gp->angle);
    if (!created)
    {
        created = new_gate(circ, gname, gp->operands, gp->creg_operands, gp->duration, gp->angle);
        if (!created)
        {
            FATAL("MakePrimtive: failed creating gate " << prim_gname << " or " << gname);
        }
    }
    DOUT("... MakePrimtive: new gate created for: " << prim_gname << " or " << gname);
}

size_t MaxFreeCycle()
{
    return fc.Max();
}

// nonq and q gates follow separate flows through Past:
// - q gates are put in waitinglg when added and then scheduled; and then ordered by cycle into lg
//      in lg they are waiting to be inspected and scheduled, until [too many are there,] a nonq comes or end-of-circuit
// - nonq gates first cause lg to be flushed/cleared to output before the nonq gate is output
// all gates in outlg are out of view for scheduling/mapping optimization and can be taken out to elsewhere
void FlushAll()
{
    for( auto & gp : lg )
    {
        outlg.push_back(gp);
    }
    lg.clear();         // so effectively, lg's content was moved to outlg

    // fc.Init(platformp); // needed?
    // cycle.clear();      // needed?
                        // cycle is initialized to empty map
                        // is ok without windowing, but with window, just delete the ones outside the window
}

// gp as nonq gate immediately goes to outlg
void ByPass(ql::gate* gp)
{
    if (!lg.empty())
    {
        FlushAll();
    }
    outlg.push_back(gp);
}

// mainPast flushes outlg to parameter oc
void Out(ql::circuit& oc)
{
    for( auto & gp : outlg )
    {
        oc.push_back(gp);
    }
    outlg.clear();
}

};  // end class Past


// =========================================================================================
// Alter: one alternative way to make two real qbits (operands of a 2-qubit gate) nearest neighbor (NN);
// of these two qubits, the first qubit is called the source, the second is called the target qubit.
// The Alter stores a series of real qubit indices; qubits/indices are equivalent to the nodes in the grid.
// An Alter represents a 2-qubit gate and a path through the grid from source to target qubit, with each hop between
// qubits/nodes only between neighboring nodes in the grid; the intention is that all but one hops
// translate into swaps and that one hop remains that will be the place to do the 2-qubit gate.
//
// Actually, the Alter goes through several stages:
// - first, for the given 2-qubit gate that is stored in targetgp,
//   while finding a path from its source to its target, the current path is kept in total;
//   fromSource, fromTarget, past and score are not used; past is a clone of the main past
// - paths are found starting from the source node, and aiming to reach the target node,
//   each time adding one additional hop to the path
//   fromSource, fromTarget, and score are still empty and not used
// - each time another continuation of the path is found, the current Alter is cloned
//   and the difference continuation represented in the total attribute; it all starts with an empty Alter
//   fromSource, fromTarget, and score are still empty and not used
// - once all alternative total paths for the 2-qubit gate from source to target have been found
//   each of these is split again in all possible ways (to ILP overlap swaps from source and target);
//   the split is the place where the two-qubit gate is put
// - the alternative splits are made separate Alters and for each
//   of these the two partial paths are stored in fromSource and fromTarget;
//   a partial path stores its starting and end nodes (so contains 1 hop less than its length);
//   the partial path of the target operand is reversed, so starts at the target qubit
// - then we add swaps to past following the recipee in fromSource and fromTarget; this extends past;
//   also we compute score as the latency extension caused by these swaps
//
// At the end, we have a list of Alters, each with a private Past, and a private latency extension.
// The partial paths represent lists of swaps to be inserted.
// The initial two-qubit gate gets the qubits at the ends of the partial paths as operands.
// The main selection criterium from the Alters is to select the one with the minimum latency extension.
// Having done that, the other Alters can be discarded and the selected one committed to the main Past.
class Alter
{
public:
    const ql::quantum_platform   *platformp;  // descriptions of resources for scheduling
    ql::quantum_kernel     *kernelp;    // kernel class pointer to allow calling kernel private methods
    size_t                  nq;         // width of Past and Virt2Real map is number of real qubits
    size_t                  ct;         // cycle time, multiplier from cycles to nano-seconds

    ql::gate               *targetgp;   // gate that this variation aims to make NN
    std::vector<size_t>     total;      // full path, including source and target nodes
    std::vector<size_t>     fromSource; // partial path after split, starting at source
    std::vector<size_t>     fromTarget; // partial path after split, starting at target, backward

    Past                    past;       // cloned main past, extended with swaps from this path
    double                  score;      // e.g. latency extension caused by the path
    bool                    didscore;   // initially false, true after assignment to score

// explicit Alter constructor
// needed for virgin construction
Alter()
{
    DOUT("Constructing Alter");
}

// Alter initializer
// This should only be called after a virgin construction and not after cloning a path.
void Init(const ql::quantum_platform* p, ql::quantum_kernel* k)
{
    DOUT("Alter::Init(number of qubits=" << nq);
    platformp = p;
    kernelp = k;

    nq = platformp->qubit_number;
    ct = platformp->cycle_time;
    // total, fromSource and fromTarget start as empty vectors
    past.Init(platformp, kernelp);      // initializes past to empty
    didscore = false;                   // will not print score for now
}

// printing facilities of Paths
// print path as hd followed by [0->1->2]
// and then followed by "implying" swap(q0,q1) swap(q1,q2)
void partialPrint(std::string hd, std::vector<size_t> & pp)
{
    if (!pp.empty())
    {
        int started = 0;
        for (auto & ppe : pp)
        {
            if (started == 0)
            {
                started = 1;
                std::cout << hd << "[";
            }
            else
            {
                std::cout << "->";
            }
            std::cout << ppe;
        }
        if (started == 1)
        {
            std::cout << "]";
//          if (pp.size() >= 2)
//          {
//              std::cout << " implying:";
//              for (size_t i = 0; i < pp.size()-1; i++)
//              {
//                  std::cout << " swap(q" << pp[i] << ",q" << pp[i+1] << ")";
//              }
//          }
        }
    }
}

void DPRINT(std::string s)
{
    if ( ql::utils::logger::LOG_LEVEL >= ql::utils::logger::log_level_t::LOG_DEBUG )
        Print(s);
}
void Print(std::string s)
{
    // std::cout << s << "- " << targetgp->qasm();
    std::cout << s << "- " << targetgp->qasm();
    if (fromSource.empty() && fromTarget.empty())
    {
        partialPrint(", total path:", total);
    }
    else
    {
        partialPrint(", path from source:", fromSource);
        partialPrint(", from target:", fromTarget);
    }
    if (didscore)
    {
        std::cout << ", score=" << score;
    }
    // past.Print("past in Alter");
    std::cout << std::endl;
}

static
void DPRINT(std::string s, std::vector<Alter> & va)
{
    if ( ql::utils::logger::LOG_LEVEL >= ql::utils::logger::log_level_t::LOG_DEBUG )
        Print(s, va);
}
static
void Print(std::string s, std::vector<Alter> & va)
{
    int started = 0;
    for (auto & a : va)
    {
        if (started == 0)
        {
            started = 1;
            std::cout << s << "[" << va.size() << "]={" << std::endl;
        }
        a.Print("");
    }
    if (started == 1)
    {
        std::cout << "}" << std::endl;
    }
}

static
void DPRINT(std::string s, std::list<Alter> & la)
{
    if ( ql::utils::logger::LOG_LEVEL >= ql::utils::logger::log_level_t::LOG_DEBUG )
        Print(s, la);
}
static
void Print(std::string s, std::list<Alter> & la)
{
    int started = 0;
    for (auto & a : la)
    {
        if (started == 0)
        {
            started = 1;
            std::cout << s << "[" << la.size() << "]={" << std::endl;
        }
        a.Print("");
    }
    if (started == 1)
    {
        std::cout << "}" << std::endl;
    }
}

// add a node to the path in front, extending its length with one
void Add2Front(size_t q)
{
    total.insert(total.begin(), q); // hopelessly inefficient
}

// add to a max of maxnumbertoadd swap gates for the current path to the given past
// this past can be a path-local one or the main past
// after having added them, schedule the result into that past
void AddSwaps(Past & past, std::string mapselectswapsopt)
{
    // DOUT("Addswaps " << mapselectswapsopt);
    if ("one"==mapselectswapsopt || "all"==mapselectswapsopt)
    {
        size_t  numberadded = 0;
        size_t  maxnumbertoadd = ("one"==mapselectswapsopt ? 1 : MAX_CYCLE);

        size_t  fromSourceQ;
        size_t  toSourceQ;
        fromSourceQ = fromSource[0];
        for ( size_t i = 1; i < fromSource.size() && numberadded < maxnumbertoadd; i++ )
        {
            toSourceQ = fromSource[i];
            past.AddSwap(fromSourceQ, toSourceQ);
            fromSourceQ = toSourceQ;
            numberadded++;
        }
    
        size_t  fromTargetQ;
        size_t  toTargetQ;
        fromTargetQ = fromTarget[0];
        for ( size_t i = 1; i < fromTarget.size() && numberadded < maxnumbertoadd; i++ )
        {
            toTargetQ = fromTarget[i];
            past.AddSwap(fromTargetQ, toTargetQ);
            fromTargetQ = toTargetQ;
            numberadded++;
        }
    }
    else
    {
        MapperAssert("earliest"==mapselectswapsopt);
        if (fromSource.size() >= 2 && fromTarget.size() >= 2)
        {
            if (past.IsFirstSwapEarliest(fromSource[0], fromSource[1], fromTarget[0], fromTarget[1]))
            {
                past.AddSwap(fromSource[0], fromSource[1]);
            }
            else
            {
                past.AddSwap(fromTarget[0], fromTarget[1]);
            }
        }
        else if (fromSource.size() >= 2)
        {
            past.AddSwap(fromSource[0], fromSource[1]);
        }
        else if (fromTarget.size() >= 2)
        {
            past.AddSwap(fromTarget[0], fromTarget[1]);
        }
    }

    past.Schedule();
}

// compute cycle extension of the current alternative in prevPast relative to the given base past
//
// Extend can be called in a deep exploration where pasts have been extended
// each one on top of a previous one, starting from the base past;
// the currPast here is the last extended one, i.e. on top of which this extension should be done;
// the basePast is the ultimate base past relative to which the total extension is to be computed.
//
// Do this by adding the swaps described by this alternative
// to an alternative-local copy of the current past;
// keep this resulting past in the current alternative (for later use);
// compute the total extension of all pasts relative to the base past
// and store this extension in the alternative's score for later use
void Extend(Past currPast, Past basePast)
{
    // DOUT("... clone past, add swaps, compute overall score and keep it all in current alternative");
    past = currPast;   // explicitly clone currPast to an alternative-local copy of it, Alter.past
    // DOUT("... adding swaps to alternative-local past ...");
    AddSwaps(past, "all");
    // DOUT("... done adding/scheduling swaps to alternative-local past");

    auto mapperopt = ql::options::get("mapper");
    if ("maxfidelity" == mapperopt)
    {
        FATAL("Mapper option maxfidelity has been disabled");
        // score = ql::quick_fidelity(past.lg);
    }
    else
    {
        score = past.MaxFreeCycle() - basePast.MaxFreeCycle();
    }
    didscore = true;
}

// split the path
// starting from the representation in the total attribute,
// generate all split path variations where each path is split once at any hop in it
// the intention is that the mapped two-qubit gate can be placed at the position of that hop
// all result paths are added/appended to the given result list
//
// when at the hop of a split a two-qubit gate cannot be placed, the split is not done there
// this means at the end that, when all hops are inter-core, no split is added to the result
//
// distance=5   means length=6  means 4 swaps + 1 CZ gate, e.g.
// index in total:      0           1           2           length-3        length-2        length-1
// qubit:               2   ->      5   ->      7   ->      3       ->      1       CZ      4
void Split(Grid &grid, std::list<Alter> & resla)
{
    // DOUT("Split ...");

    size_t length = total.size();
    MapperAssert (length >= 2);   // distance >= 1 so path at least: source -> target
    for (size_t rightopi = length-1; rightopi >= 1; rightopi--)
    {
        size_t leftopi = rightopi - 1;
        MapperAssert (leftopi >= 0);
        // DOUT("... leftopi=" << leftopi);
        // leftopi is the index in total that holds the qubit that becomes the left operand of the gate
        // rightopi is the index in total that holds the qubit that becomes the right operand of the gate
        // rightopi == leftopi + 1
        if (grid.IsInterCoreHop(total[leftopi], total[rightopi]))
        {
            // an inter-core hop cannot execute a two-qubit gate, so is not a valid alternative
            // DOUT("... skip inter-core hop from qubit=" << total[leftopi] << " to qubit=" << total[rightopi]);
            continue;
        }

        Alter    na = *this;      // na is local copy of the current path, including total
        // na = *this;            // na is local copy of the current path, including total
        // na.DPRINT("... copy of current alter");

        // fromSource will contain the path with qubits at indices 0 to leftopi
        // fromTarget will contain the path with qubits at indices rightopi to length-1, reversed
        //      reversal of fromTarget is done since swaps need to be generated starting at the target
        size_t fromi, toi;

        na.fromSource.resize(leftopi+1);
        // DOUT("... fromSource size=" << na.fromSource.size());
        for (fromi = 0, toi = 0; fromi <= leftopi; fromi++, toi++)
        {
            // DOUT("... fromSource: fromi=" << fromi << " toi=" << toi);
            na.fromSource[toi] = na.total[fromi];
        }

        na.fromTarget.resize(length-leftopi-1);
        // DOUT("... fromTarget size=" << na.fromTarget.size());
        for (fromi = length-1, toi = 0; fromi > leftopi; fromi--, toi++)
        {
            // DOUT("... fromTarget: fromi=" << fromi << " toi=" << toi);
            na.fromTarget[toi] = na.total[fromi];
        }

        // na.DPRINT("... copy of alter after split");
        resla.push_back(na);
        // DOUT("... added to result list");
        // DPRINT("... current alter after split");
    }
}

};  // end class Alter


// =========================================================================================
// Future: input window for mapper
//
// The future window shows the gates that still must be mapped as the availability list
// of a list scheduler that would work on a dependence graph representation of each input circuit.
// This future window is initialized once for the whole program, and gets a method call
// when it should switch to a new circuit (corresponding to a new kernel).
// In each circuit and thus each dependence graph the gates (including classical instruction) are found;
// the dependence graph models their dependences and also whether they act as barriers,
// an example of the latter being a classical branch.
// The availability list with gates (including classical instructions) is the main interface
// to the mapper, i.e. the mapper selects one or more element(s) from it to map next;
// it may even create alternatives for each combination of available gates.
// The gates in the list have attributes like criticality, which can be exploited by the mapper.
// The dependence graph and the availability list operations are provided by the Scheduler class.
//
// The future is a window because in principle it could be implemented incrementally,
// i.e. that the dependence graph would be extended when an attribute gets below a threshold,
// e.g. when successors of a gate are interrogated for a particular attribute.
// A problem might be that criticality requires having seen the end of the circuit,
// but the space overhead of this attribute is much less than that of a full dependence graph.
// The implementation below is not incremental: it creates the dep graph for a circuit completely.
//
// The implementation below just selects the most critical gate from the availability list
// as next candidate to map, the idea being that any collateral damage of mapping this gate
// will have a lower probability of increasing circuit depth
// than taking a non-critical gate as first one to map.
// Later implementations may become more sophisticated.
//
// With option maplookaheadopt=="no", the future window's dependence graph (scheduled and avlist) are not used.
// Instead a copy of the input circuit (input_gatepv) is created and iterated over (input_gatepp).

class Future
{
public:
    const ql::quantum_platform            *platformp;
    Scheduler                       *schedp;        // a pointer, since dependence graph doesn't change
    ql::circuit                     input_gatepv;   // input circuit when not using scheduler based avlist

    std::map<ql::gate*,bool>        scheduled;      // state: has gate been scheduled, here: done from future?
    std::list<ListDigraph::Node>    avlist;         // state: which nodes/gates are available for mapping now?
    ql::circuit::iterator           input_gatepp;   // state: alternative iterator in input_gatepv

// just program wide initialization
void Init( const ql::quantum_platform *p)
{
    // DOUT("Future::Init ...");
    platformp = p;
    // DOUT("Future::Init [DONE]");
}

// Set/switch input to the provided circuit
// nq and nc are parameters because nc may not be provided by platform but by kernel
// the latter should be updated when mapping multiple kernels
void SetCircuit(ql::quantum_kernel& kernel, Scheduler& sched, size_t nq, size_t nc)
{
    DOUT("Future::SetCircuit ...");
    schedp = &sched;
    std::string maplookaheadopt = ql::options::get("maplookahead");
    if ("no" == maplookaheadopt)
    {
        input_gatepv = kernel.c;                                // copy to free original circuit to allow outputing to
        input_gatepp = input_gatepv.begin();                    // iterator set to start of input circuit copy
    }
    else
    {
        schedp->init(kernel.c, *platformp, nq, nc);             // fills schedp->graph (dependence graph) from all of circuit
                                                                // and so also the original circuit can be output to after this
        for( auto & gp : kernel.c )
        {
            scheduled[gp] = false;   // none were scheduled
        }
        scheduled[schedp->instruction[schedp->s]] = false;      // also the dummy nodes not
        scheduled[schedp->instruction[schedp->t]] = false;
        avlist.clear();
        avlist.push_back(schedp->s);
        schedp->set_remaining(ql::forward_scheduling);          // to know criticality

        if (ql::options::get("print_dot_graphs") == "yes")
        {
            std::string     map_dot;
            std::stringstream fname;

            schedp->get_dot(map_dot);

            fname << ql::options::get("output_dir") << "/" << kernel.name << "_" << "mapper" << ".dot";
            IOUT("writing " << "mapper" << " dependence graph dot file to '" << fname.str() << "' ...");
            ql::utils::write_file(fname.str(), map_dot);
        }
    }
    DOUT("Future::SetCircuit [DONE]");
}

// Get from avlist all gates that are non-quantum into nonqlg
// Non-quantum gates include: classical, and dummy (SOURCE/SINK)
// Return whether some non-quantum gate was found
bool GetNonQuantumGates(std::list<ql::gate*>& nonqlg)
{
    nonqlg.clear();
    std::string maplookaheadopt = ql::options::get("maplookahead");
    if ("no" == maplookaheadopt)
    {
        ql::gate*   gp = *input_gatepp;
        if (input_gatepp != input_gatepv.end())
        {
            if (gp->type() == ql::__classical_gate__
                || gp->type() == ql::__dummy_gate__
                )
            {
                nonqlg.push_back(gp);
            }
        }
    }
    else
    {
        for ( auto n : avlist)
        {
            ql::gate*  gp = schedp->instruction[n];
            if (gp->type() == ql::__classical_gate__
                || gp->type() == ql::__dummy_gate__
                )
            {
                nonqlg.push_back(gp);
            }
        }
    }
    return nonqlg.size() != 0;
}

// Get all gates from avlist into qlg
// Return whether some gate was found
bool GetGates(std::list<ql::gate*>& qlg)
{
    qlg.clear();
    std::string maplookaheadopt = ql::options::get("maplookahead");
    if ("no" == maplookaheadopt)
    {
        if (input_gatepp != input_gatepv.end())
        {
            ql::gate*  gp = *input_gatepp;
            if (gp->operands.size() > 2)
            {
                FATAL(" gate: " << gp->qasm() << " has more than 2 operand qubits; please decompose such gates first before mapping.");
            }
            qlg.push_back(gp);
        }
    }
    else
    {
        for ( auto n : avlist)
        {
            ql::gate*  gp = schedp->instruction[n];
            if (gp->operands.size() > 2)
            {
                FATAL(" gate: " << gp->qasm() << " has more than 2 operand qubits; please decompose such gates first before mapping.");
            }
            qlg.push_back(gp);
        }
    }
    return qlg.size() != 0;
}

// Indicate that a gate currently in avlist has been mapped, can be taken out of the avlist
// and its successors can be made available
void DoneGate(ql::gate* gp)
{
    std::string maplookaheadopt = ql::options::get("maplookahead");
    if ("no" == maplookaheadopt)
    {
        input_gatepp = std::next(input_gatepp);
    }
    else
    {
        schedp->TakeAvailable(schedp->node[gp], avlist, scheduled, ql::forward_scheduling);
    }
}

// Return gp in lag that is most critical (provided lookahead is enabled)
// This is used in tiebreak, when every other option has failed to make a distinction.
ql::gate* MostCriticalIn(std::list<ql::gate*>& lag)
{
    std::string maplookaheadopt = ql::options::get("maplookahead");
    if ("no" == maplookaheadopt)
    {
        return lag.front();
    }
    else
    {
        return schedp->find_mostcritical(lag);
    }
}

};  // end class Future



// =========================================================================================
// Mapper: map operands of gates and insert swaps so that two-qubit gate operands are NN.
// All gates must be unary or two-qubit gates. The operands are virtual qubit indices.
// After mapping, all virtual qubit operands have been mapped to real qubit operands.

// For the mapper to work,
// the number of virtual qubits (nvq) must be less equal to the number of real qubits (nrq): nvq <= nrq;
// the mapper assumes that the virtual qubit operands (vqi) are encoded as a number 0 <= vqi < nvq
// and that the real qubit operands (rqi) are encoded as a number 0 <= rqi < nrq.
// The nrq is given by the platform, nvq is given by the program.
// The mapper ignores the latter (0 <= vqi < nvq was tested when creating the gates),
// and assumes vqi, nvq, rqi and nrq to be of the same type (size_t) 0<=qi<nrq.
// Because of this, it makes no difference between nvq and nrq, and refers to both as nq,
// and initializes the latter from the platform.
// All maps mapping virtual and real qubits to something are of size nq.

// Classical registers are ignored by the mapper currently. TO BE DONE.

// The mapping is done in the context of a grid of qubits defined by the given platform.
// This grid is initialized once for the whole program and constant after that.

// Each kernel in the program is independently mapped (see the Map method),
// ignoring inter-kernel control flow and thereby the requirement to pass on the current mapping.
// However, for each kernel there are two methods: initial placement and a heuristic,
// of which initial placement may do a half-hearted job, while heuristic will always be successful in finding a map;
// but what initial placement may find, it will be used by the heuristic as an initial mapping; they are in this order.

// Anticipating on the inter-kernel mapping, the mapper maintains a kernel input mapping coming from the context,
// and produces a kernel output mapping for the context; the mapper updates the kernel's circuit from virtual to real.
//
// Without inter-kernel control flow, the flow is as follows:
// - mapping starts from a 1 to 1 mapping of virtual to real qubits (the kernel input mapping)
//      in which all virtual qubits are initialized to a fixed constant state (|0>/inited), suitable for replacing swap by move
// - optionally attempt an initial placement of the circuit, starting from the kernel input mapping
//      and thus optionally updating the virtual to real map and the state of used virtuals (from inited to inuse)
// - anyhow use heuristics to map the input (or what initial placement left to do),
//      mapping the virtual gates to (sets of) real gates, and outputing the new map and the new virtuals' state
// - optionally decompose swap and/or cnot gates in the real circuit to primitives (MakePrimitive)

// Inter-kernel control flow and consequent mapping dependence between kernels is not implemented. TO BE DONE
// The design of mapping multiple kernels is as follows (TO BE ADAPTED TO NEW REALSTATE):
// The mapping is done kernel by kernel, in the order that they appear in the list of kernels:
// - initially the program wide initial mapping is a 1 to 1 mapping of virtual to real qubits
// - when start to map a kernel, there is a set of already mapped kernels, and a set of not yet mapped kernels;
//       of each mapped kernel, there is an output mapping, i.e. the mapping of virts to reals with the rs per virtual;
//       when mapping was ready, and the current kernel has a set of kernels
//       which are direct predecessor in the program's control flow;
//       a subset of those direct predecessors thus has been mapped and another subset not mapped;
//       the output mappings of the mapped predecessor kernels are input
// - unify these multiple input mappings to a single one; this may introduce swaps on the control flow edges;
//      the result is the input mapping of the current kernel; keep it for later reference
// - attempt an initial placement of the circuit, starting from the kernel input mapping
// - anyhow use heuristics to map the input (or what initial placement left to do)
// - when done:
//       keep the output mapping as the kernel's output mapping;
//       for all mapped successor kernels, compute a transition from output to their input,
//       and add it to the edge; the edge code must be optimized for:
//       - being empty: nothing needs to be done
//       - having a source with one succ; the edge code can be appended to that succ
//       - having a target with one pred; the edge code can be prepended to that pred
//       - otherwise, a separate intermediate kernel for the transition code must be created, and added
// THE ABOVE INTER-KERNEL MAPPING IS NOT IMPLEMENTED.

// The Mapper's main entry is Map which manages the input and output streams of QASM instructions,
// and does the logic between (global) initial placement mapper and the (more local) heuristic mapper.
// It selects the quantum gates from it, and maps these in the context of what was mapped before (the Past).
// Each gate is separately mapped in MapGate in the main Past's context.
class Mapper
{
private:
                                    // Initialized by Mapper::Init
                                    // OpenQL wide configuration, all constant after initialization
    const ql::quantum_platform *platformp;// current platform: topology and gate definitions
    ql::quantum_kernel   *kernelp;  // (copy of) current kernel (class) with free private circuit and methods
                                    // primarily to create gates in Past; Past is part of Mapper and of each Alter
                                    
    size_t          nq;             // number of qubits in the platform, number of real qubits
    size_t          nc;             // number of cregs in the platform, number of classical registers
    size_t          cycle_time;     // length in ns of a single cycle of the platform
                                    // is divisor of duration in ns to convert it to cycles
    Grid            grid;           // current grid

                                    // Initialized by Mapper.Map
    std::mt19937    gen;            // Standard mersenne_twister_engine, not yet seeded

public:
                                    // Passed back by Mapper::Map to caller for reporting
    size_t          nswapsadded;    // number of swaps added (including moves)
    size_t          nmovesadded;    // number of moves added
    std::vector<size_t> v2r_in;     // v2r[virtual qubit index] -> real qubit index | UNDEFINED_QUBIT
    std::vector<int>    rs_in;      // rs[real qubit index] -> {nostate|wasinited|hasstate}
    std::vector<size_t> v2r_ip;     // v2r[virtual qubit index] -> real qubit index | UNDEFINED_QUBIT
    std::vector<int>    rs_ip;      // rs[real qubit index] -> {nostate|wasinited|hasstate}
    std::vector<size_t> v2r_out;    // v2r[virtual qubit index] -> real qubit index | UNDEFINED_QUBIT
    std::vector<int>    rs_out;     // rs[real qubit index] -> {nostate|wasinited|hasstate}


// Mapper constructor is default synthesized

private:

// initial path finder
// generate paths with source src and target tgt as a list of path into resla;
// this result list resla is allocated by caller and is empty on the call;
// which indicates which paths are generated; see below the enum whichpaths;
// on top of this, the other mapper options apply
typedef
enum {
    wp_all_shortest,            // all shortest paths
    wp_left_shortest,           // only the shortest along the left side of the rectangle of src and tgt
    wp_right_shortest,          // only the shortest along the right side of the rectangle of src and tgt
    wp_leftright_shortest       // both the left and right shortest
} whichpaths_t;

// Find shortest paths between src and tgt in the grid, bounded by a particular strategy (which);
// budget is the maximum number of hops allowed in the path from src and is at least distance to tgt;
// it can be higher when not all hops qualify for doing a two-qubit gate or to find more than just the shortest paths.
void GenShortestPaths(ql::gate* gp, size_t src, size_t tgt, size_t budget, std::list<Alter> & resla, whichpaths_t which)
{
    std::list<Alter> genla;    // list that will get the result of a recursive Gen call

    // DOUT("GenShortestPaths: " << "src=" << src << " tgt=" << tgt << " budget=" << budget << " which=" << which);
    MapperAssert (resla.empty());

    if (src == tgt) {
        // found target
        // create a virgin Alter and initialize it to become an empty path
        // add src to this path (so that it becomes a distance 0 path with one qubit, src)
        // and add the Alter to the result list 
        Alter  a;
        a.Init(platformp, kernelp);
        a.targetgp = gp;
        a.Add2Front(src);
        resla.push_back(a);
        // a.DPRINT("... empty path after adding to result list");
        // Alter::DPRINT("... result list after adding empty path", resla);
        // DOUT("... will return now");
        return;
    }

    // start looking around at neighbors for serious paths
    size_t d = grid.Distance(src, tgt);
    MapperAssert (d >= 1);

    // reduce neighbors nbs to those n continuing a path within budget
    // src=>tgt is distance d, budget>=d is allowed, attempt src->n=>tgt
    // src->n is one hop, budget from n is one less so distance(n,tgt) <= budget-1 (i.e. distance < budget)
    // when budget==d, this defaults to distance(n,tgt) <= d-1
    auto nbl = grid.nbs[src];
    nbl.remove_if( [this,budget,tgt](const size_t& n) { return grid.Distance(n,tgt) >= budget; } );

    // rotate neighbor list nbl such that largest difference between angles of adjacent elements is beyond back()
    // this makes only sense when there is an underlying xy grid; when not, which can only be wp_all_shortest
    grid.Normalize(src, nbl);
    // subset to those neighbors that continue in direction(s) we want
    if (which == wp_left_shortest)
    {
        nbl.remove_if( [nbl](const size_t& n) { return n != nbl.front(); } );
    }
    else if (which == wp_right_shortest)
    {
        nbl.remove_if( [nbl](const size_t& n) { return n != nbl.back(); } );
    }
    else if (which == wp_leftright_shortest)
    {
        nbl.remove_if( [nbl](const size_t& n) { return n != nbl.front() && n != nbl.back(); } );
    }

    // std::cout << "... before iterating, nbl: ";
    // for (auto dn : nbl) { std::cout << dn << " "; } std::cout << std::endl;

    // for all resulting neighbors, find all continuations of a shortest path
    for (auto & n : nbl)
    {
        whichpaths_t newwhich = which;
        // but for each neighbor only look in desired direction, if any
        if (which == wp_leftright_shortest && nbl.size() != 1)
        {
            // when looking both left and right still, and there is a choice now, split into left and right
            if (n == nbl.front())
            {
                newwhich = wp_left_shortest;
            }
            else
            {
                newwhich = wp_right_shortest;
            }
        }
        GenShortestPaths(gp, n, tgt, budget-1, genla, newwhich);  // get list of possible paths in budget-1 from n to tgt in genla
        resla.splice(resla.end(), genla);           // moves all of genla to resla; makes genla empty
    }
    // resla contains all paths starting from a neighbor of src, to tgt

    // add src to front of all to-be-returned paths from src's neighbors to tgt
    for (auto & a : resla)
    {
        // DOUT("... GenShortestPaths, about to add src=" << src << "in front of path");
        a.Add2Front(src);
    }
    // DOUT("... GenShortestPaths: returning from call of:" << "src=" << src << " tgt=" << tgt << " budget=" << budget << " which=" << which);
}

// Generate shortest paths in the grid for making gate gp NN, from qubit src to qubit tgt, with an alternative for each one
// - compute budget; usually it is distance but it can be higher such as for multi-core
// - reduce the number of paths depending on the mappathselect option
// - when not all shortest paths found are valid, take these out
// - paths are further split because each split may give rise to a separate alternative
//      a split is a hop where the two-qubit gate is assumed to be done;
//      and after splitting each alternative contains two lists,
//      one before and one after (reversed) the envisioned two-qubit gate;
//      all result alternatives are such that a two-qubit gate can be placed at the split
// End result is a list of alternatives (in resla) suitable for being evaluated for any routing metric.
void GenShortestPaths(ql::gate* gp, size_t src, size_t tgt, std::list<Alter> & resla)
{
    std::list<Alter> directla;  // list that will hold all not-yet-split Alters directly from src to tgt

    size_t budget = grid.MinHops(src, tgt);
    std::string mappathselectopt = ql::options::get("mappathselect");
    if ("all" == mappathselectopt)
    {
        GenShortestPaths(gp, src, tgt, budget, directla, wp_all_shortest);
    }
    else if ("borders" == mappathselectopt)
    {
        GenShortestPaths(gp, src, tgt, budget, directla, wp_leftright_shortest);
    }
    else
    {
        FATAL("Unknown value of mapppathselect option " << mappathselectopt);
    }

    // DOUT("about to split the paths");
    for (auto & a : directla)
    {
        a.Split(grid, resla);
    }
    // Alter::DPRINT("... after generating and splitting the paths", resla);
}

// Generate all possible variations of making gp NN, starting from given past (with its mappings),
// and return the found variations by appending them to the given list of Alters, la
void GenAltersGate(ql::gate* gp, std::list<Alter>& la, Past& past)
{
    auto&   q = gp->operands;
    MapperAssert (q.size() == 2);
    size_t  src = past.MapQubit(q[0]);  // interpret virtual operands in past's current map
    size_t  tgt = past.MapQubit(q[1]);
    size_t d = grid.Distance(src, tgt);     // and find distance between real counterparts
    DOUT("GenAltersGate: " << gp->qasm() << " in real (q" << src << ",q" << tgt << ") at distance=" << d );
    past.DFcPrint();

    GenShortestPaths(gp, src, tgt, la);// find shortest paths from src to tgt, and split these
    MapperAssert(la.size() != 0);
    // Alter::DPRINT("... after GenShortestPaths", la);
}

// Generate all possible variations of making gates in lg NN, starting from given past (with its mappings),
// and return the found variations by appending them to the given list of Alters, la
// Depending on maplookahead only take first (most critical) gate or take all gates.
void GenAlters(std::list<ql::gate*> lg, std::list<Alter>& la, Past& past)
{
    std::string maplookaheadopt = ql::options::get("maplookahead");
    if ("all" == maplookaheadopt)
    {
        // create alternatives for each gate in lg
        // DOUT("GenAlters, " << lg.size() << " 2q gates; create an alternative for each");
        for (auto gp : lg)
        {
            // gen alternatives for gp and add these to la
            // DOUT("GenAlters: create alternatives for: " << gp->qasm());
            GenAltersGate(gp, la, past);  // gen all possible variations to make gp NN, in current v2r mapping ("past")
        }
    }
    else
    {
        // only take the first gate in avlist, the most critical one, and generate alternatives for it
        ql::gate*  gp = lg.front();
        // DOUT("GenAlters, " << lg.size() << " 2q gates; take first: " << gp->qasm());
        GenAltersGate(gp, la, past);  // gen all possible variations to make gp NN, in current v2r mapping ("past")
    }
}

// start the random generator with a seed
// that is unique to the microsecond
void RandomInit()
{
    auto ts = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    // DOUT("Seeding random generator with " << ts );
    gen.seed(ts);
}

// if the maptiebreak option indicates so,
// generate a random int number in range 0..count-1 and use
// that to index in list of alternatives and to return that one,
// otherwise return a fixed one (front, back or first most critical one
Alter ChooseAlter(std::list<Alter>& la, Future& future)
{
    if (la.size() == 1)
    {
        return la.front();
    }

    std::string maptiebreakopt = ql::options::get("maptiebreak");
    if ("critical" == maptiebreakopt)
    {
        std::list<ql::gate*> lag;
        for (auto& a : la)
        {
            lag.push_back(a.targetgp);
        }
        ql::gate*   gp = future.MostCriticalIn(lag);
        MapperAssert(gp != NULL);
        for (auto& a : la)
        {
            if (a.targetgp == gp)
            {
                // DOUT(" ... took first alternative with most critical target gate");
                return a;
            }
        }
        return la.front();
    }
    if ("random" == maptiebreakopt)
    {
        Alter res;
        std::uniform_int_distribution<> dis(0, (la.size()-1));
        size_t choice = dis(gen);
        size_t i = 0;
        for (auto& a : la)
        {
            if (i == choice)
            {
                res = a;
                break;
            }
            i++;
        }
        // DOUT(" ... took random draw " << choice << " from 0.." << (la.size()-1));
        return res;
    }
    if ("last" == maptiebreakopt)
    {
        // DOUT(" ... took last " << " from 0.." << (la.size()-1));
        return la.back();
    }
    if ("first" == maptiebreakopt)
    {
        // DOUT(" ... took first " << " from 0.." << (la.size()-1));
        return la.front();
    }
    return la.front();  // to shut up gcc
}

// Map the gate/operands of a gate that has been routed or doesn't require routing
void MapRoutedGate(ql::gate* gp, Past& past)
{
    DOUT("MapRoutedGate on virtual: " << gp->qasm() );

    // MakeReal of this gate maps its qubit operands and optionally updates its gate name
    // when the gate name was updated, a new gate with that name is created;
    // when that new gate is a composite gate, it is immediately decomposed (by gate creation)
    // the resulting gate/expansion (anyhow a sequence of gates) is collected in circ
    ql::circuit circ;   // result of MakeReal
    past.MakeReal(gp, circ);        
    for (auto newgp : circ)
    {
        DOUT(" ... new mapped real gate, about to be added to past: " << newgp->qasm() );
        past.AddAndSchedule(newgp);
    }
}

// commit Alter resa
// generating swaps in past
// and taking it out of future when done with it
void CommitAlter(Alter& resa, Future& future, Past& past)
{
    ql::gate*  resgp = resa.targetgp;   // and the 2q target gate then in resgp
    resa.DPRINT("... CommitAlter, alternative to commit, will add swaps and then map target 2q gate");

    std::string mapselectswapsopt = ql::options::get("mapselectswaps");
    resa.AddSwaps(past, mapselectswapsopt);

    // when only some swaps were added, the resgp might not yet be NN, so recheck
    auto&   q = resgp->operands;
    if (grid.Distance(past.MapQubit(q[0]), past.MapQubit(q[1])) == 1)
    {
        // resgp is NN: so done with this 2q gate
        // DOUT("... CommitAlter, target 2q is NN, map it and done: " << resgp->qasm());
        MapRoutedGate(resgp, past);     // the 2q target gate is NN now and thus can be mapped
        future.DoneGate(resgp);         // and then taken out of future
    } else {
        // DOUT("... CommitAlter, target 2q is not NN yet, keep it: " << resgp->qasm());
    }
}

// Find gates in future.avlist that do not require routing, take them out and map them.
// Ultimately, no gates remain or only gates that require routing.
// Return false when no gates remain at all.
// Return true when any gates remain; those gates are returned in lg.
//
// Behavior depends on the value of option maplookahead and alsoNN2q parameter
// alsoNN2q is true:
//   maplookahead == "no":             while (next in circuit is nonq or 1q) map gate; return when it is 2q (maybe NN)
//                                     in this case, GetNonQuantumGates only returns a nonq when it is next in circuit
//              == "1qfirst":          while (nonq or 1q) map gate; return most critical 2q (maybe NN)
//              == "noroutingfirst":   while (nonq or 1q or 2qNN) map gate; return most critical 2q (nonNN)
//              == "all":              while (nonq or 1q or 2qNN) map gate; return all 2q (nonNN)
// alsoNN2q is false:
//   maplookahead == "no":             while (next in circuit is nonq or 1q) map gate; return when it is 2q (maybe NN)
//                                     in this case, GetNonQuantumGates only returns a nonq when it is next in circuit
//              == "1qfirst":          while (nonq or 1q) map gate; return most critical 2q (nonNN or NN)
//              == "noroutingfirst":   while (nonq or 1q) map gate; return most critical 2q (nonNN or NN)
//              == "all":              while (nonq or 1q) map gate; return all 2q (nonNN or NN)
//
bool MapMappableGates(Future& future, Past& past, std::list<ql::gate*>& lg, bool alsoNN2q)
{
    std::list<ql::gate*>   nonqlg; // list of non-quantum gates in avlist
    std::list<ql::gate*>   qlg;    // list of (remaining) gates in avlist

    DOUT("MapMappableGates entry");
    while (1)
    {
        if (future.GetNonQuantumGates(nonqlg))
        {
            // avlist contains non-quantum gates
            // and GetNonQuantumGates indicates these (in nonqlg) must be done first
            DOUT("MapMappableGates, there is a set of non-quantum gates");
            for (auto gp : nonqlg)
            {
                // here add code to map qubit use of any non-quantum instruction????
                // dummy gates are nonq gates internal to OpenQL such as SOURCE/SINK; don't output them
                if ( gp->type() != ql::__dummy_gate__)
                {
                    // past only can contain quantum gates, so non-quantum gates must by-pass Past
                    past.ByPass(gp);    // this flushes past.lg first to outlg
                }
                future.DoneGate(gp); // so on avlist= nonNN2q -> NN2q -> 1q -> nonq: the nonq is done first
                DOUT("MapMappableGates, done with " << gp->qasm());
            }
            DOUT("MapMappableGates, done with set of non-quantum gates, continuing ...");
            continue;
        }
        if (!future.GetGates(qlg))
        {
            DOUT("MapMappableGates, no gates anymore, return");
            // avlist doesn't contain any gate
            lg.clear();
            return false;
        }

        // avlist contains quantum gates
        // and GetNonQuantumGates/GetGates indicate these (in qlg) must be done now
        bool foundone = false;  // whether a quantum gate was found that never requires routing
        for (auto gp : qlg)
        {
            if ( gp->type() == ql::gate_type_t::__wait_gate__
                || gp->operands.size() == 1
                )
            {
                // a quantum gate not requiring routing ever is found
                MapRoutedGate(gp, past);
                future.DoneGate(gp);
                foundone = true;    // a quantum gate was found that never requires routing
                                    // so on avlist= nonNN2q -> NN2q -> 1q: the 1q is done first
                break;
            }
        }
        if (foundone)
        {
            continue;
        }
        // qlg only contains 2q gates (that could require routing)
        if ( alsoNN2q)
        {
            // when there is a 2q in qlg that is mappable already, map it
            // when more, take most critical one first (because qlg is ordered, most critical first)
            for (auto gp : qlg)
            {
                auto&   q = gp->operands;
                size_t  src = past.MapQubit(q[0]);      // interpret virtual operands in current map
                size_t  tgt = past.MapQubit(q[1]);
                size_t  d = grid.Distance(src, tgt);    // and find distance between real counterparts
                if (d == 1)
                {
                    DOUT("MapMappableGates, NN no routing: " << gp->qasm() << " in real (q" << src << ",q" << tgt << ")");
                    MapRoutedGate(gp, past);
                    future.DoneGate(gp);
                    foundone = true;    // a 2q quantum gate was found that was mappable
                                        // so on avlist= nonNN2q -> NN2q: the NN2q is done first
                    break;
                }
            }
            if (foundone)
            {
                // found a mappable 2q one, and mapped it
                // don't map more mappable 2q ones
                // (they might not be critical, the now available 1q gates may hide a more critical 2q gate),
                // but deal with all available non-quantum, and 1q gates first,
                // and only when non of those remain, map a next mappable 2q (now most critical) one
                DOUT("MapMappableGates, found and mapped an easy quantum gate, continuing ...");
                continue;
            }
            DOUT("MapMappableGates, only nonNN 2q gates remain: ...");
        }
        else
        {
            DOUT("MapMappableGates, only 2q gates remain (nonNN and NN): ...");
        }
        // avlist (qlg) only contains 2q gates (when alsoNN2q: only non-NN ones; otherwise also perhaps NN ones)
        lg = qlg;
        if ( ql::utils::logger::LOG_LEVEL >= ql::utils::logger::log_level_t::LOG_DEBUG )
        {
            for (auto gp: lg)
            {
                DOUT("... 2q gate returned: " << gp->qasm());
            }
        }
        return true;
    }
}

// select Alter determined by strategy defined by mapper options
// - if base[rc], select from whole list of Alters, of which all 'remain'
// - if minextend[rc], select Alter from list of Alters with minimal cycle extension of given past
//   when several remain with equal minimum extension, recurse to reduce this set of remaining ones
//   - level: level of recursion at which SelectAlter is called: 0 is base, 1 is 1st, etc.
//   - option mapselectmaxlevel: max level of recursion to use, where inf indicates no maximum
// - maptiebreak option indicates which one to take when several (still) remain
// result is returned in resa
void SelectAlter(std::list<Alter>& la, Alter & resa, Future& future, Past& past, Past& basePast, int level)
{
                                // la are all alternatives we enter with
    MapperAssert(!la.empty());  // so there is always a result Alter

    std::list<Alter> gla;       // good alternative subset of la, suitable to go in recursion with
    std::list<Alter> bla;       // best alternative subset of gla, suitable to choose result from

    DOUT("SelectAlter ENTRY level=" << level << " from " << la.size() << " alternatives");
    auto mapperopt = ql::options::get("mapper");
    if (mapperopt == "base"|| mapperopt == "baserc")
    {
        Alter::DPRINT("... SelectAlter base (equally good/best) alternatives:", la);
        resa = ChooseAlter(la, future);
        resa.DPRINT("... the selected Alter is");
        // DOUT("SelectAlter DONE level=" << level << " from " << la.size() << " alternatives");
        return;
    }
    MapperAssert(mapperopt == "minextend" || mapperopt == "minextendrc" || mapperopt == "maxfidelity");

    // Compute a.score of each alternative relative to basePast, and sort la on it, minimum first
    for (auto & a : la)
    {
        a.DPRINT("Considering extension by alternative: ...");
        a.Extend(past, basePast);           // locally here, past will be cloned and kept in alter
                                            // and the extension stored into the a.score
    }
    la.sort([this](const Alter &a1, const Alter &a2) { return a1.score < a2.score; });
    Alter::DPRINT("... SelectAlter sorted all entry alternatives after extension:", la);

    // Reduce sorted list of alternatives (la) to list of good alternatives (gla)
    // suitable to find in recursion which is/are really best;
    // this need not be only those with minimum extend, it can be more.
    // With option mapselectmaxwidth="min" in which "min" stands for minimal number, we get just those minimal ones.
    // With other option values, we are more forgiving but that easily lets the number of alternatives explode.
    gla = la;
    gla.remove_if( [this,la](const Alter& a) { return a.score != la.front().score; } );
    size_t  las = la.size();
    size_t  glas = gla.size();
    auto mapselectmaxwidthopt = ql::options::get("mapselectmaxwidth");
    if ("min" != mapselectmaxwidthopt)
    {
        size_t  keep = 1;
        if ( "minplusone" == mapselectmaxwidthopt)
        {
            keep = glas+1;
        }
        else if ("minplushalfmin" == mapselectmaxwidthopt)
        {
            keep = glas+glas/2;
        }
        else if ( "minplusmin" == mapselectmaxwidthopt)
        {
            keep = glas*2;
        }
        else if ( "all" == mapselectmaxwidthopt)
        {
            keep = las;
        }
        if (keep < las)
        {
            gla = la;
            std::list<Alter>::iterator  ia;
            ia = gla.begin();
            advance(ia, keep);
            gla.erase(ia, gla.end());
        }
        else
        {
            gla = la;
        }
    }
    // DOUT("SelectAlter mapselectmaxwidth=" << mapselectmaxwidthopt << " level=" << level << " reduced la to gla");
    Alter::DPRINT("... SelectAlter good alternatives before recursion:", gla);

    // Prepare for recursion;
    // option mapselectmaxlevel indicates the maximum level of recursion (0 is no recursion)
    auto mapselectmaxlevelstring = ql::options::get("mapselectmaxlevel");
    int  mapselectmaxlevel = ("inf" == mapselectmaxlevelstring) ?  MAX_CYCLE : atoi(mapselectmaxlevelstring.c_str());

    // When maxlevel has been reached, stop the recursion, and choose from the best minextend/maxfidelity alternatives
    if (level >= mapselectmaxlevel)
    {
        // Reduce list of good alternatives (gla) to list of minextend/maxfidelity best alternatives (bla)
        // and make a choice from that list to return as result
        bla = gla;
        bla.remove_if( [this,gla](const Alter& a) { return a.score != gla.front().score; } );
        Alter::DPRINT("... SelectAlter reduced to best alternatives to choose result from:", bla);
        resa = ChooseAlter(bla, future);
        resa.DPRINT("... the selected Alter (STOPPING RECURSION) is");
        // DOUT("SelectAlter DONE level=" << level << " from " << bla.size() << " best alternatives");
        return;
    }

    // Otherwise, prepare using recursion to choose from the good alternatives,
    // i.e. make a recursion step looking ahead to decide which alternative is best
    //
    // For each alternative in gla,
    // lookahead for next non-NN2q gates, and comparing them for their alternative mappings;
    // the lookahead alternative with the least overall extension (i.e. relative to basePast) is chosen,
    // and the current alternative on top of which it was build
    // is chosen at the current level, unwinding the recursion.
    //
    // Recursion could stop above because:
    // - max level of recursion was reached
    // Recursion can stop here because of:
    // - end-of-circuit (no non-NN 2q gates remain).
    //
    // When gla.size() == 1, we still want to know its minimum extension, to compare with competitors,
    // since that is not just a local figure but the extension from basePast;
    // so indeed with only one alternative we may still go into recursion below.
    // This means that recursion always goes to maxlevel or end-of-circuit.
    // This anomaly may need correction.
    // DOUT("... SelectAlter level=" << level << " entering recursion with " << gla.size() << " good alternatives");
    for (auto& a : gla)
    {
        a.DPRINT("... ... considering alternative:");
        Future future_copy = future;            // copy!
        Past   past_copy = past;                // copy!
        CommitAlter(a, future_copy, past_copy);
        a.DPRINT("... ... committed this alternative first before recursion:");

        bool    havegates;                  // are there still non-NN 2q gates to map?
        std::list<ql::gate*> lg;            // list of non-NN 2q gates taken from avlist, as returned from MapMappableGates
        std::string maplookaheadopt = ql::options::get("maplookahead");
        std::string maprecNN2qopt = ql::options::get("maprecNN2q");
        // In recursion, look at option maprecNN2q:
        // - MapMappableGates with alsoNN2q==true is greedy and immediately maps each 1q and NN 2q gate
        // - MapMappableGates with alsoNN2q==false is not greedy, maps all 1q gates but not the (NN) 2q gates
        //
        // when yes and when maplookaheadopt is noroutingfirst or all, let MapMappableGates stop mapping only on nonNN2q
        // when no, let MapMappableGates stop mapping on any 2q
        // This creates more clear recursion: one 2q at a time instead of a possible empty set of NN2qs followed by a nonNN2q;
        // also when a NN2q is found, this is perfect; this is not seen when immediately mapping all NN2qs.
        // So goal is to prove that maprecNN2q should be no at this place, in the recursion step, but not at level 0!
        bool alsoNN2q = ("yes" == maprecNN2qopt) && ( "noroutingfirst" == maplookaheadopt || "all" == maplookaheadopt );
        havegates = MapMappableGates(future_copy, past_copy, lg, alsoNN2q); // map all easy gates; remainder returned in lg

        if (havegates)
        {
            // DOUT("... ... SelectAlter level=" << level << ", committed + mapped easy gates, now facing " << lg.size() << " 2q gates to evaluate next");
            std::list<Alter> la;                // list that will hold all variations, as returned by GenAlters
            GenAlters(lg, la, past_copy);       // gen all possible variations to make gates in lg NN, in current past.v2r mapping
            // DOUT("... ... SelectAlter level=" << level << ", generated for these 2q gates " << la.size() << " alternatives; RECURSE ... ");
            Alter resa;                         // result alternative selected and returned by next SelectAlter call
            SelectAlter(la, resa, future_copy, past_copy, basePast, level+1); // recurse, best in resa ...
            resa.DPRINT("... ... SelectAlter, generated for these 2q gates ... ; RECURSE DONE; resulting alternative ");
            a.score = resa.score;               // extension of deep recursion is treated as extension at current level,
                                                // by this an alternative started bad may be compensated by deeper alts
        }
        else
        {
            // DOUT("... ... SelectAlter level=" << level << ", no gates to evaluate next; RECURSION BOTTOM");
            auto mapperopt = ql::options::get("mapper");
            if ("maxfidelity" == mapperopt)
            {
                FATAL("Mapper option maxfidelity has been disabled");
                // a.score = ql::quick_fidelity(past_copy.lg);
            }
            else
            {
                a.score = past_copy.MaxFreeCycle() - basePast.MaxFreeCycle();
            }
            a.DPRINT("... ... SelectAlter, after committing this alternative, mapped easy gates, no gates to evaluate next; RECURSION BOTTOM");
        }
        a.DPRINT("... ... DONE considering alternative:");
    }
    // Sort list of good alternatives (gla) on score resulting after recursion
    gla.sort([this](const Alter &a1, const Alter &a2) { return a1.score < a2.score; });
    Alter::DPRINT("... SelectAlter sorted alternatives after recursion:", gla);

    // Reduce list of good alternatives (gla) of before recursion to list of equally minimal best alternatives now (bla)
    // and make a choice from that list to return as result
    bla = gla;
    bla.remove_if( [this,gla](const Alter& a) { return a.score != gla.front().score; } );
    Alter::DPRINT("... SelectAlter equally best alternatives on return of RECURSION:", bla);
    resa = ChooseAlter(bla, future);
    resa.DPRINT("... the selected Alter is");
    // DOUT("... SelectAlter level=" << level << " selecting from " << bla.size() << " equally good alternatives above DONE");
    DOUT("SelectAlter DONE level=" << level << " from " << la.size() << " alternatives");
}

// Given the states of past and future
// map all mappable gates and find the non-mappable ones
// for those evaluate what to do next and do it;
// during recursion, comparison is done with the base past (bottom of recursion stack),
// and past is the last past (top of recursion stack) relative to which the mapping is done.
void MapGates(Future& future, Past& past, Past& basePast)
{
    std::list<ql::gate*>   lg;              // list of non-mappable gates taken from avlist, as returned from MapMappableGates
    std::string maplookaheadopt = ql::options::get("maplookahead");
    bool alsoNN2q = ( "noroutingfirst" == maplookaheadopt || "all" == maplookaheadopt );
    while (MapMappableGates(future, past, lg, alsoNN2q))  // returns false when no gates remain
    {
        // all gates in lg are two-qubit quantum gates that cannot be mapped
        // select which one(s) to (partially) route, according to one of the known strategies
        // the only requirement on the code below is that at least something is done that decreases the problem

        // generate all variations
        std::list<Alter> la;                // list that will hold all variations, as returned by GenAlters
        GenAlters(lg, la, past);            // gen all possible variations to make gates in lg NN, in current past.v2r mapping
    
        // select best one
        Alter resa;
        SelectAlter(la, resa, future, past, basePast, 0);
                                            // select one according to strategy specified by options; result in resa
    
        // commit to best one
        // add all or just one swap, as described by resa, to THIS past, and schedule them/it in
        CommitAlter(resa, future, past);
    }
}

// Map the circuit's gates in the provided context (v2r maps), updating circuit and v2r maps
void MapCircuit(ql::quantum_kernel& kernel, Virt2Real& v2r)
{
    Future  future;         // future window, presents input in avlist
    Past    mainPast;       // past window, contains output schedule, storing all gates until taken out
    Scheduler sched;        // new scheduler instance (from src/scheduler.h) used for its dependence graph

    future.Init(platformp);
    future.SetCircuit(kernel, sched, nq, nc); // constructs depgraph, initializes avlist, ready for producing gates
    kernel.c.clear();       // future has copied kernel.c to private data; kernel.c ready for use by new_gate
    kernelp = &kernel;      // keep kernel to call kernelp->gate() inside Past.new_gate(), to create new gates

    mainPast.Init(platformp, kernelp);  // mainPast and Past clones inside Alters ready for generating output schedules into
    mainPast.ImportV2r(v2r);    // give it the current mapping/state
    // mainPast.DPRINT("start mapping");

    MapGates(future, mainPast, mainPast);
    mainPast.FlushAll();                // all output to mainPast.outlg, the output window of mainPast

    // mainPast.DPRINT("end mapping");

    DOUT("... retrieving outCirc from mainPast.outlg; swapping outCirc with kernel.c, kernel.c contains output circuit");
    ql::circuit outCirc;
    mainPast.Out(outCirc);                          // copy (final part of) mainPast's output window into this outCirc
    kernel.c.swap(outCirc);                         // and then to kernel.c
    kernel.cycles_valid = true;                     // decomposition was scheduled in; see Past.Add() and Past.Schedule()
    mainPast.ExportV2r(v2r);
    nswapsadded = mainPast.NumberOfSwapsAdded();
    nmovesadded = mainPast.NumberOfMovesAdded();
}

public:

// decompose all gates that have a definition with _prim appended to its name
void MakePrimitives(ql::quantum_kernel& kernel)
{
    DOUT("MakePrimitives circuit ...");

    ql::circuit     input_gatepv = kernel.c;    // copy to allow kernel.c use by Past.new_gate
    kernel.c.clear();                           // kernel.c ready for use by new_gate

    Past            mainPast;                   // output window in which gates are scheduled
    mainPast.Init(platformp, kernelp);

    for( auto & gp : input_gatepv )
    {
        ql::circuit tmpCirc;
        mainPast.MakePrimitive(gp, tmpCirc);    // decompose gp into tmpCirc; on failure, copy gp into tmpCirc
        for (auto newgp : tmpCirc)
        {
            mainPast.AddAndSchedule(newgp);     // decomposition is scheduled in gate by gate
        }
    }
    mainPast.FlushAll();

    ql::circuit     outCirc;                    // ultimate output gate stream
    mainPast.Out(outCirc);
    kernel.c.swap(outCirc);
    kernel.cycles_valid = true;                 // decomposition was scheduled in above

    DOUT("MakePrimitives circuit [DONE]");
}   // end MakePrimitives

// map kernel's circuit, main mapper entry once per kernel
// JvS: moved to mapper.cc ahead of restructuring everything else for persistent INITIALPLACE switch
void Map(ql::quantum_kernel& kernel);

// initialize mapper for whole program
// lots could be split off for the whole program, once that is needed
//
// initialization for a particular kernel is separate (in Map entry)
void Init(const ql::quantum_platform *p)
{
    // DOUT("Mapping initialization ...");
    // DOUT("... Grid initialization: platform qubits->coordinates, ->neighbors, distance ...");
    platformp = p;
    nq = p->qubit_number;
    // nc = p->creg_number;  // nc should come from platform, but doesn't; is taken from kernel in Map
    RandomInit();
    // DOUT("... platform/real number of qubits=" << nq << ");
    cycle_time = p->cycle_time;

    grid.Init(platformp);

    // DOUT("Mapping initialization [DONE]");
}   // end Init


};  // end class Mapper

#endif // ifndef QL_MAPPER_H
