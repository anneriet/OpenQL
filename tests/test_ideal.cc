#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <cassert>

#include <time.h>

#include <ql/openql.h>

// test cnot control operand commutativity
// i.e. best result is the reverse original order
void
test_cnot_controlcommute(std::string v)
{
    int n = 7;
    std::string prog_name = "test_" + v;
    std::string kernel_name = "test_" + v;
    float sweep_points[] = { 1 };

    ql::quantum_platform starmon("starmon", "test_ideal.json");
    ql::set_platform(starmon);
    ql::quantum_program prog(prog_name, starmon, n, 0);
    ql::quantum_kernel k(kernel_name, starmon, n, 0);
    prog.set_sweep_points(sweep_points, sizeof(sweep_points)/sizeof(float));

    // for (int j=0; j<7; j++) { k.gate("x", j); }

    k.gate("cnot", 3,0);
    k.gate("cnot", 3,6);
    k.gate("t", 6);
    k.gate("y", 6);
    k.gate("cnot", 3,1);
    k.gate("t", 1);
    k.gate("y", 1);
    k.gate("t", 1);
    k.gate("y", 1);
    k.gate("cnot", 3,5);
    k.gate("t", 5);
    k.gate("y", 5);
    k.gate("t", 5);
    k.gate("y", 5);
    k.gate("t", 5);
    k.gate("y", 5);

    // for (int j=0; j<7; j++) { k.gate("x", j); }

    prog.add(k);

    prog.compile( );
}

// test cnot target operand commutativity
// i.e. best result is the reverse original order
void
test_cnot_targetcommute(std::string v)
{
    int n = 7;
    std::string prog_name = "test_" + v;
    std::string kernel_name = "test_" + v;
    float sweep_points[] = { 1 };

    ql::quantum_platform starmon("starmon", "test_ideal.json");
    ql::set_platform(starmon);
    ql::quantum_program prog(prog_name, starmon, n, 0);
    ql::quantum_kernel k(kernel_name, starmon, n, 0);
    prog.set_sweep_points(sweep_points, sizeof(sweep_points)/sizeof(float));

    // for (int j=0; j<7; j++) { k.gate("x", j); }

    k.gate("cnot", 0,3);
    k.gate("cnot", 6,3);
    k.gate("t", 6);
    k.gate("y", 6);
    k.gate("cnot", 1,3);
    k.gate("t", 1);
    k.gate("y", 1);
    k.gate("t", 1);
    k.gate("y", 1);
    k.gate("cnot", 5,3);
    k.gate("t", 5);
    k.gate("y", 5);
    k.gate("t", 5);
    k.gate("y", 5);
    k.gate("t", 5);
    k.gate("y", 5);

    // for (int j=0; j<7; j++) { k.gate("x", j); }

    prog.add(k);

    prog.compile( );
}

// test cz any operand commutativity
// i.e. best result is the reverse original order
void
test_cz_anycommute(std::string v)
{
    int n = 7;
    std::string prog_name = "test_" + v;
    std::string kernel_name = "test_" + v;
    float sweep_points[] = { 1 };

    ql::quantum_platform starmon("starmon", "test_ideal.json");
    ql::set_platform(starmon);
    ql::quantum_program prog(prog_name, starmon, n, 0);
    ql::quantum_kernel k(kernel_name, starmon, n, 0);
    prog.set_sweep_points(sweep_points, sizeof(sweep_points)/sizeof(float));

    // for (int j=0; j<7; j++) { k.gate("x", j); }

    k.gate("cz", 0,3);
    k.gate("cz", 3,6);
    k.gate("t", 6);
    k.gate("y", 6);
    k.gate("cz", 1,3);
    k.gate("t", 1);
    k.gate("y", 1);
    k.gate("t", 1);
    k.gate("y", 1);
    k.gate("cz", 3,5);
    k.gate("t", 5);
    k.gate("y", 5);
    k.gate("t", 5);
    k.gate("y", 5);
    k.gate("t", 5);
    k.gate("y", 5);

    // for (int j=0; j<7; j++) { k.gate("x", j); }

    prog.add(k);

    prog.compile( );
}

// steane qec on s7 with cnots
void
test_steaneqec(std::string v)
{
    int n = 7;
    std::string prog_name = "test_" + v;
    std::string kernel_name = "test_" + v;
    float sweep_points[] = { 1 };

    ql::quantum_platform starmon("starmon", "test_ideal.json");
    ql::set_platform(starmon);
    ql::quantum_program prog(prog_name, starmon, n, 0);
    ql::quantum_kernel k(kernel_name, starmon, n, 0);
    prog.set_sweep_points(sweep_points, sizeof(sweep_points)/sizeof(float));

    k.gate("prepz", 3);
    k.gate("prepz", 5);
    k.gate("h", 5);
    k.gate("cnot", 5, 3);
    k.gate("cnot", 0, 3);
    k.gate("cnot", 1, 3);
    k.gate("cnot", 6, 3);
    k.gate("cnot", 2, 5);
    k.gate("cnot", 5, 3);
    k.gate("h", 5);
    k.gate("measure", 3);
    k.gate("measure", 5);

    prog.add(k);
    
    prog.compile( );
}

// all cnots with operands that are neighbors in s7
// no or hardly any significant difference between pre179 and post179 scheduling,
// slight differences may occur when the json file maps cnot to its constituent primitive gates
void
test_manyNN(std::string v)
{
    int n = 7;
    std::string prog_name = "test_" + v;
    std::string kernel_name = "test_" + v;
    float sweep_points[] = { 1 };

    ql::quantum_platform starmon("starmon", "test_ideal.json");
    ql::set_platform(starmon);
    ql::quantum_program prog(prog_name, starmon, n, 0);
    ql::quantum_kernel k(kernel_name, starmon, n, 0);
    prog.set_sweep_points(sweep_points, sizeof(sweep_points)/sizeof(float));

    for (int j=0; j<7; j++) { k.gate("x", j); }

    // a list of all cnots that are ok in trivial mapping
    k.gate("cnot", 0,2);
    k.gate("cnot", 0,3);
    k.gate("cnot", 1,3);
    k.gate("cnot", 1,4);
    k.gate("cnot", 2,0);
    k.gate("cnot", 2,5);
    k.gate("cnot", 3,0);
    k.gate("cnot", 3,1);
    k.gate("cnot", 3,5);
    k.gate("cnot", 3,6);
    k.gate("cnot", 4,1);
    k.gate("cnot", 4,6);
    k.gate("cnot", 5,2);
    k.gate("cnot", 5,3);
    k.gate("cnot", 6,3);
    k.gate("cnot", 6,4);

    for (int j=0; j<7; j++) { k.gate("x", j); }

    prog.add(k);

    prog.compile( );
}

int main(int argc, char ** argv)
{
    ql::utils::logger::set_log_level("LOG_DEBUG");
    ql::options::set("scheduler_uniform", "no");
    ql::options::set("scheduler", "ALAP");
    ql::options::set("scheduler_post179", "yes");

    test_cnot_controlcommute("cnot_controlcommute");
    test_cnot_targetcommute("cnot_targetcommute");
    test_cz_anycommute("cz_anycommute");
    test_steaneqec("steaneqec");
    test_manyNN("manyNN");

    return 0;
}