#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <cassert>

#include <time.h>

#include <ql/openql.h>


// test wait as gate
void
test_wait_simple(std::string v, std::string schedopt, std::string sched_post179opt)
{
    int n = 7;
    std::string prog_name = "test_" + v + "_schedopt=" + schedopt + "_sched_post179opt=" + sched_post179opt;
    std::string kernel_name = "test_" + v + "_schedopt=" + schedopt + "_sched_post179opt=" + sched_post179opt;
    float sweep_points[] = { 1 };

    ql::quantum_platform starmon("starmon", "hardware_config_cc_light.json");
    ql::set_platform(starmon);
    ql::quantum_program prog(prog_name, starmon, n, 0);
    ql::quantum_kernel k(kernel_name, starmon, n, 0);
    prog.set_sweep_points(sweep_points, sizeof(sweep_points)/sizeof(float));

    std::vector<size_t> operands = {0};

    k.gate("x", 0);
    k.wait(operands, 40);
    k.gate("x", 0);

    prog.add(k);

    ql::options::set("scheduler", schedopt);
    ql::options::set("scheduler_post179", sched_post179opt);
    prog.compile( );
}

void
test_wait_parallel(std::string v, std::string schedopt, std::string sched_post179opt)
{
    int n = 7;
    std::string prog_name = "test_" + v + "_schedopt=" + schedopt + "_sched_post179opt=" + sched_post179opt;
    std::string kernel_name = "test_" + v + "_schedopt=" + schedopt + "_sched_post179opt=" + sched_post179opt;
    float sweep_points[] = { 1 };

    ql::quantum_platform starmon("starmon", "hardware_config_cc_light.json");
    ql::set_platform(starmon);
    ql::quantum_program prog(prog_name, starmon, n, 0);
    ql::quantum_kernel k(kernel_name, starmon, n, 0);
    prog.set_sweep_points(sweep_points, sizeof(sweep_points)/sizeof(float));

    std::vector<size_t> operands = {1};

    k.gate("x", 0);
    k.wait(operands, 20);
    k.gate("x", 1);

    prog.add(k);

    ql::options::set("scheduler", schedopt);
    ql::options::set("scheduler_post179", sched_post179opt);
    prog.compile( );
}

void
test_wait_barrier(std::string v, std::string schedopt, std::string sched_post179opt)
{
    int n = 7;
    std::string prog_name = "test_" + v + "_schedopt=" + schedopt + "_sched_post179opt=" + sched_post179opt;
    std::string kernel_name = "test_" + v + "_schedopt=" + schedopt + "_sched_post179opt=" + sched_post179opt;
    float sweep_points[] = { 1 };

    ql::quantum_platform starmon("starmon", "hardware_config_cc_light.json");
    ql::set_platform(starmon);
    ql::quantum_program prog(prog_name, starmon, n, 0);
    ql::quantum_kernel k(kernel_name, starmon, n, 0);
    prog.set_sweep_points(sweep_points, sizeof(sweep_points)/sizeof(float));

    std::vector<size_t> operands = {0,1};

    k.gate("x", 0);
    k.gate("x", 1);
    k.gate("y", 0);
    k.wait(operands, 0);
    k.gate("measure", 0);
    k.gate("measure", 1);

    prog.add(k);

    ql::options::set("scheduler", schedopt);
    ql::options::set("scheduler_post179", sched_post179opt);
    prog.compile( );
}

void
test_wait_wouter(std::string v, std::string schedopt, std::string sched_post179opt)
{
    int n = 7;
    std::string prog_name = "test_" + v + "_schedopt=" + schedopt + "_sched_post179opt=" + sched_post179opt;
    std::string kernel_name = "test_" + v + "_schedopt=" + schedopt + "_sched_post179opt=" + sched_post179opt;
    float sweep_points[] = { 1 };

    ql::quantum_platform starmon("starmon", "hardware_config_cc_light.json");
    ql::set_platform(starmon);
    ql::quantum_program prog(prog_name, starmon, n, n);
    ql::quantum_kernel k(kernel_name, starmon, n, n);
    prog.set_sweep_points(sweep_points, sizeof(sweep_points)/sizeof(float));

	k.gate("ry90", 0);
	k.gate("ry90", 1);
	k.gate("ry90", 2);
	k.gate("ry90", 3);
	k.gate("ry90", 4);
	k.wait({0, 1, 2, 3, 4}, 0);
	
	k.gate("measure", std::vector<size_t> {0}, std::vector<size_t> {0});
	k.wait({0}, 0);
	
	k.gate("ym90", 5);

    prog.add(k);

    ql::options::set("scheduler", schedopt);
    ql::options::set("scheduler_post179", sched_post179opt);
    prog.compile( );
}


int main(int argc, char ** argv)
{
    ql::utils::logger::set_log_level("LOG_DEBUG");
    ql::options::set("scheduler_uniform", "no");
    ql::options::set("scheduler_commute", "yes");

//  test_wait_simple("wait_simple", "ASAP", "no");
//  test_wait_simple("wait_simple", "ASAP", "yes");
//  test_wait_simple("wait_simple", "ALAP", "no");
//  test_wait_simple("wait_simple", "ALAP", "yes");

//  test_wait_parallel("wait_parallel", "ASAP", "no");
//  test_wait_parallel("wait_parallel", "ASAP", "yes");
//  test_wait_parallel("wait_parallel", "ALAP", "no");
//  test_wait_parallel("wait_parallel", "ALAP", "yes");

//  test_wait_barrier("wait_barrier", "ASAP", "no");
//  test_wait_barrier("wait_barrier", "ASAP", "yes");
//  test_wait_barrier("wait_barrier", "ALAP", "no");
//  test_wait_barrier("wait_barrier", "ALAP", "yes");

    test_wait_wouter("wait_wouter", "ASAP", "no");
    test_wait_wouter("wait_wouter", "ASAP", "yes");
    test_wait_wouter("wait_wouter", "ALAP", "no");
    test_wait_wouter("wait_wouter", "ALAP", "yes");

    return 0;
}
