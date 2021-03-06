/** \file
 * Implementation for pass that writes circuits using the qsoverlay format.
 */

#include "qsoverlay.h"

#include <iostream>
#include "utils/filesystem.h"
#include "utils/map.h"
#include "options.h"
#include "kernel.h"

namespace ql {

using namespace utils;

//Only support for DiCarlo setup atm
void write_qsoverlay_program(
    quantum_program *programp,
    UInt num_qubits,
    const quantum_platform &platform,
    const Str &suffix,
    UInt ns_per_cycle,
    Bool compiled
) {

    // TODO remove the next line. Using this because qsoverlay has some bugs when time is explicited
    compiled = false;

    QL_IOUT("Writing scheduled QSoverlay program");
    Str qfname(options::get("output_dir") + "/" + "quantumsim_" + programp->unique_name + "_" + suffix + ".py");
    QL_DOUT("Writing scheduled QSoverlay program " << qfname);
    QL_IOUT("Writing scheduled QSoverlay program " << qfname);
    OutFile fout(qfname);

    fout << "# Quantumsim (via Qsoverlay) program generated by OpenQL\n"
         << "# Please modify at your will to obtain extra information from Quantumsim\n\n";

    fout << "import numpy as np\n"
         << "from qsoverlay import DiCarlo_setup\n"
         << "from qsoverlay.circuit_builder import Builder\n";

    fout << "import quantumsim.sparsedm as sparsedm\n"
         << "\n"
         << "# print('GPU is used:', sparsedm.using_gpu)\n"
         << "\n"
         << "\n";

    // Gate correspondence
    Map<Str, Str> gate_map = {
        {"prepz", "prepz"},
        {"x", "X"},
        {"x45", "RX"},
        {"x90", "RX"},
        {"xm45", "RX"},
        {"xm90", "RX"},
        {"y", "Y"},
        {"y45", "RY"},
        {"y90", "RY"},
        {"ym45", "RY"},
        {"ym90", "RY"},
        {"h", "H"},
        {"cz", "CZ"},
        {"measure", "Measure"},
    };

    Map<Str, Str> angles = {
        {"x45", "np.pi/4"},
        {"x90", "np.pi/2"},
        {"xm45", "-np.pi/4"},
        {"xm90", "-np.pi/2"},
        {"y45", "np.pi/4"},
        {"y90", "np.pi/2"},
        {"ym45", "-np.pi/4"},
        {"ym90", "-np.pi/2"},
    };

    if (!compiled) {
        gate_map.set("cnot") = "CNOT";
        // gate_map["t"] = "RZ";
        // angles["t"] = "np.pi/4";
        // gate_map["tdag"] = "RZ";
        // angles["t"] = "-np.pi/4";
    }

    // Create qubit list

    Str qubit_list{};
    for (UInt qubit = 0; qubit < num_qubits; qubit++) {
        qubit_list += "'";
        qubit_list += (qubit != num_qubits-1) ? (to_string(qubit) + "', ") : (to_string(qubit) + "'");
    }

    // Circuit creation
    fout << "\n#Now the circuit is created\n"

         << "\ndef circuit_generated(noise_flag, setup_name = 'DiCarlo_setup'):\n"
         << "	qubit_list = [" + qubit_list + "]\n"
         << "	if setup_name == 'DiCarlo_setup':\n"
         << "		setup = DiCarlo_setup.quick_setup(qubit_list, noise_flag = noise_flag)\n"
         << "	b = Builder(setup)\n"
         << "	b.new_circuit(circuit_title = '" << programp->kernels.front().name << "')\n";


    // Circuit creation: Add gates
    for (auto & gate: programp->kernels.front().c) {
        Str qs_name;
        try {
            qs_name = gate_map.at(gate->name);
        } catch (std::exception &e) {
            // WOUT("Next gate: " + gate->name + " .... WRONG");
            QL_EOUT("Qsoverlay: unknown gate detected!: " + gate->name);
            throw Exception("Qsoverlay: unknown gate detected!:" + gate->name, false);
        }

        // IOUT(gate->name);
        if (gate->operands.size() == 1) {
            QL_IOUT("Gate operands: " + to_string(gate->operands[0]));
        } else if (gate->operands.size() == 2) {
            QL_IOUT("Gate operands: " + to_string(gate->operands[0]) + ", " + to_string(gate->operands[1]));
        } else {
            QL_IOUT("GATE OPERANDS: Problem encountered");
        }


        fout << "	b.add_gate('" << qs_name  << "', " << "['"
             << to_string(gate->operands[0])
             << (( gate->operands.size() == 1 ) ? "']" : ("', '" + to_string(gate->operands[1]) + "']"));


        // Add angles for the gates that require it
        if (qs_name == "RX" || qs_name == "RY" || qs_name == "t" || qs_name == "tdag") {
            fout << ", angle = " << angles.dbg(gate->name);
        }

        // Add gate timing, if circuit was compiled.
        if (qs_name == "prepz") {
            if (compiled) {
                fout << ", time = " << to_string((gate->cycle - 1) * ns_per_cycle + gate->duration);
            }
        } else if (qs_name == "Measure") {
            fout << ", output_bit = " << "'" << gate->operands[0] << "_out'";
            if (compiled) {
                fout << ", time = " << to_string((gate->cycle-1)*ns_per_cycle + gate->duration/4);
            }
        } else {
            if (compiled) {
                fout << ", time = " << to_string((gate->cycle-1)*ns_per_cycle + gate->duration/2);
            }
        }
        fout << ")\n";
    }

    fout << "\n"
         << "	b.finalize()\n"
         << "	return b.circuit\n";

    fout.close();
    QL_IOUT("Writing scheduled QSoverlay program [Done]");
}

} // namespace ql
