import os
import filecmp
import unittest
from openql import openql as ql

rootDir = os.path.dirname(os.path.realpath(__file__))

curdir = os.path.dirname(__file__)
config_fn = os.path.join(curdir, 'test_cfg_cbox.json')
platf = ql.Platform("starmon", config_fn)

ql.set_output_dir("output")

class Test_qubits(unittest.TestCase):

    def test_1_qubit(self):
        ql.init()
        k = ql.Kernel("aKernel", platf)

        # populate kernel
        k.prepz(0)
        k.x(0)
        k.y(0)
        k.rx90(0)
        k.measure(0)

        nqubits = 1
        sweep_points = [2]
        num_circuits = 1
        p = ql.Program("aProgram", nqubits, platf)
        p.set_sweep_points(sweep_points, num_circuits)

        p.add_kernel(k)  # add kernel to program
        # compile  opt  verbose
        p.compile(False, False)

        gold = rootDir + '/golden/test_1_qubit.qasm'
        isSame = filecmp.cmp('output/aProgram.qasm', gold)
        self.assertTrue(isSame)

    def test_2_qubit(self):
        ql.init()

        k = ql.Kernel("aKernel", platf)

        # populate kernel
        k.prepz(0)
        k.prepz(1)
        k.prepz(2)
        k.cnot(0, 1)
        k.clifford(1, 2)
        k.measure(2)

        nqubits = 3
        sweep_points = [2]
        num_circuits = 1
        p = ql.Program("aProgram", nqubits, platf)
        p.set_sweep_points(sweep_points, num_circuits)

        p.add_kernel(k)  # add kernel to program
        # compile  opt  verbose
        p.compile(False, False)

        gold = rootDir + '/golden/test_2_qubit.qasm'
        isSame = filecmp.cmp('output/aProgram.qasm', gold)
        self.assertTrue(isSame)

    def test_3_qubit(self):
        ql.init()

        k = ql.Kernel("aKernel", platf)

        # populate kernel
        k.prepz(0)
        k.prepz(1)
        k.prepz(2)
        k.toffoli(0, 1, 2)
        k.measure(2)

        nqubits = 3
        sweep_points = [2]
        num_circuits = 1
        p = ql.Program("aProgram", nqubits, platf)
        p.set_sweep_points(sweep_points, num_circuits)

        p.add_kernel(k)  # add kernel to program
        # compile  opt  verbose
        p.compile(False, False)

        gold = rootDir + '/golden/test_3_qubit.qasm'
        isSame = filecmp.cmp('output/aProgram.qasm', gold)
        self.assertTrue(isSame)

if __name__ == '__main__':
    unittest.main()
