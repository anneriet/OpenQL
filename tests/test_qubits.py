import os
import filecmp
import unittest
from openql import openql as ql

rootDir = os.path.dirname(os.path.realpath(__file__))

curdir = os.path.dirname(__file__)
config_fn = os.path.join(curdir, 'test_config_default.json')
platf = ql.Platform("starmon", config_fn)

output_dir = os.path.join(curdir, 'test_output')
ql.set_output_dir(output_dir)


class Test_qubits(unittest.TestCase):

    def test_1_qubit(self):
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
        p = ql.Program("1_qubit_program", nqubits, platf)
        p.set_sweep_points(sweep_points, num_circuits)

        p.add_kernel(k)  # add kernel to program
        # compile  opt  verbose
        p.compile(False, False)

        gold = rootDir + '/golden/test_1_qubit.qasm'
        qasm_fn = os.path.join(output_dir, p.name+'.qasm')
        isSame = filecmp.cmp(qasm_fn, gold)
        self.assertTrue(isSame)

    def test_2_qubit(self):

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
        p = ql.Program("2_qubit_program", nqubits, platf)
        p.set_sweep_points(sweep_points, num_circuits)

        p.add_kernel(k)  # add kernel to program
        # compile  opt  verbose
        p.compile(False, False)

        gold = rootDir + '/golden/test_2_qubit.qasm'
        qasm_fn = os.path.join(output_dir, p.name+'.qasm')
        isSame = filecmp.cmp(qasm_fn, gold)
        self.assertTrue(isSame)

    def test_3_qubit(self):
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
        p = ql.Program("3_qubit_program", nqubits, platf)
        p.set_sweep_points(sweep_points, num_circuits)

        p.add_kernel(k)  # add kernel to program
        # compile  opt  verbose
        p.compile(False, False)

        gold = rootDir + '/golden/test_3_qubit.qasm'
        qasm_fn = os.path.join(output_dir, p.name+'.qasm')
        isSame = filecmp.cmp(qasm_fn, gold)
        self.assertTrue(isSame)

if __name__ == '__main__':
    unittest.main()