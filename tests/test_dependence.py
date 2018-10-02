import os
import filecmp
import unittest
from openql import openql as ql

rootDir = os.path.dirname(os.path.realpath(__file__))

curdir = os.path.dirname(__file__)
config_fn = os.path.join(curdir, 'test_config_default.json')
platf = ql.Platform("starmon", config_fn)

output_dir = os.path.join(curdir, 'test_output')

def set_ql_options():
    ql.set_option('output_dir', output_dir)
    ql.set_option('optimize', 'no')
    ql.set_option('scheduler', 'ASAP')
    ql.set_option('log_level', 'LOG_INFO')

def file_compare(fn1, fn2):
    isSame = False
    with open(fn1, 'r') as f1:
        with open(fn2, 'r') as f2:
            a = f1.read()
            b = f2.read()
            f1.close()
            f2.close()
            isSame = (a==b)
    return isSame

class Test_dependence(unittest.TestCase):

    # @unittest.expectedFailure
    # @unittest.skip
    def test_independent(self):
        set_ql_options()
        nqubits = 4
        # populate kernel
        k = ql.Kernel("aKernel", platf, nqubits)

        for i in range(4):
            k.prepz(i)

        # no dependence
        k.cz(0, 1)
        k.cz(2, 3)

        k.measure(0)
        k.measure(1)

        sweep_points = [2]

        p = ql.Program("independent", platf, nqubits)
        p.set_sweep_points(sweep_points, len(sweep_points))
        p.add_kernel(k)
        p.compile()

        gold_fn = rootDir + '/golden/test_independence.qasm'
        qasm_fn = os.path.join(output_dir, p.name+'_scheduled.qasm')
        self.assertTrue( file_compare(qasm_fn, gold_fn) )

        # load qasm
        qasm_files = []
        qasm_files.append(os.path.join(output_dir, 'independent.qasm'))
        qasm_files.append(os.path.join(output_dir, 'independent_scheduled.qasm'))

        for qasm_file in qasm_files:
           qasm_reader = ql.QASM_Loader(qasm_file)
           errors = qasm_reader.load()
           self.assertTrue(errors == 0)



    # @unittest.skip
    def test_WAW(self):
        set_ql_options()

        nqubits = 4
        # populate kernel
        k = ql.Kernel("aKernel", platf, nqubits)

        for i in range(4):
            k.prepz(i)

        # q1 dependence
        k.cz(0, 1)
        k.cz(2, 1)

        k.measure(0)
        k.measure(1)

        sweep_points = [2]

        p = ql.Program("WAW", platf, nqubits)
        p.set_sweep_points(sweep_points, len(sweep_points))
        p.add_kernel(k)
        p.compile()

        gold_fn = rootDir + '/golden/test_WAW_ASAP.qasm'
        qasm_fn = os.path.join(output_dir, p.name+'_scheduled.qasm')
        
        print(gold_fn)
        print(qasm_fn)

        self.assertTrue( file_compare(qasm_fn, gold_fn) )

        # load qasm
        qasm_files = []
        qasm_files.append(os.path.join(output_dir, p.name+'.qasm'))
        qasm_files.append(os.path.join(output_dir, p.name+'_scheduled.qasm'))

        for qasm_file in qasm_files:
           qasm_reader = ql.QASM_Loader(qasm_file)
           errors = qasm_reader.load()
           self.assertTrue(errors == 0)


    # @unittest.skip
    def test_RAR_Control(self):
        set_ql_options()
        nqubits = 4

        # populate kernel
        k = ql.Kernel("aKernel", platf, nqubits)

        for i in range(4):
            k.prepz(i)

        # q0 dependence
        k.cz(0, 1)
        k.cz(0, 2)

        k.measure(0)
        k.measure(1)

        sweep_points = [2]

        p = ql.Program("RAR", platf, nqubits)
        p.set_sweep_points(sweep_points, len(sweep_points))
        p.add_kernel(k)
        p.compile()

        gold_fn = rootDir + '/golden/test_RAR_Control_ASAP.qasm'
        qasm_fn = os.path.join(output_dir, p.name+'_scheduled.qasm')

        self.assertTrue( file_compare(qasm_fn, gold_fn) )

        # load qasm
        qasm_files = []
        qasm_files.append(os.path.join(output_dir, 'RAR.qasm'))
        qasm_files.append(os.path.join(output_dir, 'RAR_scheduled.qasm'))

        for qasm_file in qasm_files:
           qasm_reader = ql.QASM_Loader(qasm_file)
           errors = qasm_reader.load()
           self.assertTrue(errors == 0)



    # @unittest.skip
    def test_RAW(self):
        set_ql_options()
        nqubits = 4

        # populate kernel
        k = ql.Kernel("aKernel", platf, nqubits)

        for i in range(4):
            k.prepz(i)

        # q1 dependence
        k.cz(0, 1)
        k.cz(1, 2)

        k.measure(0)
        k.measure(1)

        sweep_points = [2]

        p = ql.Program("RAW", platf, nqubits)
        p.set_sweep_points(sweep_points, len(sweep_points))
        p.add_kernel(k)
        p.compile()

        gold_fn = rootDir + '/golden/test_RAW_ASAP.qasm'
        qasm_fn = os.path.join(output_dir, p.name+'_scheduled.qasm')
        
        self.assertTrue( file_compare(qasm_fn, gold_fn) )

        # load qasm
        qasm_files = []
        qasm_files.append(os.path.join(output_dir, 'RAW.qasm'))
        qasm_files.append(qasm_fn)

        for qasm_file in qasm_files:
           qasm_reader = ql.QASM_Loader(qasm_file)
           errors = qasm_reader.load()
           self.assertTrue(errors == 0)



    # @unittest.skip
    def test_WAR(self):
        set_ql_options()
        nqubits = 4

        # populate kernel
        k = ql.Kernel("aKernel", platf, nqubits)

        for i in range(4):
            k.prepz(i)

        # q0 dependence
        k.cz(0, 1)
        k.cz(2, 0)

        k.measure(0)
        k.measure(1)

        sweep_points = [2]

        p = ql.Program("WAR", platf, nqubits)
        p.set_sweep_points(sweep_points, len(sweep_points))
        p.add_kernel(k)
        p.compile()

        gold_fn = rootDir + '/golden/test_WAR_ASAP.qasm'
        qasm_fn = os.path.join(output_dir, p.name+'_scheduled.qasm')
        
        self.assertTrue( file_compare(qasm_fn, gold_fn) )

        # load qasm
        qasm_files = []
        qasm_files.append(os.path.join(output_dir, 'WAR.qasm'))
        qasm_files.append(qasm_fn)

        for qasm_file in qasm_files:
           qasm_reader = ql.QASM_Loader(qasm_file)
           errors = qasm_reader.load()
           self.assertTrue(errors == 0)


    # @unittest.skip
    def test_swap_single(self):
        set_ql_options()
        nqubits = 4

        # populate kernel
        k = ql.Kernel("aKernel", platf, nqubits)

        k.gate("x", [0]);
        k.gate("swap", [0, 1])
        k.gate("x", [0])

        sweep_points = [2]

        p = ql.Program("swap_single", platf, nqubits)
        p.set_sweep_points(sweep_points, len(sweep_points))
        p.add_kernel(k)
        p.compile()

        gold_fn = rootDir + '/golden/test_swap_single_ASAP.qasm'
        qasm_fn = os.path.join(output_dir, p.name+'_scheduled.qasm')

        self.assertTrue( file_compare(qasm_fn, gold_fn) )

        # load qasm
        qasm_files = []
        qasm_files.append(os.path.join(output_dir, 'swap_single.qasm'))
        qasm_files.append(qasm_fn)

        for qasm_file in qasm_files:
           qasm_reader = ql.QASM_Loader(qasm_file)
           errors = qasm_reader.load()
           self.assertTrue(errors == 0)

    # @unittest.skip
    def test_swap_multi(self):
        set_ql_options()
        nqubits = 5

        # populate kernel
        k = ql.Kernel("aKernel", platf, nqubits)

        # swap test with 2 qubit gates
        k.gate("x", [0])
        k.gate("x", [1])
        k.gate("swap", [0, 1])
        k.gate("cz", [0, 2])
        k.gate("cz", [1, 4])

        sweep_points = [2]

        p = ql.Program("swap_multi", platf, nqubits)
        p.set_sweep_points(sweep_points, len(sweep_points))
        p.add_kernel(k)
        p.compile()

        gold_fn = rootDir + '/golden/test_swap_multi_ASAP.qasm'
        qasm_fn = os.path.join(output_dir, p.name+'_scheduled.qasm')
        self.assertTrue( file_compare(qasm_fn, gold_fn) )

        # load qasm
        qasm_files = []
        qasm_files.append(os.path.join(output_dir, 'swap_multi.qasm'))
        qasm_files.append(qasm_fn)

        for qasm_file in qasm_files:
           qasm_reader = ql.QASM_Loader(qasm_file)
           errors = qasm_reader.load()
           self.assertTrue(errors == 0)


if __name__ == '__main__':
    unittest.main()
