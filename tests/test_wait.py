import os
import unittest
from openql import openql as ql
from utils import file_compare

rootDir = os.path.dirname(os.path.realpath(__file__))

curdir = os.path.dirname(__file__)
output_dir = os.path.join(curdir, 'test_output')


class Test_wait(unittest.TestCase):

    def setUp(self):
        ql.set_option('output_dir', output_dir)
        ql.set_option('optimize', 'no')
        ql.set_option('scheduler', 'ASAP')
        ql.set_option('log_level', 'LOG_WARNING')
        ql.set_option('scheduler_post179', 'yes')
        ql.set_option("scheduler_commute", 'yes')
        

    def test_wait_simple(self):

        config_fn = os.path.join(curdir, 'hardware_config_cc_light.json')
        platform = ql.Platform('seven_qubits_chip', config_fn)
        sweep_points = [1, 2]
        num_qubits = platform.get_qubit_number()
        p = ql.Program('test_wait_simple', platform, num_qubits)
        p.set_sweep_points(sweep_points, len(sweep_points))

        k = ql.Kernel('aKernel', platform, num_qubits)

        k.gate("x", [0])
        k.gate("wait", [0], 40)  # OR k.wait([0], 40)
        k.gate("x", [0])

        p.add_kernel(k)
        p.compile()

        QISA_fn = os.path.join(output_dir, p.name+'.qisa')
        gold_fn = rootDir + '/golden/test_wait_simple.qisa'
        self.assertTrue(file_compare(QISA_fn, gold_fn))

    def test_wait_parallel(self):

        config_fn = os.path.join(curdir, 'hardware_config_cc_light.json')
        platform = ql.Platform('seven_qubits_chip', config_fn)
        sweep_points = [1, 2]
        num_qubits = platform.get_qubit_number()
        p = ql.Program('test_wait_parallel', platform, num_qubits)
        p.set_sweep_points(sweep_points, len(sweep_points))

        k = ql.Kernel('aKernel', platform, num_qubits)

        # wait should not be in parallel with another gate
        k.gate("x", [0])
        k.gate("wait", [1], 20)  # OR k.wait([0], 20)
        k.gate("x", [1])

        p.add_kernel(k)
        p.compile()

        QISA_fn = os.path.join(output_dir, p.name+'.qisa')
        gold_fn = rootDir + '/golden/test_wait_parallel.qisa'
        self.assertTrue(file_compare(QISA_fn, gold_fn))

    def test_wait_sweep(self):

        config_fn = os.path.join(curdir, 'hardware_config_cc_light.json')
        platform = ql.Platform('seven_qubits_chip', config_fn)
        sweep_points = [1, 2]
        num_qubits = 7
        p = ql.Program('test_wait_sweep', platform, num_qubits)
        p.set_sweep_points(sweep_points, len(sweep_points))

        qubit_idx = 0
        waits = [20, 40, 60, 100, 200, 400, 800, 1000, 2000]
        for kno, wait_nanoseconds in enumerate(waits):
            k = ql.Kernel("kernel_"+str(kno), platform, num_qubits)

            k.prepz(qubit_idx)

            k.gate('rx90', [qubit_idx])
            k.gate("wait", [qubit_idx], wait_nanoseconds)

            k.gate('rx180', [qubit_idx])
            k.gate("wait", [qubit_idx], wait_nanoseconds)

            k.gate('rx90', [qubit_idx])
            k.gate("wait", [qubit_idx], wait_nanoseconds)

            k.gate('measure', [qubit_idx])

            # add the kernel to the program
            p.add_kernel(k)

        # compile the program
        p.compile()

        QISA_fn = os.path.join(output_dir, p.name+'.qisa')
        gold_fn = rootDir + '/golden/test_wait_sweep.qisa'
        self.assertTrue(file_compare(QISA_fn, gold_fn))

    def test_wait_multi(self):

        config_fn = os.path.join(curdir, 'hardware_config_cc_light.json')
        platform = ql.Platform('seven_qubits_chip', config_fn)
        sweep_points = [1, 2]
        num_qubits = platform.get_qubit_number()
        p = ql.Program('test_wait_multi', platform, num_qubits)
        p.set_sweep_points(sweep_points, len(sweep_points))

        k = ql.Kernel('aKernel', platform, num_qubits)

        for i in range(4):
            k.gate("x", [i])

        k.gate("wait", [0, 1, 2, 3], 40)
        k.wait([0, 1, 2, 3], 40)

        for i in range(4):
            k.gate("measure", [i])

        p.add_kernel(k)
        p.compile()

        QISA_fn = os.path.join(output_dir, p.name+'.qisa')
        gold_fn = rootDir + '/golden/test_wait_multi.qisa'
        self.assertTrue(file_compare(QISA_fn, gold_fn))

    def test_wait_barrier(self):

        config_fn = os.path.join(curdir, 'hardware_config_cc_light.json')
        platform = ql.Platform('seven_qubits_chip', config_fn)
        sweep_points = [1, 2]
        num_qubits = platform.get_qubit_number()
        p = ql.Program('test_wait_barrier', platform, num_qubits)
        p.set_sweep_points(sweep_points, len(sweep_points))

        k = ql.Kernel('aKernel', platform, num_qubits)

        k.gate("x", [0])
        k.gate("x", [1])
        k.gate("y", [0])
        k.gate("wait", [0, 1], 0)  # this will serve as barrier
        k.gate("measure", [0])
        k.gate("measure", [1])

        p.add_kernel(k)
        p.compile()

        QISA_fn = os.path.join(output_dir, p.name+'.qisa')
        gold_fn = rootDir + '/golden/test_wait_barrier.qisa'
        self.assertTrue(file_compare(QISA_fn, gold_fn))

    def test_barrier(self):

        config_fn = os.path.join(curdir, 'hardware_config_cc_light.json')
        platform = ql.Platform('seven_qubits_chip', config_fn)
        sweep_points = [1, 2]
        num_qubits = platform.get_qubit_number()
        p = ql.Program('test_barrier', platform, num_qubits)
        p.set_sweep_points(sweep_points, len(sweep_points))

        k = ql.Kernel('aKernel', platform, num_qubits)

        k.gate("x", [0])
        k.gate("x", [1])
        k.gate("y", [0])

        # k.barrier([0, 1])
        # OR
        k.gate("barrier", [0, 1])

        k.gate("measure", [0])
        k.gate("measure", [1])

        p.add_kernel(k)
        p.compile()

        QISA_fn = os.path.join(output_dir, p.name+'.qisa')
        gold_fn = rootDir + '/golden/test_barrier.qisa'
        self.assertTrue(file_compare(QISA_fn, gold_fn))

    def test_wait_barrier_222(self):

        ql.set_option('scheduler', 'ALAP')
        config_fn = os.path.join(curdir, 'hardware_config_cc_light.json')
        platform = ql.Platform('seven_qubits_chip', config_fn)
        sweep_points = [1]
        num_qubits = platform.get_qubit_number()
        p = ql.Program('test_wait_barrier_222', platform, num_qubits)
        p.set_sweep_points(sweep_points, len(sweep_points))

        k = ql.Kernel('aKernel', platform, num_qubits)

        x, xN, xE, xW, xS, z = 0, 1, 2, 3, 4, 5 

        k.gate('ry90', [x])
        k.gate('ry90', [xN])
        k.gate('ry90', [xE])
        k.gate('ry90', [xW])
        k.gate('ry90', [xS])
        k.gate('barrier', [x, xN, xE, xW, xS])
        # which is same as:
        # k.gate('wait', [x, xN, xE, xW, xS], 0)

        k.gate('measure', [x])
        k.gate('barrier', [x])
        # which is same as:
        # k.gate('wait', [x], 0)

        k.gate('rx90', [z])

        p.add_kernel(k)
        p.compile()

        # QISA_fn = os.path.join(output_dir, p.name+'.qisa')
        # gold_fn = rootDir + '/golden/test_barrier.qisa'
        # self.assertTrue(file_compare(QISA_fn, gold_fn))


if __name__ == '__main__':
    unittest.main()
