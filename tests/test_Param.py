import os
import unittest
from openql import openql as ql

curdir = os.path.dirname(os.path.realpath(__file__))
config_fn = os.path.join(curdir, 'hardware_config_qx.json')
platf = ql.Platform("platform", config_fn)

output_dir = os.path.join(curdir, 'test_output')
nqubits = 10
k = ql.Kernel("k", platf, nqubits)
k2 = ql.Kernel("k", platf, nqubits)
p = ql.Program('test_Param', platf, nqubits)
p2 = ql.Program('test_Param_golden', platf, nqubits)
c = ql.Compiler("testCompiler")
c.add_pass_alias("Writer", "outputIR")
c.set_pass_option("outputIR", "write_qasm_files", "yes")

TypeStrList = ["INT", "REAL", "COMPLEX", "ANGLE"]

def file_compare_no_wait(self, fn1, fn2):
    with open(fn1, 'r') as f1:
        with open(fn2, 'r') as f2:
            for linef1 in f1:
                if(linef1.startswith("    wait")):
                        break
                linef2 = f2.readline()
                if(linef2.startswith("    wait")):
                        linef2 = f2.readline()
                self.assertEqual(linef1, linef2)
    return False

class Test_Param(unittest.TestCase):

    @classmethod
    def setUp(self):
        ql.initialize()
        ql.set_option('output_dir', output_dir)
        ql.set_option('optimize', 'no')
        ql.set_option('scheduler', 'ALAP')
        ql.set_option('log_level', 'LOG_WARNING')
        ql.set_option('log_level', 'LOG_DEBUG')

    # no parameters with Param creation
    def test_param_notype(self):
        param = ql.Param()
        
    # All allowed type names for Param creation
    def test_param_typeStr(self):
        for typeStr in TypeStrList:
            with self.subTest(typeStr=typeStr):
                param = ql.Param(typeStr)
                self.assertEqual(param.typeStr, typeStr)
        
    # Checks it raises an exception about the unknown typename
    # split at ',' so the text specifiying allowed types can be extended when new types are added to the definition
    def test_param_undefinedtype(self):
        typeStr = "UNDEFINEDTYPE"
        with self.assertRaises(Exception) as cm:
            ql.Param(typeStr)
        expectedException = 'Error : Type ' + typeStr + ' is not known'
        self.assertEqual(str(cm.exception).split(',', maxsplit=1)[0], expectedException)

    
    # All allowed type names for Param creation
    # Verify that name is passed on
    def test_param_name(self):
        for typeStr in TypeStrList:
            with self.subTest(typeStr=typeStr):
                name = "p1"
                param = ql.Param(typeStr, name)
                self.assertEqual(param.name, name)
    
    def test_param_randomname(self):
        for typeStr in TypeStrList:
            with self.subTest(typeStr=typeStr):
                param = ql.Param(typeStr)
                self.assertNotEqual(param.name, "")

    # All allowed type names for Param creation
    # Verify that value is set
    def test_param_value_int(self):
        for typeStr in TypeStrList:
            with self.subTest(typeStr=typeStr):
                name = "p1"
                value = 4
                param = ql.Param(typeStr, name, value)
                self.assertEqual(param.name, name)
                self.assertEqual(param.value, value)

    def test_param_real_val(self):
        typeStr = "REAL" 
        name = "p1"
        value = 2.5
        param = ql.Param(typeStr, name, value)
        self.assertEqual(param.name, name)
        self.assertEqual(param.value, value)

    def test_param_complex_val(self):
        typeStr = "COMPLEX"
        name = "p1"
        value = 2 + 3j
        param = ql.Param(typeStr, name, value)
        self.assertEqual(param.name, name)
        self.assertEqual(param.value, value)
    
    def test_param_set_value(self):
        for typeStr in TypeStrList:
            with self.subTest(typeStr=typeStr):
                value = 4
                param = ql.Param(typeStr)
                param.set_value(value)
                self.assertEqual(param.value, value)

    def test_param_gate_int_oldcompile(self):
        k = ql.Kernel("k", platf, nqubits)
        p = ql.Program('test_Param', platf, nqubits)
        typeStr = "INT"
        value = 4
        param = ql.Param(typeStr)
        k.gate("hadamard", param)

        p.add_kernel(k)

        param.set_value(value)
        p.compile()

    def test_param_gate_int_compiler(self):
        k = ql.Kernel("k", platf, nqubits)
        k2 = ql.Kernel("k", platf, nqubits)
        p = ql.Program('test_Param', platf, nqubits)
        p2 = ql.Program('test_Param_golden', platf, nqubits)
        typeStr = "INT"
        value = 4
        param = ql.Param(typeStr)
        k.gate("h", param)
        k2.gate("h", [value])
        
        p.add_kernel(k)
        p2.add_kernel(k2)

        param.set_value(value)
        c.compile(p)
        c.compile(p2)

        PARAM_fn = os.path.join(output_dir, p.name + '.qasm')
        GOLD_fn = os.path.join(output_dir, p2.name + '.qasm')
        file_compare_no_wait(self, PARAM_fn, GOLD_fn)

    def test_param_gate_real_compiler(self):
        k = ql.Kernel("k", platf, nqubits)
        k2 = ql.Kernel("k", platf, nqubits)
        p = ql.Program('test_Param', platf, nqubits)
        p2 = ql.Program('test_Param_golden', platf, nqubits)
        typeStr = "REAL"
        value = 3.1415
        param = ql.Param(typeStr)
        k.gate("rz", [0], 0, param)
        k2.gate("rz", [0], 0, value)
        
        p.add_kernel(k)
        p2.add_kernel(k2)

        param.set_value(value)
        c.compile(p)
        c.compile(p2)

        PARAM_fn = os.path.join(output_dir, p.name + '.qasm')
        GOLD_fn = os.path.join(output_dir, p2.name + '.qasm')
        file_compare_no_wait(self, PARAM_fn, GOLD_fn)

    def test_param_gate_int_valuecompiletime(self):
        k = ql.Kernel("k", platf, nqubits)
        k2 = ql.Kernel("k", platf, nqubits)
        p = ql.Program('test_Param', platf, nqubits)
        p2 = ql.Program('test_Param_golden', platf, nqubits)
        typeStr = "INT"
        value = 3
        param = ql.Param(typeStr)
        k.gate("rz", param, 0, value*0.5)
        k2.gate("rz", [value], 0, value*0.5)
        
        p.add_kernel(k)
        p2.add_kernel(k2)

        c.compile(p, [param], [value])
        c.compile(p2) 

        PARAM_fn = os.path.join(output_dir, p.name + '.qasm')
        GOLD_fn = os.path.join(output_dir, p2.name + '.qasm')
        file_compare_no_wait(self, PARAM_fn, GOLD_fn)

    def test_param_gate_real_valuecompiletime(self):
        k = ql.Kernel("k", platf, nqubits)
        k2 = ql.Kernel("k", platf, nqubits)
        p = ql.Program('test_Param', platf, nqubits)
        p2 = ql.Program('test_Param_golden', platf, nqubits)
        typeStr = "REAL"
        value = 3.1415
        param = ql.Param(typeStr)
        k.gate("rz", [0], 0, param)
        k2.gate("rz", [0], 0, value)
        
        p.add_kernel(k)
        p2.add_kernel(k2)

        c.compile(p, [param], [value])
        c.compile(p2)

        PARAM_fn = os.path.join(output_dir, p.name + '.qasm')
        GOLD_fn = os.path.join(output_dir, p2.name + '.qasm')
        file_compare_no_wait(self, PARAM_fn, GOLD_fn)

    def test_param_gate_int_real_valuecompiletime(self):
        k = ql.Kernel("k", platf, nqubits)
        k2 = ql.Kernel("k", platf, nqubits)
        p = ql.Program('test_Param', platf, nqubits)
        p2 = ql.Program('test_Param_golden', platf, nqubits)
        typeStr = "INT"
        value = 2
        param = ql.Param(typeStr)

        typeStr = "REAL"
        value2 = 3.1415
        param2 = ql.Param(typeStr)
        pint = ql.Param("INT")
        preal = ql.Param("REAL")
        k.gate("rz", param, 0, param2)
        k2.gate("rz", [value], 0, value2)
        
        p.add_kernel(k)
        p2.add_kernel(k2)

        c.compile(p, [param, param2], [value, value2])
        c.compile(p2)

        PARAM_fn = os.path.join(output_dir, p.name + '.qasm')
        GOLD_fn = os.path.join(output_dir, p2.name + '.qasm')
        file_compare_no_wait(self, PARAM_fn, GOLD_fn)
if __name__ == '__main__':
    unittest.main()
