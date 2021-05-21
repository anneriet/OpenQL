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
    # Verify name is passed on
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
    # Verify value is set
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
        typeStr = "INT"
        value = 4
        param = ql.Param(typeStr)
        k.gate("hadamard", param)

        p.add_kernel(k)

        param.set_value(value)
        p.compile()

    def test_param_gate_int_compiler(self):
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
        self.assertListEqual(
            list(open(PARAM_fn)),
            list(open(GOLD_fn)))
        self.assertTrue(file_compare(PARAM_fn, GOLD_fn))

if __name__ == '__main__':
    unittest.main()
