version 1.0
# this file has been automatically generated by the OpenQL compiler please do not modify it manually.
qubits 7

.kernel_toffoli_composgate
    y90 q[1]
    xm45 q[1]
    { ym90 q[3] | ym90 q[1] }
    { x q[5] | cz q[1],q[3] }
    { ym90 q[2] | ym90 q[5] }
    { cz q[5],q[2] | y90 q[3] | ym90 q[1] }
    cz q[3],q[1]
    { y90 q[2] | ym90 q[5] }
    { cz q[2],q[5] | y90 q[1] | ym90 q[3] }
    cz q[1],q[3]
    { y90 q[5] | ym90 q[2] }
    { cz q[5],q[2] | y90 q[3] }
    ym90 q[3]
    { y90 q[0] | cz q[5],q[3] }
    xm45 q[0]
    { y q[0] | y90 q[2] | y90 q[3] | ym90 q[5] }
    { cz q[2],q[0] | cz q[3],q[5] }
    wait 1
    { y q[0] | y90 q[5] | ym90 q[3] }
    { x45 q[0] | cz q[5],q[3] }
    y q[0]
    { ym90 q[0] | y90 q[3] }
    cz q[3],q[0]
    wait 1
    { y90 q[0] | ym90 q[3] }
    cz q[0],q[3]
    wait 1
    ym90 q[2]
    { y90 q[3] | ym90 q[0] }
    { cz q[5],q[2] | cz q[3],q[0] }
    wait 1
    y q[2]
    { x45 q[2] | cz q[5],q[3] }
    y q[2]
    cz q[5],q[2]
    y90 q[0]
    { ym90 q[0] | y90 q[2] }
    cz q[2],q[0]
    wait 1
    { y90 q[0] | ym90 q[2] }
    cz q[0],q[2]
    wait 1
    { y90 q[2] | ym90 q[0] }
    cz q[2],q[0]
    wait 1
    y q[3]
    xm45 q[3]
    { y q[3] | y90 q[0] }
    cz q[0],q[3]
    wait 1
    y q[3]
    x45 q[3]
    y q[3]
    cz q[5],q[3]
    wait 1
    y90 q[3]
    measure q[3], b[0]
    wait 10
    y90 q[0]
    xm45 q[0]
    ym90 q[0]
    { ym90 q[2] | x q[0] }
    { ym90 q[0] | cond(b[0]) y90 q[3] | cond(b[0]) y90 q[5] }
    { cond(b[0]) x q[0] | cond(b[0]) xm45 q[3] | cond(b[0]) xm45 q[5] }
    { cond(b[0]) ym90 q[0] | cond(b[0]) y q[3] | cond(b[0]) ym90 q[5] }
    { cond(b[0]) cz q[0],q[3] | cz q[5],q[2] }
    wait 1
    { y90 q[2] | ym90 q[5] }
    { cond(b[0]) y q[3] | cz q[2],q[5] }
    cond(b[0]) x45 q[3]
    { y90 q[5] | ym90 q[2] }
    { cond(b[0]) y q[3] | cz q[5],q[2] }
    ym90 q[3]
    cz q[5],q[3]
    wait 1
    { cond(b[0]) ym90 q[0] | y90 q[2] | y90 q[3] | ym90 q[5] }
    { cond(b[0]) cz q[2],q[0] | cz q[3],q[5] }
    wait 1
    { y90 q[5] | ym90 q[3] }
    cz q[5],q[3]
    cond(b[0]) y q[0]
    { cond(b[0]) x45 q[0] | cond(b[0]) cz q[2],q[5] }
    cond(b[0]) y q[0]
    cond(b[0]) cz q[2],q[0]
    wait 1
    y90 q[3]
    { ym90 q[3] | cond(b[0]) y90 q[0] }
    cz q[0],q[3]
    wait 1
    { y90 q[3] | ym90 q[0] }
    cz q[3],q[0]
    wait 1
    { y90 q[0] | ym90 q[3] }
    { cond(b[0]) y q[5] | cz q[0],q[3] }
    cond(b[0]) xm45 q[5]
    { cond(b[0]) y q[5] | y90 q[3] }
    cond(b[0]) cz q[3],q[5]
    wait 1
    cond(b[0]) y q[5]
    { cond(b[0]) y90 q[3] | cond(b[0]) x45 q[5] }
    { cond(b[0]) xm45 q[3] | cond(b[0]) y q[5] }
    { cond(b[0]) ym90 q[3] | cond(b[0]) cz q[2],q[5] }
    cond(b[0]) x q[3]
    { cond(b[0]) ym90 q[3] | cond(b[0]) y90 q[5] }
    measure q[5], b[0]
    wait 14
    cond(b[0]) y90 q[2]
    cond(b[0]) xm45 q[2]
    { ym90 q[0] | cond(b[0]) ym90 q[2] }
    cz q[2],q[0]
    cond(b[0]) y90 q[5]
    { cond(b[0]) x q[3] | cond(b[0]) xm45 q[5] }
    { cond(b[0]) ym90 q[3] | cond(b[0]) y q[5] }
    { cond(b[0]) cz q[3],q[5] | y90 q[0] | ym90 q[2] }
    cz q[0],q[2]
    cond(b[0]) y q[5]
    { cond(b[0]) x45 q[5] | y90 q[2] | ym90 q[0] }
    { cond(b[0]) y q[5] | cz q[2],q[0] }
    ym90 q[5]
    cz q[2],q[5]
    wait 1
    { y90 q[5] | ym90 q[2] }
    cz q[5],q[2]
    cond(b[0]) ym90 q[3]
    { y90 q[0] | y90 q[2] | ym90 q[5] }
    { cond(b[0]) cz q[0],q[3] | cz q[2],q[5] }
    wait 1
    cond(b[0]) cz q[0],q[2]
    wait 1
    cond(b[0]) y q[3]
    cond(b[0]) x45 q[3]
    cond(b[0]) y q[3]
    cond(b[0]) cz q[0],q[3]
    y90 q[5]
    { ym90 q[5] | cond(b[0]) y90 q[3] }
    cz q[3],q[5]
    wait 1
    { y90 q[5] | ym90 q[3] }
    cz q[5],q[3]
    wait 1
    { y90 q[3] | ym90 q[5] }
    { cond(b[0]) y q[2] | cz q[3],q[5] }
    cond(b[0]) xm45 q[2]
    { cond(b[0]) y q[2] | y90 q[5] }
    cond(b[0]) cz q[5],q[2]
    wait 1
    cond(b[0]) y q[2]
    { cond(b[0]) y90 q[5] | cond(b[0]) x45 q[2] }
    { cond(b[0]) xm45 q[5] | cond(b[0]) y q[2] }
    { cond(b[0]) ym90 q[5] | cond(b[0]) cz q[0],q[2] }
    cond(b[0]) x q[5]
    { cond(b[0]) ym90 q[5] | cond(b[0]) y90 q[2] }
    measure q[2], b[0]
    wait 14
    cond(b[0]) y90 q[0]
    cond(b[0]) xm45 q[0]
    { ym90 q[3] | cond(b[0]) ym90 q[0] }
    cz q[0],q[3]
    wait 1
    { y90 q[3] | ym90 q[0] }
    cz q[3],q[0]
    wait 1
    cond(b[0]) y90 q[2]
    cond(b[0]) xm45 q[2]
    { cond(b[0]) x q[5] | cond(b[0]) y q[2] }
    { cond(b[0]) ym90 q[5] | y90 q[0] | ym90 q[3] }
    { cond(b[0]) cz q[5],q[2] | cz q[0],q[3] }
    wait 1
    cond(b[0]) y q[2]
    cond(b[0]) x45 q[2]
    cond(b[0]) y q[2]
    ym90 q[2]
    cz q[0],q[2]
    wait 1
    { y90 q[2] | ym90 q[0] }
    cz q[2],q[0]
    wait 1
    { cond(b[0]) ym90 q[5] | y90 q[3] }
    { cond(b[0]) cz q[3],q[5] | y90 q[0] | ym90 q[2] }
    cz q[0],q[2]
    cond(b[0]) y q[5]
    { cond(b[0]) x45 q[5] | cond(b[0]) cz q[3],q[0] }
    cond(b[0]) y q[5]
    cond(b[0]) cz q[3],q[5]
    y90 q[2]
    { ym90 q[2] | cond(b[0]) y90 q[5] }
    cz q[5],q[2]
    wait 1
    { y90 q[2] | ym90 q[5] }
    cz q[2],q[5]
    wait 1
    { y90 q[5] | ym90 q[2] }
    { cond(b[0]) y q[0] | cz q[5],q[2] }
    cond(b[0]) xm45 q[0]
    { cond(b[0]) y q[0] | y90 q[2] }
    cond(b[0]) cz q[2],q[0]
    wait 1
    cond(b[0]) y q[0]
    cond(b[0]) x45 q[0]
    cond(b[0]) y q[0]
    cond(b[0]) cz q[3],q[0]
    wait 1
    cond(b[0]) y90 q[2]
    cond(b[0]) xm45 q[2]
    cond(b[0]) ym90 q[2]
    cond(b[0]) x q[2]
    { cond(b[0]) ym90 q[2] | cond(b[0]) y90 q[0] }
