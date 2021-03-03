version 1.0
# this file has been automatically generated by the OpenQL compiler please do not modify it manually.
qubits 7

.kernel_allD
    { x q[0] | ym90 q[3] }
    cz q[0],q[3]
    wait 1
    cz q[0],q[3]
    wait 1
    ym90 q[2]
    cz q[0],q[2]
    wait 1
    { y90 q[3] | ym90 q[0] }
    cz q[3],q[0]
    wait 1
    { y90 q[0] | ym90 q[3] }
    cz q[0],q[3]
    wait 1
    { y90 q[3] | ym90 q[5] }
    cz q[3],q[5]
    ym90 q[1]
    cz q[3],q[1]
    ym90 q[6]
    cz q[3],q[6]
    y90 q[1]
    { x q[1] | ym90 q[4] | y90 q[5] }
    { cz q[1],q[4] | cz q[5],q[2] }
    wait 1
    { ym90 q[1] | y90 q[2] | ym90 q[5] }
    { cz q[3],q[1] | cz q[2],q[5] }
    wait 1
    { y90 q[1] | ym90 q[3] | y90 q[5] | ym90 q[2] }
    { cz q[1],q[3] | cz q[5],q[2] }
    wait 1
    { y90 q[3] | ym90 q[5] }
    cz q[3],q[5]
    wait 1
    cz q[3],q[6]
    wait 1
    y90 q[5]
    x q[5]
    { ym90 q[0] | cz q[5],q[2] }
    cz q[3],q[0]
    { ym90 q[1] | cz q[5],q[2] }
    cz q[3],q[1]
    { y90 q[2] | ym90 q[5] }
    { y90 q[1] | cz q[2],q[5] }
    cz q[1],q[4]
    { y90 q[5] | ym90 q[2] }
    { ym90 q[1] | cz q[5],q[2] }
    cz q[3],q[1]
    y90 q[2]
    { cz q[2],q[0] | ym90 q[5] }
    cz q[3],q[5]
    y90 q[0]
    { x q[0] | ym90 q[3] }
    cz q[0],q[3]
    wait 1
    ym90 q[2]
    cz q[0],q[2]
    wait 1
    { y90 q[3] | y90 q[2] | ym90 q[0] }
    { cz q[3],q[6] | cz q[2],q[0] }
    wait 1
    { y90 q[6] | ym90 q[3] }
    cz q[6],q[3]
    ym90 q[2]
    { y90 q[3] | ym90 q[6] | y90 q[0] }
    { cz q[3],q[6] | cz q[0],q[2] }
    wait 1
    ym90 q[3]
    cz q[0],q[3]
    wait 1
    { y90 q[6] | cz q[0],q[3] }
    cz q[6],q[4]
    { y90 q[3] | ym90 q[0] }
    { y90 q[4] | ym90 q[6] | cz q[3],q[0] }
    cz q[4],q[6]
    y90 q[2]
    { y90 q[6] | ym90 q[4] | y90 q[0] | ym90 q[3] }
    { cz q[6],q[4] | cz q[0],q[3] }
    cz q[2],q[5]
    { y90 q[3] | ym90 q[6] | ym90 q[0] }
    { cz q[3],q[6] | cz q[2],q[0] }
    wait 1
    { y90 q[6] | cz q[2],q[0] }
    cz q[6],q[4]
    { y90 q[0] | ym90 q[2] }
    { y90 q[4] | ym90 q[6] | cz q[0],q[2] }
    cz q[4],q[6]
    ym90 q[0]
    { y90 q[6] | ym90 q[4] | cz q[3],q[0] }
    cz q[6],q[4]
    cz q[3],q[1]
    { ym90 q[6] | y90 q[2] }
    { cz q[3],q[6] | cz q[2],q[0] }
    wait 1
    { y90 q[6] | y90 q[0] | ym90 q[3] }
    { cz q[6],q[4] | cz q[0],q[3] }
    wait 1
    { y90 q[4] | ym90 q[6] | y90 q[3] }
    { cz q[4],q[6] | cz q[3],q[1] }
    wait 1
    { y90 q[6] | ym90 q[4] }
    { ym90 q[0] | cz q[6],q[4] }
    cz q[3],q[0]
    ym90 q[6]
    cz q[3],q[6]
    wait 1
    cz q[3],q[0]
    y90 q[6]
    { x q[6] | ym90 q[3] }
    cz q[6],q[3]
    wait 1
    cz q[6],q[3]
    wait 1
    cz q[6],q[4]
    wait 1
    { y90 q[3] | ym90 q[6] }
    cz q[3],q[6]
    wait 1
    { y90 q[6] | ym90 q[3] }
    cz q[6],q[3]
    wait 1
    y90 q[3]
    cz q[3],q[5]
    wait 1
    y90 q[5]
    { x q[5] | ym90 q[2] }
    cz q[5],q[2]
    wait 1
    { y90 q[2] | ym90 q[5] }
    cz q[2],q[5]
    wait 1
    { y90 q[5] | ym90 q[2] }
    cz q[5],q[2]
    cz q[3],q[0]
    ym90 q[5]
    cz q[3],q[5]
    y90 q[2]
    cz q[2],q[5]
    cz q[3],q[1]
    y90 q[5]
    { x q[5] | ym90 q[3] }
    { cz q[2],q[0] | cz q[5],q[3] }
    wait 1
    { y90 q[3] | ym90 q[6] }
    { cz q[2],q[0] | cz q[3],q[6] }
    wait 1
    { y90 q[0] | ym90 q[2] | y90 q[6] | ym90 q[3] }
    { cz q[0],q[2] | cz q[6],q[3] }
    wait 1
    { y90 q[2] | ym90 q[0] | y90 q[3] | ym90 q[6] }
    { cz q[2],q[0] | cz q[3],q[6] }
    wait 1
    { y90 q[0] | ym90 q[3] | y90 q[6] }
    { cz q[0],q[3] | cz q[6],q[4] }
    wait 1
    cz q[5],q[3]
    wait 1
    { y90 q[4] | y90 q[3] }
    { ym90 q[6] | x q[3] | ym90 q[0] }
    { cz q[4],q[6] | cz q[3],q[0] }
    wait 1
    { y90 q[6] | ym90 q[4] | y90 q[0] | ym90 q[3] }
    { cz q[6],q[4] | cz q[0],q[3] }
    wait 1
    { y90 q[3] | ym90 q[6] }
    cz q[3],q[6]
    wait 1
    y90 q[6]
    cz q[6],q[4]
    wait 1
    { y90 q[4] | ym90 q[6] }
    { cz q[3],q[1] | cz q[4],q[6] }
    wait 1
    { y90 q[6] | ym90 q[4] }
    { ym90 q[0] | cz q[6],q[4] }
    cz q[3],q[0]
    { ym90 q[6] | y90 q[4] }
    { cz q[3],q[6] | cz q[4],q[1] }
    wait 1
    ym90 q[3]
    { cz q[5],q[3] | y90 q[1] | ym90 q[4] }
    cz q[1],q[4]
    { ym90 q[2] | x q[6] }
    { cz q[5],q[2] | y90 q[4] | y q[6] }
    { y90 q[3] | cz q[4],q[6] }
    { x q[3] | ym90 q[5] | ym90 q[1] }
    { cz q[3],q[5] | cz q[4],q[1] }
    wait 1
    { y90 q[6] | ym90 q[4] }
    cz q[6],q[4]
    { y90 q[5] | ym90 q[3] }
    cz q[5],q[3]
    wait 1
    { y90 q[3] | ym90 q[6] }
    cz q[3],q[6]
    y90 q[4]
    { cz q[3],q[1] | cz q[4],q[6] }
    ym90 q[5]
    cz q[3],q[5]
    wait 1
    { y90 q[0] | y90 q[1] | y90 q[2] | y90 q[5] | y90 q[6] }
    { x q[1] | x q[2] | x q[3] | x q[6] }
