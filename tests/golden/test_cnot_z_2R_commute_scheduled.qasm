version 1.0
# this file has been automatically generated by the OpenQL compiler please do not modify it manually.
qubits 7

.aKernel
    { t q[0] | z q[0] | t q[3] | z q[3] }
    wait 2
    cnot q[3],q[0]
    wait 3
