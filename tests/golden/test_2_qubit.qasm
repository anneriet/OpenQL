# this file has been automatically generated by the OpenQL compiler please do not modify it manually.
qubits 3

.aKernel
   prepz q0
   prepz q1
   prepz q2
   cnot q0,q1
   ry90 q2
   rx90 q2
   measure q2
