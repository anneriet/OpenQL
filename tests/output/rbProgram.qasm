# this file has been automatically generated by the OpenQL compiler please do not modify it manually.
qubits 3

.kernel1
   prepz q0
   prepz q1
   x q0
   y q0
   cnot q0, q1
   toffoli q0, q1, q2
   rx90 q0
   mrx90 q0
   mry90 q0
   measure q2
.cal0_1
   prepz q0
   measure q0
.cal0_2
   prepz q0
   measure q0
.cal1_1
   prepz q0
   x q0
   measure q0
.cal1_2
   prepz q0
   x q0
   measure q0
