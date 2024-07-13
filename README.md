This version of solver solves coexisting phases for a given parent concentration vector.

Please directly execute the file example.py to find phase-coexistence results of an examplary case (from Li and Jacobs, 2024).
Initial guess is included in the file. Note that a good initial guess is pivotal to the success of finding the coexisting phases.  

Users define the free energy function in free_energy.py. 

coexist.py constitutes the core part of the solver, which in principle can be applied to any physical free energy functions.
