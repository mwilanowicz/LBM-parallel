# Project Schedule

## Development Timeline

* **Stage 1 - Report I (26.03): Baseline Sequential Solver**
    * Design of the core D2Q9 solver architecture in C++20.
    * Execution of a physically correct single-threaded simulation.
    * Algorithm verification via mass conservation monitoring.
* **Stage 2 - Report II (16.04): Physical Validation & Boundary Conditions**
    * Implementation of "bounce-back" boundary conditions for walls.
    * Simulation of the Lid-Driven Cavity (LDC) test case.
    * Comparison of velocity profiles with reference literature data.
* **Stage 3 - Report III (07.05): Parallelism & Data Determinism**
    * Domain decomposition for two execution units (threads).
    * Implementation of synchronisation mechanisms to ensure data consistency.
    * Verification of determinism (sequential vs. parallel output identity).
* **Stage 4 - Report IV (21.05): Prototype Performance Analysis**
    * Execution time benchmarks for various grid sizes.
    * Evaluation of speedup achieved during the transition to multi-threading.
    * Formulation of conclusions regarding performance and future development.
* **Stage 5 - Finalization (28.05): Final Documentation**
    * Preparation of the final report summarizing research and programming results.
* **Stage 6 - Completion (11.06): Safety Buffer**
    * Time reserved for potential error corrections identified during final testing.

---

## Risk Analysis & Mitigation Strategies

#### 1. Numerical and Physical Risks
* **Computational Instability (Stage 2):**
    * **Risk:** Potential loss of simulation stability (numerical divergence) when attempting to resolve high-dynamics flows.
    * **Mitigation:** Continuous monitoring of input parameters and enforcing upper limits on maximum velocities within the system.
* **Mass Conservation Errors (Stage 2):**
    * **Risk:** Faulty boundary condition implementation resulting in unnatural mass increase or decrease.
    * **Mitigation:** Systematic verification of the total system density at every computational step.

#### 2. Implementation and Technical Risks
* **Parallel Processing Issues (Stage 3):**
    * **Risk:** Data synchronisation errors during domain decomposition leading to results diverging from the baseline version.
    * **Mitigation:** Deployment of synchronisation barriers and conducting rigorous comparison tests with the sequential version.
* **Low Code Efficiency (Stage 4):**
    * **Risk:** Unsatisfactory speedup after switching to two threads due to memory access latencies.
    * **Mitigation:** Utilisation of optimised data structures (Structure of Arrays - SoA) and hardware resource load analysis.

#### 3. Schedule Risks
* **Reference Data Divergence (Stages 2 & 4):**
    * **Risk:** Identification of significant differences between obtained results and literature data in the final project phase.
    * **Mitigation:** Inclusion of a dedicated time buffer (Stage 6) for model recalibration and additional physical accuracy verification.