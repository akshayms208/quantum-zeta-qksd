# Quantum Krylov Subspace Diagonalization (QKSD) of the Riemann Zeta Spectrum

An ab-initio quantum model using QKSD in Qiskit to compute Riemann zeta zeros. It encodes Sierra's Russian Doll Hamiltonian into a multi-qubit system and features a custom prime sieve to simulate prime-based quantum interactions without deep phase estimation.

## Theoretical Architecture
* **Ab-Initio Number Theory Engine:** Custom prime sieve generates sequence parameters from scratch, avoiding hardcoded knowledge of the Riemann spectrum.
* **Russian Doll Hamiltonian:** Encodes reciprocal prime states into a projective multi-qubit basis via Qiskit's `PauliOp` and `SparsePauliOp`.
* **Block-QKSD Algorithm:** Extracts the low-lying eigenvalue spectrum by generating Krylov basis states, applying parameterized unitary evolution, and thresholding linear dependencies to stabilize calculations.

## Installation

Clone the repository and install the required dependencies:

```bash
git clone [https://github.com/akshayms208/quantum-zeta-qksd.git](https://github.com/akshayms208/quantum-zeta-qksd.git)
cd quantum-zeta-qksd
pip install -r requirements.txt