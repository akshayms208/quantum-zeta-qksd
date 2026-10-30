# 🌌 Quantum Krylov Subspace Diagonalization (QKSD) of the Riemann Zeta Spectrum

An ab-initio quantum computational physics model built in Python (Qiskit/SciPy) that computes the non-trivial zeros of the Riemann zeta function by simulating prime-based quantum interactions.

---

## 📖 What is this project?
The **Riemann Hypothesis** is one of the biggest unsolved problems in mathematics. It deals with the "zeros" of the Riemann zeta function, which strictly dictate the distribution of prime numbers. 

The **Hilbert-Pólya conjecture** suggests a fascinating idea: what if these mathematical zeros are actually the energy levels (eigenvalues) of a quantum mechanical system? If we can find that system, we could prove the hypothesis.

This repository explores that conjecture computationally. It takes **Sierra's Russian Doll Hamiltonian** (a theoretical quantum system where particles interact based on prime numbers) and solves for its energy levels using a hybrid quantum-classical algorithm called **Block-QKSD** (Quantum Krylov Subspace Diagonalization). 

**The result?** The algorithm successfully outputs the Riemann zeros just by looking at prime numbers—without ever being explicitly told what the zeros are supposed to be!

---

## 🏗️ Theoretical Architecture 

This project is built on three main pillars:

### 1. The Ab-Initio Prime Sieve 
To make this a true simulation, we cannot hardcode the known Riemann zeros into the system. Instead, the code features a custom "prime sieve" engine. It generates prime numbers from scratch, calculates a conditionally convergent phase correction (the fluctuating part of the Riemann counting function), and uses that raw prime data to construct the energy landscape of our system.

### 2. Sierra's Russian Doll Hamiltonian
How do we turn prime numbers into quantum mechanics? We use a model proposed by physicist Germán Sierra. 
* We map the prime-corrected energy states onto a multi-qubit system.
* The Hamiltonian (the matrix defining the system's total energy) is mapped to Pauli strings (X, Y, Z gates) using Qiskit's `SparsePauliOp`.
* In this model, the quantum states represent reciprocal prime interactions.

### 3. Block-QKSD (Quantum Krylov Subspace Diagonalization)
Finding the eigenvalues of a massive quantum matrix usually requires a method called *Quantum Phase Estimation (QPE)*, which requires impossibly deep circuits that current quantum computers can't handle. 
* Instead, we use **QKSD**. 
* We time-evolve reference states using the quantum circuit ($e^{-iHt}$).
* We measure the overlap of these states (the Krylov subspace).
* Finally, we pass that data to a classical computer (via SciPy) to solve a generalized eigenvalue problem, applying SVD thresholding to stabilize the math.

---

## ⚙️ Installation

To run this simulation on your local machine, you will need Python installed. 

1. **Clone the repository:**
   ```bash
   git clone [https://github.com/akshayms208/quantum-zeta-qksd.git](https://github.com/akshayms208/quantum-zeta-qksd.git)
   cd quantum-zeta-qksd
