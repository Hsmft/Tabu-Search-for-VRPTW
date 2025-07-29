# Tabu Search for the Vehicle Routing Problem with Time Windows (VRPTW)

![Language](https://img.shields.io/badge/language-C%2B%2B-blue.svg)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains a C++ implementation of the **Tabu Search** metaheuristic designed to solve the **Vehicle Routing Problem with Time Windows (VRPTW)**. The goal of this project is to find high-quality, near-optimal solutions by minimizing the number of vehicles and the total travel distance, while respecting all problem constraints.

---

## 📋 Project Overview

The Vehicle Routing Problem with Time Windows (VRPTW) is a classic NP-hard optimization problem. This solver finds effective solutions using a Tabu Search algorithm, which is known for its ability to escape local optima by using memory structures.

The solver implements the following workflow:

1.  **Instance Parsing:** Reads complex VRPTW problem instances from standard text files.
2.  **Initial Solution:** Generates a feasible starting solution using a **greedy nearest-neighbor heuristic**.
3.  **Optimization with Tabu Search:** Improves the solution using a Tabu Search core engine that features:
    * **Short-Term Memory:** A tabu list (implemented as a queue) that prevents the search from revisiting recent solutions, helping it to avoid cycles.
    * **Long-Term Memory:** A frequency-based memory structure is used to penalize frequently used moves, encouraging diversification in the search space.
    * **Aspiration Criteria:** Allows a tabu move if it leads to a new best-ever solution.
4.  **Neighborhood Operators:** The search explores the solution space using a variety of neighborhood "moves," including:
    * 2-Opt (intra-route)
    * Relocate Customer (inter-route)
    * Swap Customers (inter-route)
    * CROSS-Exchange
5.  **Feasibility Checks:** Ensures that all generated routes adhere to vehicle capacity and customer time window constraints.

---

## 🛠️ Technologies Used

* **Language:** C++ (utilizing C++11 features like `<chrono>` and `<random>`)
* **Libraries:** C++ Standard Library only. No external optimization libraries were used.

---

## 🚀 How to Compile and Run

### Compilation
You can compile the source code using a standard C++ compiler like g++.

```bash
g++ -std=c++11 -o solver hw2.cpp
```

### Execution
The program is run from the command line with the following arguments:

```bash
./solver [instance-file-path] [max-execution-time] [max-evaluations]
```
* `instance-file-path`: The path to the problem instance file (e.g., `instances/C101.txt`).
* `max-execution-time`: The maximum run time in seconds. Use `0` for no time limit.
* `max-evaluations`: The maximum number of objective function evaluations. Use `0` for no limit.

**Example:**
```bash
./solver instances/C101.txt 60 0
```
This command runs the solver on the `C101.txt` instance for a maximum of 60 seconds.

---

## 📄 License
This project is licensed under the MIT License.