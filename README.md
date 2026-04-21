# numerical-calculator-
this is a Python-based desktop application for solving complex numerical equations and linear algebra systems. It features a modern, dark-themed graphical user interface (GUI) built with `customtkinter`, plotting capabilities with `matplotlib`, and analytical math parsing handled seamlessly via `sympy`. 

## 🚀 Features
- **Root Finding**: Solves roots using Bisection, False Position, Fixed Point Iteration, Newton-Raphson, and Secant methods.
- **Linear Systems**: Evaluates custom Augmented Matrices (3x4) utilizing Gauss-Jordan Elimination, LU Decomposition, and Cramer's Rule.
- **Dynamic Charting**: Visually plots the function curves and intercepts in real-time.
- **Report Exporting**: Native capability to export iterative computation data directly to Excel spreadsheets (`.xlsx`) or PDF reports.

## 📦 Prerequisites

Ensure you have Python 3.x installed. Then, install the required packages using the following command:

```bash
pip install customtkinter pandas fpdf numpy matplotlib sympy
```

## 🛠️ Usage

Simply run the graphical interface file from your terminal:

```bash
python gui.py
```
*(Alternatively, you can run `python main.py` if you prefer the text-based CLI terminal menu).*

---

## 📂 Codebase Explanation

The logic is separated intuitively into two core modules:

### 1. `main.py` (The Mathematical Backend)
This file contains the underlying logic and mathematical algorithms without any graphical interface.
- **Root Finding Algorithms**: Includes algorithms like `bisect`, `false_position`, `fixedpoint`, and `newton` to computationally hunt for the `f(x)` root.
- **Linear Systems Algorithms**: Implements functions to tackle $3 \times 4$ matrices using `gauss_jordan_elimination`, `lu_decomposition`, and `cramers_rule`.

### 2. `gui.py` (The Visual Interface)
This file handles user input, error parsing, chart generation, and data formatting.
- **`NumericalApp` Class**: Initializes the `customtkinter` main window, sets up the purple/dark-mode theme palettes, and handles the `tkinter` tabs layout. 
- **Dynamic Equations** (`_parse_equation`): Uses the `sympy` library to convert a physical string (e.g., `"x**2 - 4"`) into an executable math function. **It is also responsible for automatically generating the derivative (slope) of the equation!**
- **Calculation Execution** (`solve_root` & `solve_linear`): Takes the parsed properties and connects them directly with the logic implemented inside `main.py`.
- **Plotting Integration** (`_update_plot`): Ties a `FigureCanvasTkAgg` onto the Tkinter frame to map the results of the equations over a Cartesian plane continuously.
- **Data Exporting**: Hooks up `fpdf` and `pandas` to easily allow the user to save table outcomes off-screen.
