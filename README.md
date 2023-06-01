# Verification of Quantum Systems using Barrier Certificates

Artifact for "Verification of Quantum Systems using Barrier Certificates".

## Running the Code

The examples are tested on a virtual machine using Python 3.9.5.

1. Install the requirements file (if in the provided VM you only need to run (b) to activate the environment)

    a. ```python -m venv env```

    b. ```source env/bin/activate```

    c. ```pip install -r requirements.txt```
    
2. Run the examples:

    ```python ex_<name>.py``` 

3. Generating the figures:
    
    ```python media/diagrams.py```

    Figures are available in ```media/Images/```.

## Project Navigation
```
├── media
|   ├── Animations - folder where generated animations are stored in
|   |   ├── .animations
|   ├── Images - folder where generated images are stored
|   |   ├── .images
|   | animations.py - generates bonus animations
|   | diagrams.py - generates Figures 2-3 using matplotlib 
|   | utils_plot.py - utilities for generating Figures
├── src
│   ├── __init__.py
│   ├── find_barrier.py - implementation of Algorithm 1 using SciPy
│   ├── FuncClasses.py - classes for representing barriers with arbitrary coefficients
│   ├── sympy2rep.py - conversion from Sympy symbols (and lists) to classes defined in FuncClasses
│   ├── utils.py - utilities for other files in src
├── ex_cnot.py - Generates functions in Section 5.3
├── ex_hadmard.py - Generates functions in Section 5.1
├── ex_phase.py - Generates functions in Section 5.2
├── LICENSE
├── README.md
├── requirements.txt - Required packages for environment
└── .gitignore
```
