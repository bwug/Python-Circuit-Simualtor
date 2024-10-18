# Cascade Circuit Simulator in Python 3.11

### Context & Information

In electronics, a cascade circuit is a circuit that connects the output of a circuit block or component to the input of the next in a linear sequence, similar to pipelining in computer science.

This python script will take a pre-written cascade circuit and provide a calculated output. At the moment the circuit only supports the passive components: R(esistor), L(Inductor), C(apacitor), A(dmitter)

### Usage guide

The code is entirely run using the `sys.argv` library in the python stdlib, with *numpy* being used for calculation.

If intended to run this through a venv, ensure you run `pip install -r requirements.txt` or install the most recent version of numpy.

To run, use the command `python3 input_circuit.txt output_circuit.csv`

### Testing

Thorough testing was performed to ensure robustness and speed of execution, the timing code was kept to assure the user the code completes, with the average time being sub 1-second.

Testing was performed using in-built libraries and a testing table, these are kept hidden due to this being partially done as university work.