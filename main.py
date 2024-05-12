import math
import numpy as np
import sys
import os

from time import time # Used for debugging purposes

t1 = time()

def main(input_location = [], output_location = []):
    class Node:
        def __init__(self, data):
            """
            Instantiation method to create a new "Node" object to be used in matrix calculations

            Parameters  :
            data        : list -> List of data provided, formatted as ["n1=x", "n2=y", "R=z"] where x and y are integers and z is a float
            
            Returns nothing
            """
            """
            >>> node = Node(["n1=0", "n2=1", "R=1"])
            >>> node.n1_value
            1
            >>> node.n2_value
            0
            >>> node.component_type
            'R'
            >>> node.component_value
            1.0
            >>> node = Node(["n1=a", "n2=b", "R=1"])
            Traceback (most recent call last):
            ...
            TypeError: Node values must be integers
            >>> node = Node(["n1=0", "n2=0", "R=1"])
            Traceback (most recent call last):
            ...
            ValueError: Node values cannot be the same
            >>> node = Node(["n1=-1", "n2=1", "R=1"])
            Traceback (most recent call last):
            ...
            ValueError: Node value cannot be negative
            >>> node = Node(["n1=0", "n2=1", "Z=1"])
            Traceback (most recent call last):
            ...
            ValueError: Invalid component type
            >>> node = Node([])
            Traceback (most recent call last):
            ...
            Exception: No data provided to Node Class
            """
            if not data: raise Exception("No data provided to Node Class") # Checks if data is empty
            self.impedances = []
            self.n1 = str(data[0]).split("=")
            # Sets up the first node, direction is not decided yet
            self.n2 = str(data[1]).split("=")
            # Sets up the second node, direction is not decided yet
            self.component = str(data[2]).split("=")
            # Splits the component into the type, seen below, and value, seen later
            self.component_type = self.component[0].upper()
            if self.component_type not in ["R", "G", "C", "L"]: raise ValueError("Invalid component type")
            # Sets the component type and checks if it is valid

            try:
                start = int(self.n1[1])
                end = int(self.n2[1])
                # Creates two temporary variables to store the node values so they can be swapped if needed
                if start < 0 or end < 0: raise ValueError("Node value cannot be negative")
                if start == end: raise ValueError("Node values cannot be the same")
                if start == 0 and end > 0:
                    # Checks if a node is 0 or if the first node is greater than the second
                    self.n1_value = end
                    self.n2_value = start
                else:
                    # Does not swap the nodes
                    self.n1_value = start
                    self.n2_value = end
                # Swaps the nodes so A -> B will always have B = 0 or B > A
                self.component_value = float(self.component[1])
                self.type = "series"
                if not self.n2_value: self.type = "shunt"
            except TypeError:
                raise Exception("Invalid data type parsed as node value, node value must be integer")

    class Source:
        def __init__(self, source_parameters: list = []) -> None:
            """
            Instantiation method to create a new "Source" object to be used in matrix calculations

            Parameters  :
            source_parameters : list -> List of source parameters provided, formatted as [["VT=1", "RS=1"], ["RL=1"], ["FSTART=0", "FEND=1", "NFREQS=100"]]

            Returns nothing
            """
            if not source_parameters: raise Exception("No source parameters provided") # Checks if source parameters are empty
            try:
                self._source: list = [item.split("=") for item in source_parameters[0]]
                self._load: list = [item.split("=") for item in source_parameters[1]]
                self._sim: list = [item.split("=") for item in source_parameters[2]]

                self.source_type: str = self._source[0][0]
                self.source_value: int = float(self._source[0][1])

                self.source_res: str = self._source[1][0]
                self.source_res_value: float = float(self._source[1][1])
                if self.source_res.__contains__("G"): self.source_res = float(self._source[1][1])

                self.load_type: str = self._load[0][0]
                self.load_value: float = float(self._load[0][1])            

                self.sim_start: float = float(self._sim[0][1])
                self.sim_end: float = float(self._sim[1][1])
                self.sim_freq: int = int(self._sim[2][1])
                if self.sim_start == "LFSTART": self.np_sim = np.logspace(self.sim_start, self.sim_end, self.sim_freq)
                else: self.np_sim = np.linspace(self.sim_start, self.sim_end, self.sim_freq)
                self.column = np.linspace(self.sim_start, self.sim_end, self.sim_freq)
                

                self._contents = [self.source_type, self.source_value,
                            self.source_res, self.source_res_value,
                            self.load_type, self.load_value,
                            self.sim_start, self.sim_end, self.sim_freq]
            except: raise Exception("Invalid source input array, missing data")
        def _debug(self): return self._contents # Use only to test bulk value reading

    objs = []
    nodes = []
    tags = {"<CIRCUIT>", "</CIRCUIT>", "<TERMS>", "</TERMS>", "<OUTPUT>", "</OUTPUT>", "#", "<", ">"}
    # List of tags that occur in .net files

    src = {"VT", "IN", "NFREQS", "FEND", "FSTART", "LFSTART", "LFEND", "RL", "GS"}
    # List of source parameters that can appear in a .net file, log frequencies included

    with open(input_location, "rt") as net_file: contents = net_file.read().splitlines()
    # After testing this is the fastest way to read the file

    source_values = []
    # List of source values that will be used to create the source object

    for line in contents:
        if not line: continue
        if not set(line).intersection(tags):    # Creates an intersection between the line and the tags, if its blank it continues
            if line.__contains__("n1"):         # Checks if the line contains n1, if it does it is a node
                nodes.append(Node(line.split(" ")))
                continue
            if set(line.upper().replace("="," ").split(" ")).intersection(src):
                source_values.append(line.split(" "))
            else: objs.append(line.split(" "))
    # After testing this is the fastest way to split the lines

    nodes = sorted(nodes, key=lambda x: (x.n1_value, x.n2_value)) # Sorts the nodes into A -> 0, A -> B, B -> 0, etc.
    src = Source(source_values)
    # Creates a new source object containing all information pertaining to the source

    matrices = []
    # Blank list for the combined matrices

    previous_n1, previous_n2 = 0, 0

    for node in nodes:
        # Calculates impedances for each node connection for each frequency using a linear algebra approach
        match node.component_type:
            case "R":
                node.impedances.append([node.component_value] * len(src.np_sim))
            case "G":
                node.impedances.append([1/node.component_value] * len(src.np_sim))
            case "C":
                node.impedances.append(1/(2*math.pi*1j*node.component_value*src.np_sim))
            case "L":
                node.impedances.append(2*math.pi*1j*node.component_value*src.np_sim)
            case _:
                raise ValueError("Invalid component type")
                # Should not occur as the component type is checked in the Node class

    for node in nodes:
        node.impedances = node.impedances[0] # Formats the lists to be correct

    for node in nodes:
        if node.n1_value == previous_n1 and node.n2_value == previous_n2:
            new_impedance = []
            recip_1 = np.multiply(1, 1/(np.array(last_impedance, dtype=np.clongdouble)))
            recip_2 = np.multiply(1, 1/(np.array(node.impedances, dtype=np.clongdouble)))
            for value in (1/(recip_1 + recip_2)):
                new_impedance.append(value)
            matrices[-1] = [node.n1_value, node.n2_value, new_impedance]
            last_impedance = new_impedance
            previous_n1 = node.n1_value
            previous_n2 = node.n2_value
        else:
            previous_n1 = node.n1_value
            previous_n2 = node.n2_value
            last_impedance = np.array(node.impedances, dtype=np.clongdouble)
            matrices.append([node.n1_value, node.n2_value, node.impedances])

    final_matrix = []
    for count in range(len(src.np_sim)):
        temp_matrix = np.eye(2, dtype=complex)
        for mat in matrices:
            direction = mat[1] > 0
            if direction:
                # [1 Z][0 1] if dir else [1 0][Z 1]
                new_matrix = np.array([[1, mat[2][count]], [0, 1]], dtype=np.clongdouble)
            else:
                new_matrix = np.array([[1, 0], [1/mat[2][count], 1]], dtype=np.clongdouble)
            temp_matrix = np.matmul(temp_matrix, new_matrix)
        final_matrix.append(temp_matrix)

    matrices = final_matrix
    # Puts final_matrix into matrices, to iterate through again
    # Solely used for formatting purposes

    class output_data:
        # Class to store output data using real and imaginary column
        def __init__(self):
            self.real = []
            self.imag = []

    output = {
        "Freq": output_data(),
        "Vin": output_data(),
        "Vout": output_data(),
        "Iin": output_data(),
        "Iout": output_data(),
        "Pin": output_data(),
        "Pout": output_data(),
        "Zin": output_data(),
        "Zout": output_data(),
        "Av": output_data(),
        "Ap": output_data(),
        "Ai": output_data()
    } # All headers seen in .csv examples

    output["Freq"].real = src.np_sim

    zin_values = []
    zout_values = []
    vin_values = []
    vout_values = []
    iin_values = []
    iout_values = []
    pin_values = []
    pout_values = []
    av_values = []
    ap_values = []
    ai_values = []

    for matrix in matrices:
        A = matrix[0][0]
        B = matrix[0][1]
        C = matrix[1][0]
        D = matrix[1][1]
        # Creates an ABCD matrix for each frequency

        current_zin = (A*src.load_value + B) / (C*src.load_value + D)
        if src.source_type == "IN":
            current_vin = current_zin * src.source_value
        else:
            current_vin = src.source_value * (current_zin * 1/(src.source_res_value + current_zin))
        current_cin = current_vin / current_zin

        current_voltage_gain = (src.load_value / (A * src.load_value + B))
        current_vout = current_vin * current_voltage_gain
        current_cgain = 1/(C * src.load_value + D)
        current_cout = current_cgain * current_cin
        conj = np.conjugate(current_cgain)
        current_pgain = current_voltage_gain * conj
        current_pin = current_vin * np.conjugate(current_cin)
        current_zout = (D * src.source_res_value + B) / (C * src.source_res_value + A)
        current_pout = current_pgain * current_pin



        zin_values.append(current_zin)
        zout_values.append(current_zout)
        vin_values.append(current_vin)
        vout_values.append(current_vout)
        iin_values.append(current_cin)
        iout_values.append(current_cout)
        pin_values.append(current_pin)
        pout_values.append(current_pout)
        av_values.append(current_voltage_gain)
        ap_values.append(current_pgain)
        ai_values.append(current_cgain)
        # Calculates values according to spec sheet then appends them to the output

    order = [["Freq", "Hz"]] # Keep track of the order the csv file will be written in

    for obj in objs:
        """
        This loop will iterate through each object in the objs list and calculate log and linear values for each object
        It will also split objects into header and unit, i.e [["Vin", "dBV"]] is split into "Vin" and "dBV"
        Also splits each value into its real and imaginary components, with "decible" being used as a boolean flag for logarithmic calculations
        """
        header = obj[0]
        unit = None
        decibel = False
        if len(obj) > 1:
            decibel = obj[1].__contains__("dB")
            unit = obj[1]
        order.append([header, unit])
        match header:
            case "Vin":
                value = 20*np.log10(np.absolute(vin_values)) if decibel else vin_values
                real = np.real(value) # This will always be used
                imag = np.imag(value) # This will be used if ¬decibel
                angle = np.angle(vin_values) # only use if decibel otherwise ignore
                output["Vin"].real.append(real)
                output["Vin"].imag.append(angle if decibel else imag)
                continue
            case "Vout":
                value = 20*np.log10(np.absolute(vout_values)) if decibel else vout_values
                real = np.real(value)
                imag = np.imag(value)
                angle = np.angle(vout_values)
                output["Vout"].real.append(real)
                output["Vout"].imag.append(angle if decibel else imag)
                continue
            case "Iin":
                value = 20*np.log10(np.absolute(iin_values)) if decibel else iin_values
                real = np.real(value)
                imag = np.imag(value)
                angle = np.angle(iin_values)
                output["Iin"].real.append(real)
                output["Iin"].imag.append(angle if decibel else imag)
                continue
            case "Iout":
                value = 20*np.log10(np.absolute(iout_values)) if decibel else iout_values
                real = np.real(value)
                imag = np.imag(value)
                angle = np.angle(iout_values)
                output["Iout"].real.append(real)
                output["Iout"].imag.append(angle if decibel else imag)
                continue
            case "Pin":
                value = 10*np.log10(np.absolute(pin_values)) if decibel else pin_values
                real = np.real(value)
                imag = np.imag(value)
                angle = np.angle(pin_values) 
                output["Pin"].real.append(real)
                output["Pin"].imag.append(angle if decibel else imag)
                continue
            case "Pout":
                value = 10*np.log10(np.absolute(pout_values)) if decibel else pout_values
                real = np.real(value)
                imag = np.imag(value)
                angle = np.angle(pout_values)
                output["Pout"].real.append(real)
                output["Pout"].imag.append(angle if decibel else imag)
                continue
            case "Zin":
                value = 10*np.log10(np.absolute(zin_values)) if decibel else zin_values
                real = np.real(value)
                imag = np.imag(value)
                angle = np.angle(zin_values)
                output["Zin"].real.append(real)
                output["Zin"].imag.append(angle if decibel else imag)
                continue
            case "Zout":
                value = 10*np.log10(np.absolute(zout_values)) if decibel else zout_values
                real = np.real(value)
                imag = np.imag(value)
                angle = np.angle(zout_values)
                output["Zout"].real.append(real)
                output["Zout"].imag.append(angle if decibel else imag)
                continue
            case "Av":
                value = 20*np.log10(np.absolute(av_values)) if decibel else av_values
                real = np.real(value)
                imag = np.imag(value)
                angle = np.angle(av_values)
                output["Av"].real.append(real)
                output["Av"].imag.append(angle if decibel else imag)
                continue
            case "Ap":
                value = 20*np.log10(np.absolute(ap_values)) if decibel else ap_values
                real = np.real(value)
                imag = np.imag(value)
                angle = np.angle(ap_values) 
                output["Ap"].real.append(real)
                output["Ap"].imag.append(angle if decibel else imag)
                continue
            case "Ai":
                value = 20*np.log10(np.absolute(ai_values)) if decibel else ai_values
                real = np.real(value)
                imag = np.imag(value)
                angle = np.angle(ai_values)
                output["Ai"].real.append(real)
                output["Ai"].imag.append(angle if decibel else imag)
                continue
            case _:
                raise ValueError("Invalid header")

    output["Freq"] = np.array([output["Freq"].real])

    headers = []
    units = []

    for i in range(len(order)):
        """
        This loop iterates through each header and unit and will append them to the output .csv file
        The two loops are separate as soluble values have their own formatting needed
        """
        header = order[i][0]
        unit = order[i][1]
        if not unit: unit = "L"
        if header == "Freq":
            headers.append("{:>11}".format(header + ","))
            units.append("{:>11}".format(unit + ","))
            continue
            # Frequency is a special case
        # These if statements can be combined into one, but for readability they are separate
        if unit.__contains__("dB"):
            real_header = "{:>12}".format("|" + header + "|,")
            imag_header = "{:>12}".format("/_" + header + ",")
            headers.append(real_header)
            headers.append(imag_header)
            real_unit = "{:>12}".format(unit + ",")
            imag_unit = "{:>12}".format("Rads,")
            units.append(real_unit)
            units.append(imag_unit)
        else:
            real_header = "{:>12}".format("Re(" + header + "),")
            imag_header = "{:>12}".format("Im(" + header + "),")
            headers.append(real_header)
            headers.append(imag_header)
            real_unit = "{:>12}".format(unit + ",")
            imag_unit = "{:>12}".format(unit + ",")
            units.append(real_unit)
            units.append(imag_unit)

    row_count = len(output["Freq"].real[0])

    buffer = ""

    with open(output_location, "wt") as csv_file:
        """
        - Opens the file
        - Iterates through the output data
        - Appends the formatted real and imaginary components to the row
        - Writes the row to the file
        - Repeats until all data is written
        """
        csv_file.write("".join(headers)[:-1] + "\n")
        csv_file.write("".join(units)[:-1] + "\n")
        for index in range(row_count):
            for header in range(len(order)):
                output_data = order[header][0]
                data_real = output[output_data].real
                data_imag = output[output_data].imag
                if output_data == "Freq":
                    frequency = "{:.3e}".format(data_real[0][index])
                    buffer += ("{:>11}".format(frequency + ",")) #csv_file.write("{:>11}".format(frequency + ","))
                    continue
                else:
                    if data_real == [] or data_imag == []:
                        continue
                    data_real = "{:.3e}".format(data_real[0][index])
                    data_imag = "{:.3e}".format(data_imag[0][index])
                    buffer += ("{:>12}".format(data_real + ","))
                    buffer += ("{:>12}".format(data_imag + ","))
                    #csv_file.write("{:>12}".format(data_real + ","))
                    #csv_file.write("{:>12}".format(data_imag + ","))
            buffer += "\n" #csv_file.write("\n")
        csv_file.write(buffer)


if __name__ == "__main__":
    if len(sys.argv) != 3: raise Exception("Incorrect number of arguments provided.")
    infile = sys.argv[1]
    outfile = sys.argv[2]

    main(infile, outfile)
    print(time() - t1)
