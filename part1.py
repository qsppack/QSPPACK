# Importing standard Qiskit libraries and configuring account
import qiskit as qk
from qiskit import QuantumCircuit, Aer, IBMQ
from qiskit import transpile, assemble
from qiskit.tools.jupyter import *
from qiskit.visualization import plot_histogram
from math import pi
import numpy as np
from qiskit import QuantumRegister, ClassicalRegister
from PIL import Image
import tensorflow as tf

def binary_encode(image, threshold=1e-5):
    # pixel value is 1 if it's greater than threshold or else zero
    encoded_image = [1 if j>threshold else 0 for j in image[0]]
    return np.array([encoded_image])

def encode(image):
# Normalize image values between 0 and 1
    image = image/255.0

# Reshape the image
    image = image.reshape(image.shape[0], *(28,1))
    image = tf.image.resize(image, (2,2)).numpy()
    
    idx = QuantumRegister(2, 'idx')
# grayscale pixel intensity value
    intensity = QuantumRegister(4,'intensity')
# classical register
    cr = ClassicalRegister(6, 'cr')

# create the quantum circuit for the image
    qc_image = QuantumCircuit(intensity, idx, cr)

# set the total number of qubits
    num_qubits = qc_image.num_qubits
    for idx in range(intensity.size):
        qc_image.i(idx)

# Add Hadamard gates to the pixel positions    
    qc_image.h(4)
    qc_image.h(5)

# Separate with barrier so it is easy to read later.
    qc_image.barrier()
    bin_val = binary_encode(image)
    # Use join() to concatenate the elements of the array into a string
    result_string = ''.join(map(str, bin_val))
    value01 = result_string

# Add the NOT gate to set the position at 01:
    qc_image.x(qc_image.num_qubits-1)

# We'll reverse order the value so it is in the same order when measured.
    for idx, px_value in enumerate(value01[::-1]):
        if(px_value=='1'):
            qc_image.ccx(num_qubits-1, num_qubits-2, idx)

# Reset the NOT gate
    qc_image.x(num_qubits-1)

    qc_image.barrier()
    return qc_image

def decode(histogram):
    hist_val = list(histogram.keys())
    binary_values = hist_val
    # reshape to 2x2 array
    binary_values = np.array(binary_values, dtype=np.uint8)
    image = binary_values.reshape((2, 2))
    plt.imshow(image)
    return plt.show()

def simulate(qc_image):
    qc_image.measure(range(6),range(6))
    aer_sim = Aer.get_backend('aer_simulator')
    t_qc_image = transpile(qc_image, aer_sim)
    qobj = assemble(t_qc_image, shots=1000)
    job_neqr = aer_sim.run(qobj)
    result_neqr = job_neqr.result()
    counts_neqr = result_neqr.get_counts()
    return counts_neqr

def run_part1(image):
    #encode image into a circuit
    circuit=encode(image)

    #simulate circuit
    histogram=simulate(circuit)

    #reconstruct the image
    image_re=decode(histogram)

    return circuit,image_re