import matplotlib.pyplot as plt
from sympy import factorint
import numpy as np
import argparse
import pexpect
import time
import os 

def brute_force(a, b, p): 
    x = 0
    t = 1

    for i in range(0, p):
        if t == b:
            return x

        t = (t * a) % p
        x += 1
    
    return float('nan')

def chinese_remainder_theorem(remainders, modules):
    product = 1
    for module in modules:
        product *= module
    
    result_sum = 0
    for index, module in enumerate(modules):
        partial_product = product // module
        modular_inverse = pow(partial_product, -1, module)
        result_sum += remainders[index] * modular_inverse * partial_product

    return result_sum % product

def silver_pohlig_hellman(base, exponent, prime):
    phi = prime - 1
    prime_factors = factorint(phi)
    residue_table = {}

    for prime_factor, exp_power in prime_factors.items():
        for power in range(prime_factor):
            residue = pow(base, (phi * power) // prime_factor, prime)
            residue_table[prime_factor, residue] = power

    x_results = []
    factor_exponents = [prime_factor ** exp_power for prime_factor, exp_power in prime_factors.items()]

    for prime_factor, exp_power in prime_factors.items():
        mod_exponent = prime_factor ** exp_power
        modular_x = 1
        accumulated_x = 0
        for degree in range(exp_power):
            term = pow(exponent * pow(modular_x, -1, prime), phi // (prime_factor ** (degree + 1)), prime)
            x_component = residue_table.get((prime_factor, term), 0)
            accumulated_x += x_component * (prime_factor ** degree)
            modular_x = (modular_x * pow(base, x_component * (prime_factor ** degree), prime)) % prime
        x_results.append(accumulated_x % mod_exponent)

    return chinese_remainder_theorem(x_results, factor_exponents)

def initiate_docker(digit_length):
    child = pexpect.spawn("docker run -it salo1d/nta_cp2_helper:2.0")
    child.expect("Enter prime number decimal digit length:")
    child.sendline(str(digit_length))
    return child

def get_test_numbers(child):
    child.expect("a = (\d+);")
    a = int(child.match.group(1))
    child.expect("b = (\d+);")
    b = int(child.match.group(1))
    child.expect("p = (\d+).")
    p = int(child.match.group(1))
    return a, b, p

def test_sph_algorithm(digit_length, results):
    for i in range(2): 
        child = initiate_docker(digit_length)
        a, b, p = get_test_numbers(child)
        print(f"\nTesting SPH with a = {a}, b = {b}, p = {p}")

        start_time = time.time()
        x_sph = silver_pohlig_hellman(a, b, p)
        sph_duration = time.time() - start_time
        print(f"digit length: {digit_length}, task type: {i + 1}, SPH result: x = {x_sph} (took {sph_duration:.2f} seconds)")
        
        results.append((digit_length, i+1, sph_duration))
        
        child.close()

def test_brute_force_algorithm(digit_length, results):
    for i in range(2): 
        child = initiate_docker(digit_length)
        a, b, p = get_test_numbers(child)
        print(f"\nTesting Brute Force with a = {a}, b = {b}, p = {p}")

        start_time = time.time()
        x_bf = brute_force(a, b, p)
        bf_duration = time.time() - start_time
        print(f"digit length: {digit_length}, task type: {i + 1}, Brute Force result: x = {x_bf} (took {bf_duration:.2f} seconds)")
        
        results.append((digit_length, i+1, bf_duration))
        
        child.close()

def run_tests():
    n = 15
    sph_results = []  
    bf_results = []  
    for digit_length in range(3, n):
        test_sph_algorithm(digit_length, sph_results)
    
    for digit_length in range(3, 8):
        test_brute_force_algorithm(digit_length, bf_results)
    
    digit_lengths = range(3, n)
    
    sph_times_type1 = [np.mean([result[2] for result in sph_results if result[0] == dl and result[1] == 1]) for dl in digit_lengths]
    sph_times_type2 = [np.mean([result[2] for result in sph_results if result[0] == dl and result[1] == 2]) for dl in digit_lengths]
    bf_times_type1 = [np.mean([result[2] for result in bf_results if result[0] == dl and result[1] == 1]) for dl in digit_lengths]
    bf_times_type2 = [np.mean([result[2] for result in bf_results if result[0] == dl and result[1] == 2]) for dl in digit_lengths]
    
    n_groups = len(digit_lengths)
    bar_width = 0.2 
    
    fig, ax = plt.subplots()
    
    index = np.arange(n_groups)
    
    bar1 = ax.bar(index - bar_width, sph_times_type1, bar_width, alpha=0.8, color='darkgreen', label='SPH Type 1')
    bar2 = ax.bar(index, sph_times_type2, bar_width, alpha=0.8, color='lightgreen', label='SPH Type 2')
    bar3 = ax.bar(index + bar_width, bf_times_type1, bar_width, alpha=0.8, color='navy', label='Brute Force Type 1')
    bar4 = ax.bar(index + 2*bar_width, bf_times_type2, bar_width, alpha=0.8, color='skyblue', label='Brute Force Type 2')
    
    ax.set_xlabel('Digit Length')
    ax.set_ylabel('Time Taken (seconds)')
    ax.set_title('Performance Comparison by Algorithm and Type')
    ax.set_xticks(index + bar_width / 2)
    ax.set_xticklabels([str(dl) for dl in digit_lengths])
    ax.legend()
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run different types of tests on algorithms.")
    parser.add_argument('--tests', action='store_true', help='Run the testing suite')
    args = parser.parse_args()

    if args.tests:
        run_tests()
    else:
        a = 44573028072
        b = 32221422123
        p = 78624916757

        x = silver_pohlig_hellman(a, b, p)
        print(f"Silver-Pohlig-Hellman result: x = {x}")
        
        x = brute_force(a, b, p)
        print(f"Brute force: x = {x}")
