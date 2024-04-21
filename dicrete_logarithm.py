from sympy import factorint  
import math

def brute_force(a, b, p): 
    x = 0
    t = 1

    for i in range(0, p):
        if t == b:
            return x

        t = (t * a) % p
        x = x + 1
    
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

if __name__ == "__main__":
    a = 44573028072
    b = 32221422123
    p = 78624916757

    x = silver_pohlig_hellman(a, b, p)
    print(x)

    x = brute_force(a, b, p)
    print(x)