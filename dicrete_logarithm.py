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

def chinese_remainder_theorem(remainders, modulus):
    product_of_modulus = 1
    for mod in modulus:
        product_of_modulus *= mod

    total = 0
    for i, (remainder, modulus) in enumerate(zip(remainders, modulus)):
        partial_product = product_of_modulus // modulus
        modular_inverse = pow(partial_product, -1, modulus)
        contribution = remainder * modular_inverse * partial_product
        total += contribution

    return total % product_of_modulus


def silver_pohlig_hellman(base, target, prime):
    phi = prime - 1
    prime_factors = factorint(phi)
    exponent_table = {}

    for factor, count in prime_factors.items():
        exponent_table.update(
            {(factor, pow(base, (phi * exp) // factor, prime)): exp
             for exp in range(factor)}
        )

    results = []
    base_powers = [factor ** count for factor, count in prime_factors.items()]

    for factor, exponent in prime_factors.items():
        current = 1
        partial_result = 0
        for exponent_index in range(exponent):
            power = phi // (factor ** (exponent_index + 1))
            adjusted_target = pow(target * pow(current, -1, prime), power, prime)
            partial_exponent = exponent_table.get((factor, adjusted_target), 0)
            partial_result += (partial_exponent * (factor ** exponent_index)) % (factor ** exponent)
            current = (current * pow(base, partial_exponent * (factor ** exponent_index), prime)) % prime
        results.append(partial_result)

    return chinese_remainder_theorem(results, base_powers)

if __name__ == "__main__":
    a = 798831748718
    b = 323490752685
    p = 818066545501

    x = silver_pohlig_hellman(a, b, p)
    print(x)

    x = brute_force(a, b, p)
    print(x)