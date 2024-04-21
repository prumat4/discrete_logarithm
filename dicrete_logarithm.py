import time 

def brute_force(a, b, p): 
    x = 0
    t = 1

    for i in range(0, p):
        if t == b:
            return x

        t = (t * a) % p
        x = x + 1
    
    return float('nan')


if __name__ == "__main__":
    a = 3
    b = 4
    p = 13

    x = brute_force(a, b, p)
    print(x)