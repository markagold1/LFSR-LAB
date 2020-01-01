import numpy as np

def boxmuller(u):
    M = u.size
    if M % 2:
        u = np.append(u,u[-1])

    if max(u) >= 1:
        denom = 2**np.float64(np.ceil(np.log2(np.max(u))))
        u = u / denom

    eps = np.finfo(float).eps
    u[u==0] = eps
    u1 = u[0::2]
    u2 = u[1::2]

    n1 = np.multiply(np.sqrt(-2*np.log(u1)),np.cos(2*np.pi*u2))
    n2 = np.multiply(np.sqrt(-2*np.log(u1)),np.sin(2*np.pi*u2))

    n = np.hstack([n1,n2])
    if M - u.size:
      n = n[:-1]

    
    return (u,n)


def boxmuller_example(N=65536):
    M = N
    if M % 2:
        M += 1

    u = np.random.uniform(0,1,M)
    u1 = u[0::2]
    u2 = u[1::2]

    n1 = np.multiply(np.sqrt(-2*np.log(u1)),np.cos(2*np.pi*u2))
    n2 = np.multiply(np.sqrt(-2*np.log(u1)),np.sin(2*np.pi*u2))

    n = np.hstack([n1,n2])
    if M - N:
      n = n[:-1]

    return (u,n)

