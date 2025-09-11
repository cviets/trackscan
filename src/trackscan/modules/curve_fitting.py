import numpy as np
from numpy.linalg import inv
from typing import List, Tuple
from numpy.typing import NDArray

def get_lsq_solution(n: int, f: int, t: List) -> Tuple[List, NDArray]:
    """
    Returns least squares n-degree polynomial fit to y = f(t).
    Parameters
    ----------
    n : int
        polynomial degree 
    f : list
        y values
    t : list
        x values

    Returns
    -------
    beta : array
        fitted parameters
    cov : 2d array
        covariances of parameters

    """
    A = np.array([[j**i for i in range(n, -1, -1)] for j in t])
    AT = np.transpose(A)
    cov = inv(np.matmul(AT, A))
    beta = np.matmul(np.matmul(cov, AT), f)
    beta = list(reversed(beta))
    return beta, cov

def least_squares(n: int, f: List[int], t: List[int], return_error: bool=False):
    """
    fit f, t to a n-degree polynomial. 
    Returns best-fit polynomial as a function with input x
    """
    
    A = np.array([[j**i for i in range(n, -1, -1)] for j in t])
    AT = np.transpose(A)
    cov = inv(np.matmul(AT, A))
    beta = np.matmul(np.matmul(cov, AT), f)
    beta = list(reversed(beta))
    
    beta, _ = get_lsq_solution(n, f, t)
    
    return lambda x: sum(beta[i]*x**i for i in range(len(beta)))

def least_squares_error(n, f, t):
    """
    Returns functions corresponding to +/- 1 sigma in the error of the least-squares
    n-degree polynomial fit to y = f(t). In other words, f(t) should fit within inf(t)
    and sup(t) most of the time. May not be super reliable.
    
    inf: function corresponding to 1 sigma below the best fit. inf(t) < fit(t) < sup(t)
    sup: function corresponding to 1 sigma above best fit.
    """
    
    beta, cov = get_lsq_solution(n, f, t)
    F = lambda x: sum(beta[i]*x**i for i in range(len(beta)))
    
    S = sum((F(T) - f[i])**2 for i, T in enumerate(t))
    var_beta = list(reversed([cov[i][i]*S/(len(f)-n) for i in range(len(beta))]))
    sup = lambda x: sum((beta[i]+np.sqrt(var_beta[i]))*x**i 
                        for i in range(len(beta)))
    inf = lambda x: sum((beta[i]-np.sqrt(var_beta[i]))*x**i 
                        for i in range(len(beta)))
    
    return inf, sup

def fit_recent(f, t):
    """
    Returns a list that fits y = f(t) based on its previous values. Fit at time 
    i is determined by fitting a parabola to the previous four values 
    (f(t=i-5) to f(t=i-1)). I.e., fits f(t) using its most recent values.
    """
    
    assert len(f) >= 5, "f must have length >= 5 for accurate predictions"
    prediction = f[0:4]
    for i in range(len(f)-4):
        sub_f = f[i:i+4]
        sub_t = t[i:i+4]
        F = least_squares(2, sub_f, sub_t)
        prediction.append(F(t[i+4]))
        
    return prediction
    

def find_breakpoints(f, t):
    """
    Finds where f(t) changes abruptly, as determined by when it deviates significantly
    from fit determined by fit_recent.
    
    Returns idx: set of indices where f changes abruptly.
    """
    F = fit_recent(f, t)
    residuals = [F[i]-f[i] for i in range(len(f))]
    for i in range(4, len(residuals)):
        if residuals[i] == 0:
            return set()
    res2 = [r**2 for r in residuals]
    partial_sums = [res2[i] + res2[i+1] + res2[i+2] for i in range(len(res2) - 2)]
    ssr = sum(res2)
    expected_res2 = ssr/len(res2)
    exp_partial_sums = [0, 0, expected_res2, expected_res2*2]
    exp_partial_sums.extend([expected_res2*3 for _ in range(len(partial_sums)-4)])
    # standard deviation of each residual:
    rmsd = np.sqrt(ssr/(len(t) - 2))
    # standard deviations of the squared residuals:
    d_res2 = [2*abs(x)*rmsd for x in residuals]
    d_partial_sums = [np.sqrt(d_res2[i]**2 + d_res2[i+1]**2 + d_res2[i+2]**2) for i in range(len(res2) - 2)]
    z_scores = [(partial_sums[i] - exp_partial_sums[i])/d_partial_sums[i] for i in range(2, len(partial_sums))]
    res2_z_scores = [(res2[i] - expected_res2)/d_res2[i] for i in range(4, len(res2))]
    idx = {i+1 for i in range(2, len(z_scores)) if z_scores[i]>2 and res2_z_scores[i-2] > 2}
    
    removals = {i for i in idx if i-1 in idx or i-2 in idx}
    idx -= removals
    
    return idx

def get_squared_distance_between_curves(x, f, g):
    """
    For array-like x and functions f, g of x, returns sum([(f(x) - g(x))]**2) 
    over all x in the specified interval

    """
    return sum([f(X) - g(X)]**2 for X in x)