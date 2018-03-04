from scipy.stats import multivariate_normal
import operator
import numpy as np

def correlated_significance(v1_i,v2_i, v1_var, v2_var, **kwargs):
    """
    :param v1_i: mean of intensity of the first
    :param v2_i: mean of intensity of the second
    :param v1_var: variance of the first
    :param v2_var: variance of the second
    :param where: function default: zip(i1,i2) ; i1= 12 = np.linspace(imin, imax,200)
    :return: (where on measured), signficance_value
    """
    cov = [[v1_var, 0], [0, v2_var]]
    mean = [v1_i, v2_i]

    #default sampling line
    imin = np.min([v1_i,v2_i])
    imax = np.max([v1_i,v2_i])
    distance = np.max([v1_i,v2_i]) - np.min([v1_i,v2_i])

    line = np.linspace(imin-distance, imax+distance, 500, endpoint=False)
    where=kwargs.get("where",zip(line, line))

    if np.isnan(cov).any():
        return (np.nan,np.nan), np.nan

    significance = multivariate_normal.pdf(where, mean=mean, cov=cov)

    #maxima
    where_max, significance_max = max(zip(where,significance), key=operator.itemgetter(1))



    sig_int = ((where_max[0]-v1_i)**2)/(np.sqrt(2)*v1_var) + ((where_max[1]-v2_i)**2)/(np.sqrt(2)*v2_var)
    #integral
    #sig_int = significance.sum()#*(np.linalg.norm(np.array(where[1])-np.array(where[0])))


    return where_max, sig_int



