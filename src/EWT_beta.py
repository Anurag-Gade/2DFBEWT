def EWT_beta(X):

    '''
    Given X, this function helps in the construction of Meyer's wavelet
    '''

    if X<0:
        bm = 0

    elif X>1:
        bm=1

    else:
        bm = (x**4)*(35 - 84*X + 70*(X**2) - 20*(X**3))

    return bm
