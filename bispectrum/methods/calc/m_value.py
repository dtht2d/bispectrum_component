def generate_m_val(j1, j2, j):
    """
    This function generates (m1, m2, m, m1p, m2p, mp) from input set (j1,j2,j)
    and only keep set that satisfy the condition
    """
    # Generate m values
    m = np.linspace(-j, j, int(2 * j + 1))
    mp = m.copy()
    m1 = np.linspace(-j1, j1, int(2 * j1 + 1))
    m1p = m1.copy()
    m2 = np.linspace(-j2, j2, int(2 * j2 + 1))
    m2p = m2.copy()
    s = product(m1, m2, m, m1p, m2p, mp)
    keep_list = []
    for i in s:
        m1, m2, m, m1p, m2p, mp = i[0], i[1], i[2], i[3], i[4], i[5]
        # Check input parameter conditions
        # Condition 1 & 2 & 5
        if not (abs(j1 - j2) <= j <= j1 + j2 and m1 + m2 == m and j1 + j2 - j % 1 != 0.5):
            pass
        if not (abs(j1 - j2) <= j <= j1 + j2 and m1p + m2p == mp and
                j1 + j2 - j % 1 != 0.5):
            pass
        # Condition 3 & 6
        if not all(abs(x) <= y and (x % 0.5 == 0 or x % 1 == 0) for x, y in
                   zip([m1, m2, m], [j1, j2, j])):
            pass
        if not all(abs(x) <= y and (x % 0.5 == 0 or x % 1 == 0) for x, y in
                   zip([m1p, m2p, mp], [j1, j2, j])):
            pass
        # Condition 4
        J = (j1 + j2 + j)
        if J < (int(j1 + j2 + j)) and J < 0:
            raise ValueError("Invalid input parameters: j1, j2, j \
                              must not exceed a positive integer J")
        # Condition 7
        if not all(isinstance(x, (int, float, np.integer, np.floating)) and x >= 0
                   and (x % 0.5 == 0 or x % 1 == 0) for x in [j1, j2, j]):
            raise ValueError("Invalid input parameters: j1, j2, j must be integer \
                             or half-integer non-negative numbers")
        # Condition 8
        if not all(isinstance(x, (int, float, np.integer, np.floating)) and x >= 0
                   and x % 1 == 0 for x in [j1 + m1, j2 + m2, j + m, j1 + j2 + j]):
            pass
        if not all(isinstance(x, (int, float, np.integer, np.floating)) and x >= 0
                   and x % 1 == 0 for x in [j1 + m1p, j2 + m2p, j + mp, j1 + j2 + j]):
            pass
        else:
            keep_list.append(i)
    return keep_list