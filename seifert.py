def seifert_matrix(word, strands):
    if set(map(abs, word)) != set(range(1, strands)):
        raise ValueError("The standard Seifert surface is not connected")
    loops = []
    for index in range(1, strands):
        positions = [p for p, value in enumerate(word) if abs(value) == index]
        loops.extend((index, p, q) for p, q in zip(positions, positions[1:]))
    result = [[0] * len(loops) for _ in loops]
    sign = lambda p: 1 if word[p] > 0 else -1
    for row, (i, p, q) in enumerate(loops):
        for col, (j, r, s) in enumerate(loops):
            if row == col:
                result[row][col] = -(sign(p) + sign(q)) // 2
            elif i == j and q == r:
                result[row][col] = (sign(q) + 1) // 2
            elif i == j and s == p:
                result[row][col] = (sign(p) - 1) // 2
            elif j == i + 1 and p < r < q < s:
                result[row][col] = -1
            elif j == i + 1 and r < p < s < q:
                result[row][col] = 1
    return result


def realification(v, slope=4):
    n = len(v)
    s = [[slope * (v[i][j] + v[j][i]) for j in range(n)] for i in range(n)]
    a = [[v[j][i] - v[i][j] for j in range(n)] for i in range(n)]
    return ([s[i] + [-x for x in a[i]] for i in range(n)]
            + [a[i] + s[i] for i in range(n)])
