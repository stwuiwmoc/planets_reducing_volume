# %%
if __name__ == "__main__":
    # zernike近似作用行列の書き下し
    A_index_list = [
        1, None, 5, 6, 7, None, 11,
        12, None, 16, 17, 18, None, 22,
        23, None, 27, 28, 29, None, 33
    ]

    for i in A_index_list:
        if i is None:
            print(r"\vdots", end="")
            print(r" & & ", end="")
            print(r"\vdots", end="")
            print(r" \\", end="\n")

        else:
            print(r"A _{4, " + str(i) + r"}", end="")
            print(r" & ", end="")
            print(r"\cdots", end="")
            print(r" & ", end="")
            print(r"A _{11, " + str(i) + r"}", end="")
            print(r" \\", end="\n")

    # %%
    a_scalar_index_list = [
        1, None, 5, 6, 7, None, 11,
        13, None, 17, 18, 19, None, 23,
        25, None, 29, 30, 31, None, 35
    ]

    for i in a_scalar_index_list:
        if i is None:
            print(r"\vdots", end="")
            print(r" & & ", end="")
            print(r"\vdots", end="")
            print(r" \\", end="\n")

        elif i in [6, 18, 30]:
            print(r"a _{4, " + str(i) + r"} - a _{4, " + str(i + 6) + r"}", end="")
            print(r" & ", end="")
            print(r"\cdots", end="")
            print(r" & ", end="")
            print(r"a _{11, " + str(i) + r"} - a _{11, " + str(i + 6) + r"}", end="")
            print(r" \\", end="\n")

        else:
            print(r"a _{4, " + str(i) + r"}", end="")
            print(r" & ", end="")
            print(r"\cdots", end="")
            print(r" & ", end="")
            print(r"a _{11, " + str(i) + r"}", end="")
            print(r" \\", end="\n")

    # %%
    a_vector_index_list = [
        1, None, 5, 6, 7, None, 11,
        13, None, 17, 18, 19, None, 23,
        25, None, 29, 30, 31, None, 35
    ]

    for i in a_vector_index_list:
        if i is None:
            print(r":", end="")
            print(r" \\", end="\n")

        elif i in [6, 18, 30]:
            print(r"\overrightarrow{a _{4 \sim 11, " + str(i) + r"}} \, ^T - \overrightarrow{a _{4 \sim 11, " + str(i + 6) + r"}} \, ^T", end="")
            print(r" \\", end="\n")

        else:
            print(r"\overrightarrow{a _{4 \sim 11, " + str(i) + r"}} \, ^T", end="")
            print(r" \\", end="\n")
