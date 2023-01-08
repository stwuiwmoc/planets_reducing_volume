# %%
if __name__ == "__main__":
    # cell 1
    A_index_list = [
        1, None, 5, 6, 7, None, 11,
        12, None, 16, 17, 18, None, 22,
        23, None, 27, 28, 29, None, 33
    ]

    for i in A_index_list:
        if i is None:
            print(r"    ", end="")
            print(r"\vdots", end="")
            print(r" & & ", end="")
            print(r"\vdots", end="")
            print(r" \\", end="\n")

        else:
            print(r"    ", end="")
            print(r"A _{1, " + str(i) + r"}", end="")
            print(r" & ", end="")
            print(r"\cdots", end="")
            print(r" & ", end="")
            print(r"A _{11, " + str(i) + r"}", end="")
            print(r" \\", end="\n")

    # %%
    # cell 2
    a_scalar_index_list = [
        1, None, 5, 6, 7, None, 11,
        13, None, 17, 18, 19, None, 23,
        25, None, 29, 30, 31, None, 35
    ]

    for i in a_scalar_index_list:
        if i is None:
            print(r"    ", end="")
            print(r"\vdots", end="")
            print(r" & & ", end="")
            print(r"\vdots", end="")
            print(r" \\", end="\n")

        elif i in [6, 18, 30]:
            print(r"    ", end="")
            print(r"a _{1, " + str(i) + r"} - a _{1, " + str(i + 6) + r"}", end="")
            print(r" & ", end="")
            print(r"\cdots", end="")
            print(r" & ", end="")
            print(r"a _{11, " + str(i) + r"} - a _{11, " + str(i + 6) + r"}", end="")
            print(r" \\", end="\n")

        else:
            print(r"    ", end="")
            print(r"a _{1, " + str(i) + r"}", end="")
            print(r" & ", end="")
            print(r"\cdots", end="")
            print(r" & ", end="")
            print(r"a _{11, " + str(i) + r"}", end="")
            print(r" \\", end="\n")

    # %%
    # cell 3
    a_vector_index_list = [
        1, None, 5, 6, 7, None, 11,
        13, None, 17, 18, 19, None, 23,
        25, None, 29, 30, 31, None, 35
    ]

    for i in a_vector_index_list:
        if i is None:
            print(r"    ", end="")
            print(r":", end="")
            print(r" \\", end="\n")

        elif i in [6, 18, 30]:
            print(r"    ", end="")
            print(r"\vec{a _{" + str(i) + r"}} \, ^T - \vec{a _{" + str(i + 6) + r"}} \, ^T", end="")
            print(r" \\", end="\n")

        else:
            print(r"    ", end="")
            print(r"\vec{a _{" + str(i) + r"}} \, ^T", end="")
            print(r" \\", end="\n")

    # %%
    # cell 4
    x_index_list = [
        1, None, 5, 6, 7, None, 11,
        12, None, 16, 17, 18, None, 22,
        23, None, 27, 28, 29, None, 33
    ]

    for i in x_index_list:
        if i is None:
            print(r"    ", end="")
            print(r"\vdots", end="")
            print(r" \\", end="\n")

        else:
            print(r"    ", end="")
            print(r"x _{" + str(i) + r"}", end="")
            print(r" \\", end="\n")

    # %%
    # cell 5
    y_index_list = [
        1, None, 5, 6, 7, None, 11,
        13, None, 17, 18, 19, None, 23,
        25, None, 29, 30, 31, None, 35
    ]

    for i in y_index_list:
        if i is None:
            print(r"    ", end="")
            print(r"\vdots", end="")
            print(r" \\", end="\n")

        elif i in [6, 18, 30]:
            print(r"    ", end="")
            print(r"(y _{" + str(i) + r"} + y _{" + str(i + 6) + r"}) / 2", end="")
            print(r" \\", end="\n")

        else:
            print(r"    ", end="")
            print(r"y _{" + str(i) + r"}", end="")
            print(r" \\", end="\n")

    # %%
    # cell 6
    x_index_list = [
        1, None, 5, 6, 7, None, 11,
        12, None, 16, 17, 18, None, 22,
        23, None, 27, 28, 29, None, 33
    ]
    y_index_list = [
        1, None, 5, 6, 7, None, 11,
        13, None, 17, 18, 19, None, 23,
        25, None, 29, 30, 31, None, 35
    ]

    for i in x_index_list:
        if i is None:
            print(r"| $\cdots$ ", end="")

        else:
            print(r"| $x _{" + str(i) + r"}$ ", end="")
    print("|")

    for i in range(len(x_index_list)):
        print(r"| - ", end="")
    print("|")

    for i in y_index_list:
        if i is None:
            print(r"| $\cdots$ ", end="")

        elif i in [6, 18, 30]:
            print(r"| $(y _{" + str(i) + r"} + y _{" + str(i + 6) + r"}) / 2$ ", end="")

        else:
            print(r"| $y _{" + str(i) + r"}$ ", end="")
    print("|")

    # %%
    # cell 7
    for i in range(11):
        num = i + 1
        print(r"    ", end="")
        print(r"x _{" + str(num) + r"}", end="")
        print(r" \\", end="\n")

    print(r"    \vdots \\", end="\n")

    print()

    for i in range(11):
        num = i + 1
        if num == 6:
            print(r"    ", end="")
            print(r"(y _{" + str(num) + r"} + y _{" + str(num + 6) + r"}) / 2", end="")
            print(r" \\", end="\n")

        else:
            print(r"    ", end="")
            print(r"y _{" + str(num) + r"}", end="")
            print(r" \\", end="\n")

    print(r"    \vdots \\", end="\n")

    # %%
    # cell 8
    for i in range(11):
        num = i + 1
        print(r"    ", end="")
        print(r"\vec{\vphantom{d} a _{" + str(num) + r"}}", end="")
        print(r" \\", end="\n")

    print(r"    \vdots \\", end="\n")

    print()

    for i in range(11):
        num = i + 1
        if num == 6:
            print(r"    ", end="")
            print(r"\mathcal{Z}(\vec{d _{" + str(num) + r"}} - \vec{d _{" + str(num + 6) + r"}})", end="")
            print(r" \\", end="\n")

        else:
            print(r"    ", end="")
            print(r"\mathcal{Z}(\vec{d _{" + str(num) + r"}})", end="")
            print(r" \\", end="\n")

    print(r"    \vdots \\", end="\n")

    # %%
    # cell 9
    a_index_list = [
        1, None, 5, 6, 7, None, 11,
        12, None, 16, 17, 18, None, 22,
        23, None, 27, 28, 29, None, 33
    ]
    d_index_list = [
        1, None, 5, 6, 7, None, 11,
        13, None, 17, 18, 19, None, 23,
        25, None, 29, 30, 31, None, 35
    ]

    for i in a_index_list:
        if i is None:
            print(r"| $\cdots$ ", end="")

        else:
            print(r"| $\vec{a _{" + str(i) + r"}}$ ", end="")
    print("|")

    for i in range(len(x_index_list)):
        print(r"| - ", end="")
    print("|")

    for i in d_index_list:
        if i is None:
            print(r"| $\cdots$ ", end="")

        elif i in [6, 18, 30]:
            print(r"| $\mathcal{Z}( \vec{d _{" + str(i) + r"}} - \vec{d _{" + str(i + 6) + r"}})$ ", end="")

        else:
            print(r"| $\mathcal{Z}( \vec{d _{" + str(i) + r"}})$ ", end="")
    print("|")

    # %%
    # cell 10
    print(r"||||")
    print(r"| - | - | - |")

    for i in range(11):
        for j in range(3):
            x_num = (i + 1) + 11 * j
            y_num = (i + 1) + 12 * j

            if y_num in [6, 18, 30]:
                print(r"| ", end="")
                print(r"$x _{" + str(x_num) + r"}$", end="")
                print(r" = ", end="")
                print(r"$(y _{" + str(y_num) + r"} + y _{" + str(y_num + 6) + r"}) / 2$ ", end="")
            else:
                print(r"| ", end="")
                print(r"$x _{" + str(x_num) + r"}$", end="")
                print(r" = ", end="")
                print(r"$y _{" + str(y_num) + r"}$ ", end="")
        print("|")

    # %%
    # cell 11
    print(r"||||")
    print(r"| - | - | - |")

    for i in range(11):
        for j in range(3):
            a_num = (i + 1) + 11 * j
            d_num = (i + 1) + 12 * j

            if d_num in [6, 18, 30]:
                print(r"| ", end="")
                print(r"$\vec{a _{" + str(a_num) + r"}}$", end="")
                print(r" = ", end="")
                print(r"$\mathcal{Z}(\vec{d _{" + str(d_num) + r"}} - \vec{d _{" + str(d_num + 6) + r"}})$ ", end="")
            else:
                print(r"| ", end="")
                print(r"$\vec{a _{" + str(a_num) + r"}}$", end="")
                print(r" = ", end="")
                print(r"$\mathcal{Z}(\vec{d _{" + str(d_num) + r"}})$ ", end="")
        print("|")
