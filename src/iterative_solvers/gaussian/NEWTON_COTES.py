"""Created on Nov 30 16:37:17 2023"""

import enum


class _DataTable:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def print_table(self, precision=100):
        # Find the maximum length of strings in x and y
        n_len = 6
        max_len_x = max(len(str(element)) for element in self.x)
        max_len_y = max(len(str(element)) for element in self.y)

        # Add extra space for padding
        max_len_x += 4
        max_len_y += 4
        n_len += 2

        print('-' * (n_len + max_len_x + max_len_y + 4))
        # Print the table header
        print(f"|{'n':^{n_len}}|{'x':^{max_len_x}}|{'y':^{max_len_y}}|")

        # Print the table border
        print('-' * (n_len + max_len_x + max_len_y + 4))

        # Print the table rows
        for i in range(len(self.x)):
            print(f"|{str(i + 1):^{n_len}}|{str(self.x[i]):^{max_len_x}}|"
                  f"{str(round(self.y[i], precision)):^{max_len_y}}|")

        # Print the table border
        print('-' * (n_len + max_len_x + max_len_y + 4))


@enum.unique
class Iterators(enum.Enum):
    TRAPEZOID = 1
    SIMPSON_13 = 2
    SIMPSON_38 = 3
    BOOLE = 4
    WEDDLE_3H10 = 5.1
    WEDDLE_41H140 = 5.2
