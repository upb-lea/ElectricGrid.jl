import numpy as np
import matplotlib.pyplot as plt
import copy
import random
import networkx as nx

class NodeConstructor:
    def __init__(self, num_source, num_loads, parameter, S2S_p=0.1, S2L_p=0.8, CM=None):
        self.num_source = num_source
        self.num_loads = num_loads
        self.tot_ele = num_source + num_loads
        self.S2S_p = S2S_p
        self.S2L_p = S2L_p
        self.cntr = 0
        self.num_connections = 0

        # unpack parameters
        self.parameter = parameter
        self.R_source = parameter['R_source']
        self.L_source = parameter['L_source']
        self.C_source = parameter['C_source']
        self.R_cabel = parameter['R_cabel']
        self.L_cabel = parameter['L_cabel']
        self.R_load = parameter['R_load']

        if isinstance(CM, np.ndarray):
            assert CM.shape[0] == self.tot_ele, "Expect CM to have the same number of elements as tot_ele."
            self.CM = CM
            self.num_connections = np.amax(CM)
        elif CM == None:
            self.generate_CM()
        else:
            raise f"Expect CM to be an np.ndarray or None not {type(CM)}."

    def tobe_or_n2b(self, x, p):

        # To count up the connection, cntr is returned.
        # If only one type of cabel is used this is not necessary an can be replaced by 1

        if x < p:
            self.cntr += 1
            return self.cntr
        else:
            x = 0
            return x

    def count_up(self):
        self.cntr += 1
        return self.cntr

    def generate_CM(self):
        '''
        Constructs the CM
        '''

        self.cntr = 0

        # Get a upper triangular matrix
        mask = np.tri(self.tot_ele).T
        CM = np.random.rand(self.tot_ele, self.tot_ele) * mask  # fill matrix with random entries between [0,1]
        CM = CM - np.eye(CM.shape[0]) * np.diag(CM)  # delet diagonal bc no connection with itself

        # Go throught the matrix
        # -1 bc last entrie is 0 anyway
        # start at i, bc we need to check only upper triangle

        for i in range(self.tot_ele - 1):
            for j in range(i, self.tot_ele - 1):
                if j >= self.num_source - 1:  # select propability according to column
                    CM[i, j + 1] = self.tobe_or_n2b(CM[i, j + 1], self.S2L_p)
                else:
                    CM[i, j + 1] = self.tobe_or_n2b(CM[i, j + 1], self.S2S_p)

        # Check if CM is valid
        # First: row -> col

        for i in range(self.tot_ele - 1):
            row = CM[i, i + 1:]  # get upper triangle row

            if np.count_nonzero(row) < 2:  # check if at least 2 entries are present in row
                values_to_set = (row == 0).sum().clip(0, 2)  # safe how many connections we could set
                idx = np.where(0 == row)  # save indicies of 0 entries
                idx_ = idx[0].tolist()
                idx_corr = copy.copy(idx_)  # copy the list bc python

                for k, ix in enumerate(idx_):  # enumerate colums and check if a row contains at least 2 entries
                    if np.count_nonzero(CM[:i + 1 + ix,
                                        i + 1 + ix]) >= 2:  # remove all coloums from the list which have more then 2 entries
                        idx_corr.remove(ix)

                if len(idx_corr) == 0:  # check if entries a left
                    continue
                else:
                    samples = values_to_set.clip(0,
                                                 len(idx_corr))  # get possible entries // can't be greater then idx list
                    idx_rnd = random.sample(range(0, len(idx_corr)),
                                            samples)  # choose a random index from the remaining entries
                    idx_rnd = np.array(idx_rnd)
                    col_idx = idx_rnd + i + 1  # col idx in the CM-Matrix
                    for _, col_ix in enumerate(col_idx):
                        CM[i, col_ix] = self.count_up()  # set entry and inc counter

        # and the other way: col -> row

        for i in range(0, self.tot_ele - 1):
            col = CM[:i + 1, i + 1]  # get upper triangle col

            if np.count_nonzero(col) < 2:  # check if a connection with at least 2 entries are present
                values_to_set = (col == 0).sum().clip(0, 2)  # safe how many connections we could set
                idx = np.where(0 == col)  # save indicies of 0 entries
                idx_ = idx[0].tolist()
                idx_corr = copy.copy(idx_)

                for k, ix in enumerate(idx_):  # enumerate colums and check if a row contains at least 2 entries
                    row_idx = ix
                    col_start = ix + 1
                    if np.count_nonzero(CM[row_idx,
                                        col_start:]) >= 2:  # remove all rows from the list which have more then 2 entries

                        idx_corr.remove(ix)

                if len(idx_corr) == 0:  # check if entries a left
                    continue
                else:
                    samples = values_to_set.clip(0,
                                                 len(idx_corr))  # get possible entries // can't be greater then idx list
                    idx_rnd = random.sample(range(0, len(idx_corr)),
                                            samples)  # choose a random index from the remaining 0 entries
                    idx_rnd = np.array(idx_rnd)
                    row_idx = idx_rnd + 1
                    for _, row_ix in enumerate(row_idx):
                        CM[row_ix, i + 1] = self.count_up()  # set entry

        CM = CM - CM.T  # copy with negative sign to lower triangle

        self.CM = CM

        self.num_connections = self.cntr
        pass

    def get_A_source(self):
        # this matrix is always a 2x2 for inverter
        A_source = np.zeros((2, 2))
        A_source[0, 0] = -self.R_source / self.L_source
        A_source[0, 1] = -1 / self.L_source
        A_source[1, 0] = 1 / self.C_source
        return A_source

    def get_B_source(self):
        B_source = np.zeros((2, 1))
        B_source[0, 0] = 1 / self.L_source
        return B_source

    def get_A_col(self, source_x):
        # this matrix is (2 x num_connections)
        # for this case self.C_source is assumed to be just an int.
        # Later self.C_source could be an array with the diffrent paramters and would be indexed via self.C_source[source_x]

        A_col = np.zeros((2, self.num_connections))

        CM_row = self.CM[source_x - 1]

        indizes = list(CM_row[CM_row != 0])  # get entries unequal 0
        signs = np.sign(indizes)  # get signs
        indizes_ = indizes * signs  # delet signs from indices
        indizes_.astype(dtype=np.int32)

        for i, (idx, sign) in enumerate(zip(indizes_, signs)):
            idx = int(idx)

            A_col[1, idx - 1] = sign * -1 / self.C_source

        return A_col

    def get_A_row(self, source_x):

        A_row = np.zeros((2, self.num_connections))

        CM_col = self.CM[source_x - 1]

        indizes = list(CM_col[CM_col != 0])  # get entries unequal 0

        signs = np.sign(indizes)  # get signs
        indizes_ = indizes * signs  # delet signs from indices

        for i, (idx, sign) in enumerate(zip(indizes_, signs)):
            idx = int(idx)
            A_row[1, idx - 1] = sign * 1 / self.L_cabel

        return A_row.T

    def get_A_transitions(self):
        A_transitions = np.zeros((self.num_connections, self.num_connections))
        for i in range(1, self.num_connections + 1):
            (row, col) = np.where(self.CM == i)
            (row_idx, col_idx) = (row[0], col[0])

            # check if its a S2S connection
            if col_idx < self.num_source:  # row_idx < self.num_source and

                A_transitions[i - 1, i - 1] = -self.R_cabel / self.L_cabel  # self.R_cabel[i] and self.L_cabel[i]

            # Then it has to be S2L
            else:
                # easy diagonal entry
                A_transitions[i - 1, i - 1] = -(
                            self.R_cabel + self.R_load) / self.L_cabel  # (self.R_cabel[i] + self.R_load[col_idx])/self.L_cabel[i] -> self.R_load[col_idx]? not sure

                # search for other connections to this specific load in the colum
                CM_col = self.CM[:, col_idx]

                mask = np.logical_and(CM_col > 0, CM_col != i)  # i bc we already cover this case
                indizes = list(CM_col[mask])

                # cross entries for the other connections to this load
                for j, idx in enumerate(indizes):
                    idx = int(idx)
                    A_transitions[
                        i - 1, idx - 1] = -self.R_load / self.L_cabel  # self.L_cabel[i] if LT is an arry with diffrent values and self.R_load[col_idx]?

        return A_transitions

    def generate_A(self):
        # get A_source
        A_source = np.zeros((2 * self.num_source, 2 * self.num_source))  # construct matrix of zeros
        A_source_list = [self.get_A_source() for i in range(self.num_source)]

        for i, ele in enumerate(A_source_list):
            start = 2 * i
            stop = 2 * i + 2
            A_source[start:stop, start:stop] = ele

        # get A_col
        A_col = np.zeros((2 * self.num_source, self.num_connections))
        A_col_list = [self.get_A_col(i) for i in range(1, self.num_source + 1)]  # start at 1 bc Source 1 ...

        for i, ele in enumerate(A_col_list):
            start = 2 * i
            stop = 2 * i + 2
            A_col[start:stop, :] = ele

        # get A_row
        A_row = np.zeros((self.num_connections, 2 * self.num_source))
        A_row_list = [self.get_A_row(i) for i in range(1, self.num_source + 1)]  # start at 1 bc Source 1 ...

        for i, ele in enumerate(A_row_list):
            start = 2 * i
            stop = 2 * i + 2
            A_row[:, start:stop] = ele

        A_transitions = self.get_A_transitions()

        A = np.block([[A_source, A_col],
                      [A_row, A_transitions]])

        return A

    def generate_B(self):

        B = np.zeros((2 * self.num_source + self.num_connections, self.num_source))

        B_source_list = [self.get_B_source() for i in range(1, self.num_source + 1)]  # start at 1 bc Source 1 ...
        for i, ele in enumerate(B_source_list):
            #             start_c = i
            #             stop_c = i+1
            start_r = 2 * i
            stop_r = 2 * i + 2
            B[start_r:stop_r, i:i + 1] = ele
        return B

    def generate_C(self):
        return np.eye(2 * self.num_source + self.num_connections)

    def generate_D(self):
        return 0

    def get_sys(self):
        A = self.generate_A()
        B = self.generate_B()
        C = self.generate_C()
        D = self.generate_D()
        return (A, B, C, D)

    def draw_graph(self):

        edges = []
        for i in range(1, self.num_connections + 1):
            (row, col) = np.where(self.CM == i)
            (row_idx, col_idx) = (row[0] + 1, col[0] + 1)
            edges.append((row_idx, col_idx))

        G = nx.Graph(edges)
        nx.draw(G, with_labels=True)
        plt.show()

        pass