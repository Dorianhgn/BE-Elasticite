import numpy as np
import matplotlib.pyplot as plt
import scipy.io
import matplotlib.tri as mtri
import matplotlib.tri.triangulation as tri
from mpl_toolkits.mplot3d import Axes3D

class ElasticitySolver:
    def __init__(self, maille=1, E=1000, v=0.3, fvx=100, fvy=0):
        self.maille = maille
        self.E = E
        self.v = v
        self.fv = np.array([fvx, fvy])
        self.parameters = {}
        self.results = {}

    def load_mesh(self):
        if self.maille == 1:
            mat = scipy.io.loadmat('./maillage_1.mat')
        elif self.maille == 2:
            mat = scipy.io.loadmat('./maillage_2.mat')
        elif self.maille == 3:
            mat = scipy.io.loadmat('./maillage_3.mat')
        elif self.maille == 4:
            mat = scipy.io.loadmat('./maillage_trou.mat')
        else:
            raise ValueError('Maille 1,2,3 ou 4')

        self.X = mat['p'].T
        self.C = mat['t'].T[:, :-1] - 1
        self.n_nodes = self.X.shape[0]
        self.n_elems = self.C.shape[0]
        self.n_nodes_elem = 3
        self.ndofe = 2 * self.n_nodes_elem
        self.ndof = 2 * self.n_nodes

        self.parameters.update({
            "n_nodes": self.n_nodes,
            "n_elems": self.n_elems,
            "n_nodes_elem": self.n_nodes_elem,
            "ndofe": self.ndofe,
            "ndof": self.ndof,
        })

    def fun_tri_P1_lag(self, x, y, x_nodes, y_nodes):
        Te = 1 / 2 * np.abs((x_nodes[1] - x_nodes[0]) * (y_nodes[2] - y_nodes[0])
                            - (y_nodes[1] - y_nodes[0]) * (x_nodes[2] - x_nodes[0]))

        N1 = 1 / (2 * Te) * (
            x_nodes[1] * y_nodes[2]
            - x_nodes[2] * y_nodes[1]
            + (y_nodes[1] - y_nodes[2]) * x
            + (x_nodes[2] - x_nodes[1]) * y
        )
        N2 = 1 / (2 * Te) * (
            x_nodes[2] * y_nodes[0]
            - y_nodes[2] * x_nodes[0]
            + (y_nodes[2] - y_nodes[0]) * x
            + (x_nodes[0] - x_nodes[2]) * y
        )
        N3 = 1 / (2 * Te) * (
            x_nodes[0] * y_nodes[1]
            - y_nodes[0] * x_nodes[1]
            + (y_nodes[0] - y_nodes[1]) * x
            + (x_nodes[1] - x_nodes[0]) * y
        )
        N = np.array([N1, N2, N3])
        dNdx = np.array([
            (y_nodes[1] - y_nodes[2]) / (2 * Te),
            (y_nodes[2] - y_nodes[0]) / (2 * Te),
            (y_nodes[0] - y_nodes[1]) / (2 * Te)
        ])

        dNdy = np.array([
            (x_nodes[2] - x_nodes[1]) / (2 * Te),
            (x_nodes[0] - x_nodes[2]) / (2 * Te),
            (x_nodes[1] - x_nodes[0]) / (2 * Te)
        ])
        return [N, dNdx, dNdy]

    def GetBe(self, dNdx, dNdy):
        Be = np.zeros((3, 6))
        Be[0, :3] = dNdx
        Be[1, 3:] = dNdy
        Be[2] = np.hstack((dNdy, dNdx))
        return Be

    def GetNe(self, N):
        Ne_matrix = np.zeros((2, 6))
        Ne_matrix[0, :3] = N
        Ne_matrix[1, 3:] = N
        return Ne_matrix

    def solve(self):
        self.load_mesh()

        # Initialize stiffness matrix and force vector
        K = np.zeros((self.ndof, self.ndof))
        F = np.zeros(self.ndof)

        # Material property matrix
        H = np.array([[self.E / (1 - self.v ** 2), self.v * self.E / (1 - self.v ** 2), 0],
                      [self.v * self.E / (1 - self.v ** 2), self.E / (1 - self.v ** 2), 0],
                      [0, 0, self.E / (2 * (1 + self.v))]])

        # Assemble stiffness matrix and force vector
        for e in range(self.n_elems):
            x_nodes = self.X[self.C[e], 0]
            y_nodes = self.X[self.C[e], 1]

            N, dNdx, dNdy = self.fun_tri_P1_lag(x_nodes[0], y_nodes[0], x_nodes, y_nodes)
            Be = self.GetBe(dNdx, dNdy)
            Te = 1 / 2 * np.abs((x_nodes[1] - x_nodes[0]) * (y_nodes[2] - y_nodes[0])
                                - (y_nodes[1] - y_nodes[0]) * (x_nodes[2] - x_nodes[0]))
            Ke = Te * (Be.T @ H @ Be)

            for i in range(self.n_nodes_elem):
                for j in range(self.n_nodes_elem):
                    K[self.C[e, i], self.C[e, j]] += Ke[i, j]
                    K[self.C[e, i] + self.n_nodes, self.C[e, j]] += Ke[i + self.n_nodes_elem, j]
                    K[self.C[e, i], self.C[e, j] + self.n_nodes] += Ke[i, j + self.n_nodes_elem]
                    K[self.C[e, i] + self.n_nodes, self.C[e, j] + self.n_nodes] += Ke[i + self.n_nodes_elem, j + self.n_nodes_elem]

        for e in range(self.n_elems):
            x_nodes = self.X[self.C[e], 0]
            y_nodes = self.X[self.C[e], 1]

            N1, _, _ = self.fun_tri_P1_lag(x_nodes[0], y_nodes[0], x_nodes, y_nodes)
            N2, _, _ = self.fun_tri_P1_lag(x_nodes[1], y_nodes[1], x_nodes, y_nodes)
            N3, _, _ = self.fun_tri_P1_lag(x_nodes[2], y_nodes[2], x_nodes, y_nodes)
            Ne1 = self.GetNe(N1)
            Ne2 = self.GetNe(N2)
            Ne3 = self.GetNe(N3)

            Te = 1 / 2 * np.abs((x_nodes[1] - x_nodes[0]) * (y_nodes[2] - y_nodes[0])
                                - (y_nodes[1] - y_nodes[0]) * (x_nodes[2] - x_nodes[0]))
            Fe = Te / 3 * (Ne1.T @ self.fv + Ne2.T @ self.fv + Ne3.T @ self.fv)

            for i in range(self.n_nodes_elem):
                F[self.C[e, i]] += Fe[i]
                F[self.C[e, i] + self.n_nodes] += Fe[i + self.n_nodes_elem]

        # Apply Dirichlet boundary conditions
        for n in range(self.n_nodes):
            if self.X[n, 0] == 0:
                K[n, :] = 0
                K[:, n] = 0
                F[n] = 0
                K[n, n] = 1

                K[n + self.n_nodes, :] = 0
                K[:, n + self.n_nodes] = 0
                F[n + self.n_nodes] = 0
                K[n + self.n_nodes, n + self.n_nodes] = 1

        # Solve the linear system
        U = np.linalg.solve(K, F)
        self.results["U"] = U

        # Compute deformed coordinates
        x = self.X[:, 0] + U[:self.n_nodes]
        y = self.X[:, 1] + U[self.n_nodes:]
        self.results["deformed_coordinates"] = (x, y)

        # Compute Von Mises stresses
        T = np.zeros(self.n_nodes)
        SVM = np.zeros(self.n_nodes)

        for e in range(self.n_elems):
            x_nodes = self.X[self.C[e], 0]
            y_nodes = self.X[self.C[e], 1]

            x1, x2, x3 = x_nodes
            y1, y2, y3 = y_nodes

            Te = 0.5 * np.abs((x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1))

            Nm, dNmdx, dNmdy = self.fun_tri_P1_lag(x1, y1, x_nodes, y_nodes)
            Be = self.GetBe(dNmdx, dNmdy)
            Ue = np.r_[U[self.C[e, :]], U[self.C[e, :] + self.n_nodes]]
            sigma = H @ (Be @ Ue)
            sxx, syy, sxy = sigma
            svm = (sxx ** 2 + syy ** 2 + 3 * sxy ** 2 - sxx * syy) ** 0.5

            for i in range(self.n_nodes_elem):
                T[self.C[e, i]] += Te
                SVM[self.C[e, i]] += Te * svm

        SVM /= T
        self.results["SVM"] = SVM

    def plot_results(self):
        x, y = self.results["deformed_coordinates"]
        C = self.C

        # Plot initial and deformed mesh
        Initial_triangulation = mtri.Triangulation(self.X[:, 0], self.X[:, 1], C)
        Deformed_triangulation = mtri.Triangulation(x, y, C)

        plt.figure(1)
        plt.triplot(Initial_triangulation, color='black', label='Initial Mesh')
        plt.triplot(Deformed_triangulation, color='blue', label='Deformed Mesh')
        plt.legend()
        plt.title('Initial and Deformed Mesh')
        plt.show()
        t = mtri.Triangulation(x, y, self.C)

        plt.figure()
        plt.triplot(t)
        plt.tricontourf(t, self.results["SVM"])
        plt.colorbar()
        plt.title('Contraintes de Von Mises')
        plt.show()
