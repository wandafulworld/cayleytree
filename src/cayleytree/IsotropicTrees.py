"""
Cayley tree inspired models that are isotropic, i.e. edges are undirected.
This allows them to be solved using Mahan's approach.

Functions:
1. For each model, we can construct its networkx graph object and adjacency matrix and draw the model using the networkx package
2. Each model can be evaluated for its adjacency spectrum using exact diagonalization or Mahans approach.

Note that in a tight-binding (TB) model, the adjacency spectrum corresponds to the spectrum of the TB-Hamiltonian.
"""
from logging import raiseExceptions

from AbstractTree import IsotropicAbstractTree
from AbstractTree import AbstractTree
import numpy as np
import scipy as sp
import networkx as nx
import matplotlib.pyplot as plt
import math
import logging
import itertools
from numpy.polynomial.polynomial import Polynomial

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)



class CayleyTree(IsotropicAbstractTree):
    def __init__(self,M,k,force_graph_object_creation=False, next_nearest_neighbors=False, t2=1):
        """
        A simple Cayley tree with M shells and k children per node. If N > 4000, the graph object will not be instantiated and certain
        methods will not be available. To change this, you need to set force_graph_object_creation to true.
        :param M: int, Number of shells of the Cayley Tree
        :param k: int, Number of children per node
        :param force_graph_object_creation: bool, if True the object will instantiate a graph-object even if there are more than 4000 nodes.
        """
        self.k = k
        self.M = M
        self.N = int(1 + (k+1)*(self.k**M -1)/(self.k-1)) # Number of nodes on the cayley tree
        self.forced_graph = force_graph_object_creation
        self.next_nearest_neighbors = next_nearest_neighbors
        self.t2 = t2 # for tuning of NNN-hopping

        print(self.N)

        if self.N < 4000 or force_graph_object_creation:
            self.G = IsotropicAbstractTree._tree_creator(self.N,self.k,CayleyTree._tree_edges) # networkx graph object

            if next_nearest_neighbors:
                # Now we construct the next-nearest neighbour hopping
                # First construct the lower-shell nn hoping
                for node in self.G.nodes:
                    nextnearestneighbours = AbstractTree.sub_nnneighbours(self.G, node)
                    for nnnode in nextnearestneighbours:
                        self.G.add_edge(node, nnnode, weight=self.t2)

                for node in self.G.nodes:
                    if (node + (k-2)) % k == 0 and node > k+1:
                        for combination in itertools.combinations(range(node, node+k), 2):
                            self.G.add_edge(combination[0], combination[1], weight=self.t2)
                    else:
                        continue

                # Manual add the inner circle
                for combination in itertools.combinations(range(1,k+2),2):
                    self.G.add_edge(combination[0], combination[1], weight=self.t2)

            self._A = nx.adjacency_matrix(self.G) # Sparse Matrix

    @staticmethod
    def _tree_edges(n, r):
        if n == 0:
            return
        # helper function for trees
        # yields edges in rooted tree at 0 with n nodes and branching ratio r
        nodes = iter(range(n))
        parents = [next(nodes)]  # stack of max length r

        first_time = True
        r = r + 1

        while parents:
            source = parents.pop(0)
            for i in range(r):
                try:
                    target = next(nodes)
                    parents.append(target)
                    yield source, target, {'weight': 1}
                except StopIteration:
                    break
            if first_time:
                r = r - 1
            first_time = False

    @property
    def A(self):
        return self._A.todense()

    def polynomial_solver(self):
        '''
        Constructs the eigenvalue polynomials according to Ostilli (2024) and solves for their roots.
        Only finds the geometric degeneracies.
        :return: eigenvalues of the cayley tree
        '''
        domain = [-2*np.sqrt(self.k),2*np.sqrt(self.k)]
        poly = [Polynomial([1],domain=domain,window=domain),Polynomial([0,1],domain=domain,window=domain)]
        x = Polynomial([0,1],domain=domain,window=domain) # For multiplication later on

        eigenval = poly[1].roots().tolist()  # List where we store the eigenvalues

        for l in range(2,self.M): # ToDo: If M < 3 this will not be triggered -> need case response for M = 1
            poly.append(x*poly[l-1] - (self.k)*poly[l-2])
            eigenval.extend(poly[l].roots().tolist())

        # Last polynomial P_L that isn't produced by the above
        poly.append(x*poly[self.M-1] - (self.k)*poly[self.M-2])
        eigenval.extend(poly[self.M].roots().tolist())
        asymmetric_polynomial = x*poly[self.M] - (self.k + 1)*poly[self.M - 1]
        eigenval.extend(asymmetric_polynomial.roots().tolist())

        return eigenval

    def _peak_finder(self,polynomial,n):
        x, y = polynomial.linspace(n)
        peaks, _ =sp.signal.find_peaks(np.abs(-y)) # Invert b/c we want minimas
        return x[peaks].tolist()

    def polynomial_minima_searcher(self,n): #ToDO: Fix to find correct eigenvalues
        domain = [-2*np.sqrt(self.k),2*np.sqrt(self.k)]
        poly = [np.polynomial.polynomial.Polynomial([1],domain=domain,window=domain),np.polynomial.polynomial.Polynomial([0,1],domain=domain,window=domain)]
        x = np.polynomial.polynomial.Polynomial([0,1],domain=domain,window=domain) # For multiplication later on

        eigenvals = []
        eigenvals.extend(self._peak_finder(poly[1],n))


        for l in range(2,self.M): # ToDo: If M < 3 this will not be triggered -> need case response for M = 1
            poly.append(x*poly[l-1] - (self.k)*poly[l-2])
            eigenvals.extend(self._peak_finder(poly[l], n))

        # Last polynomial P_L that isn't produced by the above
        poly.append(x*poly[self.M-1] - (self.k)*poly[self.M-2])
        eigenvals.extend(self._peak_finder(poly[self.M], n))
        asymmetric_polynomial = x * poly[self.M] - (self.k + 1)*poly[self.M - 1]
        eigenvals.extend(self._peak_finder(asymmetric_polynomial, n))

        return eigenvals

    def eff_hamiltonian_constructor(self,d,l=None):
        """
        Constructs a one-dimensional Hamiltonian with the off-diagonal elements being sqrt(r)
        :param d: int, Dimension of the H Matrix
        :return: h: (d+1,d+1)-dim np-array, matrix of a 1-D chain hamiltonian of dimension d
        """

        if self.next_nearest_neighbors:
            if l == 0: #shell symmetric
                h = np.eye(N=d)*(self.k-1)*self.t2 + np.eye(N=d,k=1)*np.sqrt(self.k) + np.eye(N=d,k=-1)*np.sqrt(self.k) + np.eye(N=d,k=2)*self.k*self.t2 + np.eye(N=d,k=-2)*self.k*self.t2
                # Boundary Terms
                h[0][0] = 0
                h[1][1] = self.k*self.t2
                h[1][0] = np.sqrt(self.k + 1)
                h[0][1] = np.sqrt(self.k + 1)
                h[2][0] = np.sqrt((self.k+1)*self.k)*self.t2
                h[0][2] = np.sqrt((self.k + 1) * self.k)*self.t2

            else: # Non-Symmetric
                h = np.eye(N=d)*(self.k-1)*self.t2 + np.eye(N=d,k=1)*np.sqrt(self.k) + np.eye(N=d,k=-1)*np.sqrt(self.k) + np.eye(N=d,k=2)*self.k*self.t2 + np.eye(N=d,k=-2)*self.k*self.t2
                # Boundary Terms
                h[0][0] = -1*self.t2

        else:
            h = np.eye(N=d, k=1) * np.sqrt(self.k) + np.eye(N=d, k=-1) * np.sqrt(self.k)

        return h


    def _eff_hamiltonian_list(self):

        Hs = []
        degeneracies = []

        if self.next_nearest_neighbors:
            Hs.append(self.eff_hamiltonian_constructor(self.M+1,0)) # Shell symmetric
            degeneracies.append(1)

            Hs.append(self.eff_hamiltonian_constructor(self.M,1))
            degeneracies.append(self.k) # degeneracy of k according to kth-root degeneracy

            for l in range(2,self.M + 1):
                Hs.append(self.eff_hamiltonian_constructor(self.M + 1 -l,l))
                degeneracies.append((self.k-1)*(self.k+1)*self.k**(l-2))


        else:
            # symmetric psi_0 != 0
            h = self.eff_hamiltonian_constructor(self.M+1)
            h[0][1] = np.sqrt(self.k+1)
            h[1][0] = np.sqrt(self.k+1)
            Hs.append(h)
            degeneracies.append(1)

            # symmetric psi_0 = 0
            # Consists of (M-1) states per branch r and an additional state 0. No transitions between the
            # different matrix blocks
            Hs.append(self.eff_hamiltonian_constructor(self.M))
            degeneracies.append(self.k) # degeneracy of k

            #antisymmetric states with degeneracy of each state in shell l being (K+1)*K^(l-1) degenerate
            for l in range(2,self.M + 1):
                Hs.append(self.eff_hamiltonian_constructor(self.M + 1 -l))
                degeneracies.append((self.k - 1)*(self.k+1)*self.k**(l-2))


        return Hs, degeneracies

    def shell_list(self):
        """Returns a list of lists nodes, with each list containing the nodes of the same shell
            Makes use of the fact that the construction algorithm constructs each shell after the other

            Returns
            -------
            shell_lists : Multiple lists corresponding to the number of Shells M
                each list containing the nodes of that shell
            """
        nodes = list(self.G.nodes)
        shell_lists = [[0]]  # already contains the 0th node
        l = 1  # position of last added node
        for shell_number in range(1, self.M + 1):
            n = (self.k + 1) * (self.k ** (shell_number - 1))
            shell_lists.append(nodes[l:(l + n)])
            l += n

        return shell_lists

    def _asymmetric_states(self,node_alpha,return_eigenstates = False):
        """
        For a given node alpha, constructs all asymmetric states orginating from the node in position basis.
        When supplied a list of eigenvectors of the effective hamiltonian connecting these states, it provides
        the eigenstates in position basis.
        :param node_alpha: The node from which the asymmetric states branch of
        :param return_eigenstates: If True, will return the eigenstates of that block of asymmetric states in position basis
        :return: N-dim array of asymmetric states or eigenstates of these asymmetric states in Position basis
        """
        subtree = self._subbrancher(node_alpha,flatten=True)
        number_of_shells = np.round(len(subtree)/(self.k) - 1).astype(int)

        branches = [[] for _ in range(self.k)]

        # From the Subtree, construct the subbranches with each subbranch consisting of lists for each shell
        i = 0 # position of last node added
        for l in range(number_of_shells):
            n = self.k**l # number of nodes per branch per shell
            for branch_number in range(self.k):
                branches[branch_number].append(subtree[i:i+n])
                i += n

        # Construct AntiSymmetric States
        asymmetric_states = []
        nth_roots = [self.nth_root_of_unity(self.k,i) for i in range(self.k)]
        nodes = [[] for _ in range(self.k)]

        for shell in range(number_of_shells):
            for branch in range(self.k):
                nodes[branch] = [item for row in branches[branch][:shell+1] for item in row] # Append a flattened array of same branch
                #print('Node of ', shell, ' is ', nodes)
            asymmetric_states.append(self._linear_combination_vector(self.N,nodes,nth_roots))

        if return_eigenstates:
            eigenvalue, eigenvectors = sp.linalg.eig(self.eff_hamiltonian_constructor(number_of_shells))
            mat = np.stack(asymmetric_states,axis=1)
            eigenstates = []
            for state in eigenvectors:
                eigenstates.append(mat.dot(state))

            return eigenstates


        return asymmetric_states




class LiebCayley(IsotropicAbstractTree):
    def __init__(self,M,k,force_graph_object_creation=False):
        """
        Constructs a Lieb-Cayley tree. This is a Cayley tree with an additional node placed in the middle of each
        edge of the tree. This model has been shown to host in-gap states of topological origin. Note that the behavior
        differs depending on whether M is chosen to be odd or even. If N > 4000, the graph object will not be instantiated and certain
        methods will not be available. To change this, you need to set force_graph_object_creation to true.
        :param M: int, the number of shells of your tree.
        :param k: int, number of children per node
        :param save_ram:
        """
        self.M = M
        self.mc = math.floor(M/2) # Number of Cayley Shells
        self.ml = math.ceil(M/2) # Number of Lieb Shells
        self.k = k # Connectivity (Degree of each node = k + 1)
        self.forced_graph = force_graph_object_creation
        logger.info('Initiating LiebCayley Tree')
        self.N = int(1 + ((k + 1)/(k-1))*(k**self.mc + k**self.ml -2)) # Number of Nodes

        if self.N < 4000 or force_graph_object_creation:
            self.G = IsotropicAbstractTree._tree_creator(self.N,self.k,LiebCayley._tree_edges) # networkx graph object
            self._A = nx.adjacency_matrix(self.G,weight='weight') # Sparse Matrix


    @staticmethod
    def _tree_edges(n,r):
        """
        Iteratively defines the tree structure
        :param n: Number of nodes
        :param r: Connectivity (k)
        :return: A list of tuples that define the edges of our tree
        """
        if n == 0:
            return
        # helper function for trees
        # yields edges in rooted tree at 0 with n nodes and branching ratio r
        nodes = iter(range(n))
        cayley_shell = [next(nodes)]  # stack of max length r
        lieb_shell = []
        r = r + 1
        first_run = True
        # Iterative filling of shells
        while cayley_shell:
            source = cayley_shell.pop(0)
            for i in range(r):
                try:
                    target = next(nodes)
                    lieb_shell.append(target)
                    yield source, target
                except StopIteration:
                    break
            if first_run:
                r = r - 1
                first_run = False

            if not cayley_shell:
                for lieb_node in lieb_shell:
                    try:
                        target = next(nodes)
                        cayley_shell.append(target)
                        yield lieb_node,target
                    except StopIteration:
                        break
                lieb_shell = []

    @property
    def A(self):
        return self._A.todense()


    def shell_list(self):
        """Returns a list of lists nodes, with each list containing the nodes of the same shell
        Makes use of the fact that the construction algorithm constructs each shell after the other

        Returns
        -------
        shell_lists : Multiple lists corresponding to the number of Shells M + 1 (incl. the center)
            each list containing the nodes of that shell
        """
        nodes = list(self.G.nodes)
        shell_lists = [[0]]
        l = 1
        for shell_number in range(1, self.M + 1):
            n = (self.k + 1) * (self.k ** (math.ceil(shell_number/2) - 1))
            shell_lists.append(nodes[l:(l + n)])
            l += n
        return shell_lists




    def _eff_hamiltonian_constructor(self,l,J = None):
        """
        Constructs the effective Hamiltonian of the Lieb-Cayley tree for a given shell number.
        Automatically adapts to M even or M odd cases. For M even, the matrix will have dimension M - l + 1.
        :param l: Number of the shell we're starting our construction from. For l = 0 we start from |0>
        :return: 2D-nparray of effective hamiltonian h and scalar of degeneracy d of the hamiltonian
        """
        if not J and J != 0:
            J = 1
        a = [np.sqrt(self.k)*J,1]
        offdiag = np.tile(a,reps=int(np.ceil(self.M/2)))
        offdiag[0] = np.sqrt(self.k + 1)*J
        # print('Length Offdiagonal:', len(offdiag))

        h = np.diag(offdiag[l:self.M],k=1) + np.diag(offdiag[l:self.M],k=-1) #The M automatically cuts off the last element of the offdiagonal if M is odd
        # print('Shape of H: ',np.shape(h))
        # print(l,',h:',h,self.M%2)
        # print('-----------------------------------------------------------')
        return h

    def _eff_hamiltonian_list(self,J = None):
        """
        Returns a list of all effective hamiltonians and a second list with the degeneracies of the eigenvalues
        of these hamiltonians. The ordering of these to lists must be 1-to-1.
        :param J: Hopping parameter from Cayley Shells too outer Lieb shells (scale sqrt(k) values)
        :return:
        """

        Hs = []
        degeneracies = []
        # Symm States (once for l = 0 and l = 1)
        Hs.append(self._eff_hamiltonian_constructor(0,J))
        degeneracies.append(1)
        Hs.append(self._eff_hamiltonian_constructor(1,J))
        degeneracies.append(self.k) # Why is this not k+1?
        # Anti-Symm States (runs from 1 to mc with lc = 2l)
        for lc in range(1,self.mc + self.M%2):
            Hs.append(self._eff_hamiltonian_constructor(2*lc + 1,J))
            degeneracies.append((self.k -1)*(self.k + 1)*self.k**(lc-1))
        # print('degen:',degeneracies)
        return Hs, degeneracies



class DoubleLiebCayley(IsotropicAbstractTree):
    def __init__(self,M,k,J1=1,J2=1,J3=1,force_graph_object_creation=False):
        self.M = M
        self.mc = math.floor(M/3) # Number of Cayley Shells
        self.ml1 = math.floor(M/3) + math.ceil(M%3 /2) # Number of Lieb-1 Shells
        self.ml2 = math.floor(M/3) + math.floor(M%3 /2)
        self.k = k # Connectivity (Degree of each node = k + 1)
        self.forced_graph = force_graph_object_creation

        self.J1 = J1
        self.J2 = J2
        self.J3 = J3

        self.N = int(1 + ((k + 1)/(k-1))*(k**self.mc + k**self.ml1 + k**self.ml2 -3)) # Number of Nodes
        if self.N < 4000 or force_graph_object_creation:
            self.G = IsotropicAbstractTree._tree_creator(self.N,self.k,DoubleLiebCayley._tree_edges,J1=J1,J2=J2,J3=J3) # networkx graph object
            self._A = nx.adjacency_matrix(self.G) # Sparse Matrix


    @staticmethod
    def _tree_edges(n,r,J1=1,J2=1,J3=1):
        """
        Iteratively defines the tree structure
        :param n: Number of nodes
        :param r: Connectivity (k)
        :return: A list of tuples that define the edges of our tree
        """
        if n == 0:
            return
        # helper function for trees
        # yields edges in rooted tree at 0 with n nodes and branching ratio r
        nodes = iter(range(n))
        cayley_shell = [next(nodes)]  # stack of max length r
        lieb_shell1 = []
        lieb_shell2 = []
        r = r + 1
        first_run = True
        # Iterative filling of shells
        while cayley_shell:
            source = cayley_shell.pop(0)
            for i in range(r):
                try:
                    target = next(nodes)
                    lieb_shell1.append(target)
                    yield source, target, {'weight': J3}
                except StopIteration:
                    break
            if first_run:
                r = r - 1
                first_run = False

            if not cayley_shell:
                for lieb_node in lieb_shell1:
                    try:
                        target = next(nodes)
                        lieb_shell2.append(target)
                        yield lieb_node,target, {'weight': J1}
                    except StopIteration:
                        break
                lieb_shell1 = []

            if not cayley_shell:
                for lieb_node in lieb_shell2:
                    try:
                        target = next(nodes)
                        cayley_shell.append(target)
                        yield lieb_node,target, {'weight': J2}
                    except StopIteration:
                        break
                lieb_shell2 = []

    @property
    def A(self):
        return self._A.todense()

    def shell_list(self):
        """Returns a list of lists nodes, with each list containing the nodes of the same shell
        Makes use of the fact that the construction algorithm constructs each shell after the other

        Returns
        -------
        shell_lists : Multiple lists corresponding to the number of Shells M + 1 (incl. the center)
            each list containing the nodes of that shell
        """
        nodes = list(self.G.nodes)
        shell_lists = [[0]]
        l = 1
        for shell_number in range(1, self.M + 1):
            n = (self.k + 1) * (self.k ** (math.ceil(shell_number/3) - 1))
            shell_lists.append(nodes[l:(l + n)])
            l += n
        return shell_lists



    def _eff_hamiltonian_constructor(self,l,J = None,J_1=1,J_2=1):
        """
        Constructs the effective Hamiltonian of the Lieb-Cayley tree for a given shell number.
        Automatically adapts to M even or M odd cases. For M even, the matrix will have dimension M - l + 1.
        :param l: Number of the shell we're starting our construction from. For l = 0 we start from |0>
        :return: 2D-nparray of effective hamiltonian h and scalar of degeneracy d of the hamiltonian
        """
        if not J and J != 0:
            J = 1
        a = [np.sqrt(self.k)*J,1*J_1,1*J_2]
        offdiag = np.tile(a,reps=int(np.ceil(self.M/3)))
        offdiag[0] = np.sqrt(self.k + 1)*J

        h = np.diag(offdiag[l:self.M],k=1) + np.diag(offdiag[l:self.M],k=-1) #The M automatically cuts off the last element of the offdiagonal if M is odd

        return h

    def _eff_hamiltonian_list(self, J = None, J_1 = 1, J_2 = 1):
        """
        Returns a list of all effective hamiltonians and a second list with the degeneracies of the eigenvalues
        of these hamiltonians. The ordering of these to lists must be 1-to-1.
        :param J: Hopping parameter from Cayley Shells too outer Lieb shells (scale sqrt(k) values)
        :return:
        """

        Hs = []
        degeneracies = []
        # Symm States (once for l = 0 and l = 1)
        Hs.append(self._eff_hamiltonian_constructor(0,J,J_1,J_2))
        degeneracies.append(1)
        Hs.append(self._eff_hamiltonian_constructor(1,J,J_1,J_2))
        degeneracies.append(self.k) # Why is this not k+1?
        # Anti-Symm States (runs from 1 to mc with lc = 2l)
        for lc in range(1,self.mc + math.ceil(self.M%3 / 2)):
            Hs.append(self._eff_hamiltonian_constructor(3*lc + 1,J,J_1,J_2))
            degeneracies.append((self.k -1)*(self.k + 1)*self.k**(lc-1))
        # print('degen:',degeneracies)
        return Hs, degeneracies




class HusimiCayley(IsotropicAbstractTree):
    def __init__(self,M,k,force_graph_object_creation=False,circle=False):
        """
        A Husimi-Cayley tree with M shells and k children per node. A Husimi tree looks like a Cayley tree but the children  If N > 4000, the graph object will not be instantiated and certain
        of each node are fully connected. This imposes loops on the tree and at k=2 leads to a tree that tiles triangles.
        methods will not be available. To change this, you need to set force_graph_object_creation to true.
        :param M: int, number of shells
        :param k: int, number of children per node
        :param force_graph_object_creation: bool, If true, forces graph-object creation even if tree has more than 4000 nodes.
        :param circle: bool, if true places a circle at the center instead of a single node.
        """
        self.M = M
        self.k = k # Connectivity (Degree of each node = k + 1)
        if not circle:
            self.N = int(1 + k*(k + 1) * (k ** M - 1) /(k-1))
        if circle:
            self.N = int((k + 1) * (k ** M - 1) / (k - 1))
        self.circle = circle

        self.forced_graph = force_graph_object_creation

        if self.N < 4000 or force_graph_object_creation:
            self.G = IsotropicAbstractTree._tree_creator(self.N,self.k,HusimiCayley._tree_edges,circle=circle) # networkx graph object
            self._A = nx.adjacency_matrix(self.G) # Sparse Matrix


    @staticmethod
    def _tree_edges(n,r,circle=False):
        """
        Iteratively defines the tree structure
        :param n: Number of nodes
        :param r: Connectivity (k)
        :return: A list of tuples that define the edges of our tree
        """
        if n == 0:
            return
        # helper function for trees
        # yields edges in rooted tree at 0 with n nodes and branching ratio r
        nodes = iter(range(n)) # stack of max length r
        connector = []
        k = r
        r = int((r + 1)*r)

        if not circle:
            cayley_shell = [next(nodes)]  # stack of max length r
            connector = []
            first_run = True
        elif circle:
            cayley_shell = []
            first_run = False
            r = k
            for i in range(r + 1):
                cayley_shell.append(next(nodes))
            for connection in itertools.combinations(cayley_shell, 2):
                yield connection[0], connection[1]


        # Iterative filling of shells
        while cayley_shell:
            source = cayley_shell.pop(0)
            for i in range(r):
                try:
                    target = next(nodes)
                    connector.append(target)
                    cayley_shell.append(target)
                    yield source, target
                except StopIteration:
                    break
            if first_run: # The first node is connected to more nodes due to the triangle form but we only want the individual triangles to connect with each other
                r = k
                for sub_connector in [connector[i:i+r] for i in range(0, len(connector), r)]: #We only want the triangles to connect
                    for connection in itertools.combinations(sub_connector,2):
                        yield connection[0], connection[1]
                connector = []
                first_run = False
            for connection in itertools.combinations(connector,2):
                yield connection[0],connection[1]
            connector = []

    @property
    def A(self):
        return self._A.todense()

    def shell_list(self):
        """Returns a list of lists nodes, with each list containing the nodes of the same shell
        Makes use of the fact that the construction algorithm constructs each shell after the other

        Returns
        -------
        shell_lists : Multiple lists corresponding to the number of Shells M + 1 (incl. the center)
            each list containing the nodes of that shell
        """
        nodes = list(self.G.nodes)
        if not self.circle:
            shell_lists = [[0]]  # already contains the 0th node
            l = 1  # position of last added node
            for shell_number in range(1, self.M + 1):
                n = self.k * (self.k + 1) * (self.k ** (shell_number - 1))
                shell_lists.append(nodes[l:(l + n)])
                l += n
        elif self.circle:
            shell_lists = [[]]
            l = 0
            for shell_number in range(1, self.M + 1):
                n = (self.k + 1) * (self.k ** (shell_number - 1))
                shell_lists.append(nodes[l:(l + n)])
                l += n

        return shell_lists



    def _eff_hamiltonian_constructor(self,l,J = None):
        """
        Constructs the effective Hamiltonian of the Lieb-Cayley tree for a given shell number.
        Automatically adapts to M even or M odd cases. For M even, the matrix will have dimension M - l + 1.
        :param l: Number of the shell we're starting our construction from. For l = 0 we start from |0>
        :return: 2D-nparray of effective hamiltonian h and scalar of degeneracy d of the hamiltonian
        """
        if not J and J != 0:
            J = 1

        if not self.circle:
            if l == 0:
                h = np.eye(self.M+1,k=0)*(self.k-1) + np.eye(self.M+1,k=1)*np.sqrt(self.k)*J + np.eye(self.M+1,k=-1)*np.sqrt(self.k)*J
                h[0][0] = 0
                h[1][0] = np.sqrt(self.k*(self.k+1))*J
                h[0][1] = np.sqrt(self.k * (self.k + 1))*J

            elif l == 1:
                h = [np.eye(self.M+1-l,k=0)*(self.k-1) + np.eye(self.M+1-l,k=1)*np.sqrt(self.k)*J + np.eye(self.M+1-l,k=-1)*np.sqrt(self.k)*J]
                h[0][0][0] = -1
                h.append(np.eye(self.M+1-l,k=0)*(self.k-1) + np.eye(self.M+1-l,k=1)*np.sqrt(self.k)*J + np.eye(self.M+1-l,k=-1)*np.sqrt(self.k)*J)
            else:
                h = np.eye(self.M+1-l,k=0)*(self.k-1) + np.eye(self.M+1-l,k=1)*np.sqrt(self.k)*J + np.eye(self.M+1-l,k=-1)*np.sqrt(self.k)*J
                h[0][0] = -1
            return h

        else: # Circle Center -> Slightly different Hamiltonians
            if l == 0:
                h = np.eye(self.M, k=0) * (self.k - 1) + np.eye(self.M, k=1) * np.sqrt(self.k) * J + np.eye(
                    self.M, k=-1) * np.sqrt(self.k) * J
                h[0][0] = self.k

            else:
                h = np.eye(self.M + 1 - l, k=0) * (self.k - 1) + np.eye(self.M + 1 - l, k=1) * np.sqrt(
                    self.k) * J + np.eye(self.M + 1 - l, k=-1) * np.sqrt(self.k) * J
                h[0][0] = -1
            return h


    def _eff_hamiltonian_list(self,J = None):
        """
        Returns a list of all effective hamiltonians and a second list with the degeneracies of the eigenvalues
        of these hamiltonians. The ordering of these to lists must be 1-to-1.
        :param J: Hopping parameter from Cayley Shells too outer Lieb shells (scale sqrt(k) values)
        :return:
        """

        Hs = []
        degeneracies = []

        # Symm States (once for l = 0 and l = 1)
        Hs.append(self._eff_hamiltonian_constructor(0,J))
        degeneracies.append(1)
        if not self.circle:
            Hs.extend(self._eff_hamiltonian_constructor(1, J))
            degeneracies.append((self.k-1)*(self.k+1)) # anti-symm states that are on one branch only but their alpha is the central node
            degeneracies.append(self.k)# Additional degeneracy for anti-symm combinations of symmetric shell states
        else:
            Hs.append(self._eff_hamiltonian_constructor(1, J))
            degeneracies.append(self.k)
        # Anti-Symm States (runs from 1 to mc with lc = 2l)
        for l in range(2,self.M+1):
            Hs.append(self._eff_hamiltonian_constructor(l,J))
            degeneracies.append((self.k -1)*(self.k + 1)*self.k**(l-2))
        return Hs, degeneracies



class LiebHusimi(IsotropicAbstractTree):
    def __init__(self,M,k,J1=1,J2=1,force_graph_object_creation=False):
        self.M = M
        self.mc = math.floor(M/2) # Number of Cayley Shells
        self.ml = math.ceil(M/2) # Number of Lieb Shells
        self.k = k # Connectivity (Degree of each node = k + 1)
        logger.info('Initiating LiebHusimi Tree')
        self.N = int(((k + 1)/(k-1))*(k**self.mc + k**self.ml -2)) # Number of Nodes
        self.J1 = J1
        self.J2 = J2

        self.forced_graph = force_graph_object_creation

        if self.N < 4000 or force_graph_object_creation:
            self.G = IsotropicAbstractTree._tree_creator(self.N,self.k,LiebHusimi._tree_edges,J1=J1,J2=J2) # networkx graph object
            self._A = nx.adjacency_matrix(self.G) # Sparse Matrix


    @staticmethod
    def _tree_edges(n,r,J1,J2):
        """
        Iteratively defines the tree structure
        :param n: Number of nodes
        :param r: Connectivity (k)
        :return: A list of tuples that define the edges of our tree
        """
        if n == 0:
            return
        # helper function for trees
        # yields edges in rooted tree at 0 with n nodes and branching ratio r (Lieb-Husimi Structure)
        nodes = iter(range(n))

        cayley_shell = []  # stack of max length r
        lieb_shell = []

        # Construct the inner circle, treating it as a lieb-shell
        for i in range(r+1):
            lieb_shell.append(next(nodes))
        for connection in itertools.combinations(lieb_shell, 2):
            yield connection[0], connection[1], {'weight': J1}


        # Iterative filling of shells
        while lieb_shell:
            source = lieb_shell.pop(0)
            try:
                target = next(nodes)
                cayley_shell.append(target)
                yield source, target, {'weight': 1}
            except StopIteration:
                break

            if not lieb_shell:
                for cayley_node in cayley_shell:
                    triangle = []
                    for i in range(r):
                        try:
                            target = next(nodes)
                            triangle.append(target)
                            lieb_shell.append(target)
                            yield cayley_node,target, {'weight': J2}
                        except StopIteration:
                            break
                    for connection in itertools.combinations(triangle, 2):
                        yield connection[0], connection[1], {'weight': J1}
                cayley_shell = []

    @property
    def A(self):
        return self._A.todense()

    def shell_list(self):
        """Returns a list of lists nodes, with each list containing the nodes of the same shell
        Makes use of the fact that the construction algorithm constructs each shell after the other

        Returns
        -------
        shell_lists : Multiple lists corresponding to the number of Shells M + 1 (incl. the center)
            each list containing the nodes of that shell
        """
        nodes = list(self.G.nodes)
        shell_lists = [[]]
        l = 0
        for shell_number in range(1, self.M + 1):
            n = (self.k + 1) * (self.k ** (math.ceil(shell_number/2) - 1))
            shell_lists.append(nodes[l:(l + n)])
            l += n
        return shell_lists



    def _eff_hamiltonian_constructor(self,l):
        """
        Constructs the effective Hamiltonian of the Lieb-Cayley tree for a given shell number.
        Automatically adapts to M even or M odd cases. For M even, the matrix will have dimension M - l + 1.
        :param l: Number of the shell we're starting our construction from. For l = 0 we start from |0>
        :return: 2D-nparray of effective hamiltonian h and scalar of degeneracy d of the hamiltonian
        """

        a = [1,self.J2*np.sqrt(self.k)]
        offdiag = np.tile(a,reps=int(np.floor(self.M/2)))

        b = [self.J1*(self.k-1),0]
        diag = np.tile(b,reps=int(np.ceil(self.M/2)))

        if l == 0:
            diag[0] = self.k*self.J1
            h = [np.diag(diag[l:self.M]) + np.diag(offdiag[l:self.M-1],k=1) + np.diag(offdiag[l:self.M-1],k=-1)] #The M automatically cuts off the last element of the offdiagonal if M is odd
            diag[0] = -1*self.J1
            h.append(np.diag(diag[l:self.M]) + np.diag(offdiag[l:self.M-1],k=1) + np.diag(offdiag[l:self.M-1],k=-1))
        else:
            diag[0] = -1*self.J1
            h = np.diag(diag[0:self.M-l]) + np.diag(offdiag[0:self.M-1-l],k=1) + np.diag(offdiag[0:self.M-1-l],k=-1) #The M automatically cuts off the last element of the offdiagonal if M is odd

        return h

    def _eff_hamiltonian_list(self):
        """
        Returns a list of all effective hamiltonians and a second list with the degeneracies of the eigenvalues
        of these hamiltonians. The ordering of these to lists must be 1-to-1.
        :param J: Hopping parameter from Cayley Shells too outer Lieb shells (scale sqrt(k) values)
        :return:
        """

        Hs = []
        degeneracies = []
        # Symm States (there are two for lc = 0)
        Hs.extend(self._eff_hamiltonian_constructor(0))
        degeneracies.append(1), degeneracies.append((self.k)) # Uncertain about the second one here

        # Anti-Symm States (runs from 1 to mc with lc = 2l)
        for lc in range(1,self.mc + self.M%2):
            Hs.append(self._eff_hamiltonian_constructor(2*lc))
            degeneracies.append((self.k -1)*(self.k + 1)*self.k**(lc-1))
        return Hs, degeneracies





if __name__ == "__main__":
    HC = LiebHusimi(5,2,J1=2,J2=1)
    fig, ax_list = plt.subplots(2,1,sharex=True)

    eval, evec = HC.exact_diagonalization()
    ax_list[0].hist(eval,bins=201)
    ax_list[0].set_ylabel('D')
    ax_list[0].set_xlabel('E/t')
    ax_list[0].set_title('Exact Diagonalization Spectrum')

    eval2,weights2 = HC.effective_diagonalization()
    ax_list[1].hist(eval2,bins=201,weights=weights2)
    ax_list[1].set_ylabel('D')
    ax_list[1].set_xlabel('E/t')
    #ax_list[1].semilogy()
    ax_list[1].set_title('Effective Hamiltonian Diagonalization Spectrum')

    plt.show()