import random
import numpy as np
import fuzzynumber as fn

class Particle:

    def __init__(self, fitness, dim, seed, minx, maxx):
        self.rnd = random.Random(seed)
        self.dim = dim
        self.location = [0 for i in range(dim)]
        self.velocity = [0 for i in range(dim)]

        for i in range(dim):
            self.location[i] = (maxx[i] - minx[i]) * self.rnd.random() + minx[i]
            self.velocity[i] = (maxx[i] - minx[i]) * self.rnd.random() + minx[i]

        self.fitness = fitness(self.location)
        self.bestPos = [self.location[i] for i in range(dim)]
        self.bestFitness = self.fitness


class DBPSO:

    def __init__(self, f, dim, N, minx, maxx, Maxgen=100, w=0.7298, mu1=1.49618, mu2=1.49618):
        self.mu1 = mu1
        self.mu2 = mu2
        self.f = f
        self.N = N
        self.w = w
        self.Maxgen = Maxgen
        self.dim = dim
        self.minx = minx
        self.maxx = maxx
        if not (len(minx) == len(maxx) == dim):
            raise Exception("Dimensions of minimum and maximum arrays must equal dimension of solution space")

    def solve(self):
        w = self.w
        fitness = self.f
        c1 = self.mu1
        c2 = self.mu2
        N = self.N
        dim = self.dim
        Maxgen = self.Maxgen
        minx = self.minx
        maxx = self.maxx

        swarm = [Particle(fitness, dim, i, minx, maxx) for i in range(N)]
        bestSwarmPos = [0 for i in range(dim)]
        bestSwarmFitness = 100000000000

        for i in range(N):
            if swarm[i].fitness < bestSwarmFitness:
                bestSwarmFitness = swarm[i].fitness
                bestSwarmPos = [swarm[i].location[j] for j in range(dim)]

        for Iter in range(Maxgen):
            for i in range(N):
                for k in range(dim):
                    r = np.random.rand(2)
                    swarm[i].velocity[k] = (w * swarm[i].velocity[k] +
                                           (c1 * r[0] * (swarm[i].bestPos[k] - swarm[i].location[k])) +
                                            (c2 * r[1] * (bestSwarmPos[k] - swarm[i].location[k])))

                    if swarm[i].velocity[k] < minx[k]:
                        swarm[i].velocity[k] = minx[k]
                    elif swarm[i].velocity[k] > maxx[k]:
                        swarm[i].velocity[k] = maxx[k]

                for k in range(dim):
                    swarm[i].location[k] += swarm[i].velocity[k]

                swarm[i].fitness = fitness(swarm[i].location)

                if swarm[i].fitness < swarm[i].bestFitness:
                    swarm[i].bestFitness = swarm[i].fitness
                    swarm[i].bestPos = [swarm[i].location[j] for j in range(dim)]

                if swarm[i].fitness < bestSwarmFitness:
                    bestSwarmPos = [swarm[i].location[j] for j in range(dim)]
                    bestSwarmFitness = swarm[i].fitness

        return bestSwarmPos, bestSwarmFitness


