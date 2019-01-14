from graph import Graph

from numpy import zeros
from numpy.random import binomial, choice, gamma, normal, rand, uniform

from time import time


# FIXME sostituire i '-1' con un mapping sugli id di utenti e topics


def _degree_to_normalized_vector(users, degree_dict):
    degree = zeros(shape=(len(users),))
    maximum = -1

    for user in users:
        degree[user - 1] = degree_dict[user]

        if maximum < degree[user - 1]:
            maximum = degree[user - 1]

    if maximum > 0:
        degree /= maximum

    return degree


def build_matrices(g: Graph, bias, blshape, blscale, tshape, tscale, n_words, prob=0.9):
    t = time()

    print('Building matrices...')

    in_degrees = _degree_to_normalized_vector(g.users, g.in_degree)
    out_degrees = _degree_to_normalized_vector(g.users, g.out_degree)

    max_val = min(min(1 - out_degrees[out_degrees < 1]), min(in_degrees[in_degrees > 0])) / 10
    min_val = max_val / 100

    s = rand(g.n_users, g.n_factors) * (max_val - min_val) + min_val
    a = rand(g.n_users, g.n_factors) * (max_val - min_val) + min_val

    for user in g.users:
        community = g.user_community[user]

        if g.out_degree[user] > 0:
            p = binomial(1, prob)

            if p == 1:
                s[user - 1, community - 1] = max(0.01, abs(normal(1 - out_degrees[user - 1], 0.05)))
            else:
                s[user - 1, community - 1] = uniform(low=0.1, high=1)

        if g.in_degree[user] > 0:
            p = binomial(1, prob)

            if p == 1:
                a[user - 1, community - 1] = max(0.01, abs(normal(in_degrees[user - 1], 0.05)))
            else:
                a[user - 1, community - 1] = uniform(low=0.1, high=1)

    p = gamma(blshape, blscale, size=n_words * g.n_factors).reshape((g.n_factors, n_words))
    probs = zeros(shape=(g.n_factors + 1))
    probs[0] = bias
    probs[1:] = (1 - bias) / g.n_factors
    assignments = choice(g.n_factors + 1, n_words, p=probs)

    for factor in g.latent_factors:
        active_words = assignments == factor  # factor - 1 + 1 per via del bias
        p[factor - 1, active_words] = gamma(tshape, tscale, size=sum(active_words))

    print('\tA.shape =', a.shape)
    print('\tS.shape =', s.shape)
    print('\tPhi.shape =', p.shape)

    print('Building matrices... complete! Elapsed time: {0:.3f} seconds'.format(time() - t))

    return a, s, p
