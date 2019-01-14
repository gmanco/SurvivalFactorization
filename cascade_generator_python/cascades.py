from graph import Graph

from numpy import arange, zeros
from numpy.random import choice, exponential, poisson, randint, uniform

from time import time
from tqdm import tqdm


def generate_cascades(n_cascades, g: Graph, a, s, phi, t_min, t_max, debug=False, min_size=5):
    t = time()

    print('Building cascades...')

    assignments = {}
    cascades = {}
    content = {}

    tot_users = arange(g.n_users, dtype=int)

    norm_a = zeros(shape=a.shape)

    for k in range(g.n_factors):
        norm_a[:, k] = a[:, k] / sum(a[:, k])

    built_cascades = 0

    with tqdm(total=n_cascades) as progress_bar:
        for c in range(n_cascades):
            t2 = time()

#            if debug:
#                print('\rBuilding cascade', c + 1, end='', flush=True)

            cascade = {}
            triggers = set()

            k = randint(low=0, high=g.n_factors)
            assignments[c] = k

            cur_user = choice(g.n_users, p=norm_a[:, k])
            n_active_nodes = 1

            timestamp = uniform(low=0, high=t_min)
            cascade[cur_user] = timestamp
            can_expand = True

            while can_expand:
                cur_users = tot_users[list(set(cascade.keys()) - triggers)]

                if len(cur_users) > 0:
                    if len(cur_users) > 1:
                        p = a[cur_users, k] / sum(a[cur_users, k])
                        cur_user = choice(cur_users, p=p)
                    else:
                        cur_user = cur_users[0]

                    triggers.add(cur_user)
                    nb = g.in_net[cur_user + 1] - 1

                    for next_user in nb:
                        try:
                            x = cascade[next_user]
                        except KeyError:
                            x = None

                        if x is None:
                            rate = s[next_user, k] * a[cur_user, k]
                            timestamp = cascade[cur_user] + exponential(1 / rate)

                            if timestamp < t_max:
                                cascade[next_user] = timestamp
                                n_active_nodes += 1
                else:
                    can_expand = False

            size = len(cascade)

            if size >= min_size:
                cascades[c] = cascade
                content[c] = poisson(lam=phi[k])
                built_cascades += 1


                strOut = 'size {} - elapsed time {:.3f} seconds [{}/{}]'.format(len(cascade),time() - t2,built_cascades,c)
            else:
                strOut = 'too short - elapsed time {:.3f} seconds[{}/{}]'.format(time() - t2,built_cascades,c)

            progress_bar.set_postfix(status=strOut)
            progress_bar.update()
    progress_bar.close()

    print('Building cascades... complete! built ', built_cascades,
          ' cascades. Elapsed time: {0:.3f} seconds'.format(time() - t))

    return assignments, cascades, content
