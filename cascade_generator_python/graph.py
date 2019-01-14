from csv import reader

from pandas import read_csv

from time import time


class Graph:
    def __init__(self):
        self.users = set()
        self.user_community = {}
        self.latent_factors = set()
        self.in_net = None
        self.in_degree = {}
        self.out_net = None
        self.out_degree = {}
        self.n_users = 0
        self.n_factors = 0

    def _build_direct_net(self, groups, id_on_first_column):
        net = {}

        for group in groups:
            key = group[0]

            if id_on_first_column:
                value = group[1][1].values
            else:
                value = group[1][0].values

            self.users.add(key)
            net[key] = value

        return net

    def read_data(self, _edge_list_file, _community_file):
        t = time()

        print('Building graph...')

        print('\tReading data...')

        df = read_csv(_edge_list_file, sep='\t', header=None)

        self.out_net = self._build_direct_net(df.groupby(by=0), True)
        self.in_net = self._build_direct_net(df.groupby(by=1), False)

        del df

        print('\tLoading communities...')

        with open(_community_file, 'rt') as f:
            reader_obj = reader(f, delimiter="\t")

            for line in reader_obj:
                community = int(line[1])
                self.latent_factors.add(community)
                self.user_community[int(line[0])] = community

        del reader_obj

        print('\tComputing degrees...')

        for user in self.users:
            self.in_degree[user] = len(self.in_net[user])
            self.out_degree[user] = len(self.out_net[user])

        self.n_users = len(self.users)
        self.n_factors = len(self.latent_factors)

        print('Building graph... complete! Elapsed time: {0:.3f} seconds'.format(time() - t))
