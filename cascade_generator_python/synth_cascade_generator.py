from cascades import generate_cascades
from graph import Graph
from matrices import build_matrices
from numpy import nonzero
from pandas import DataFrame
from time import time
from optparse import OptionParser


def print_cascades(_cascades, output_file):
    data = {'User': [], 'Item': [], 'TimeStamp': []}

    for cascade_id in _cascades:
        cascade = _cascades[cascade_id]

        for user in cascade:
            data['User'].append(user + 1)
            data['Item'].append(cascade_id + 1)
            data['TimeStamp'].append(cascade[user])

    DataFrame(data).sort_values(by=['Item', 'TimeStamp']).to_csv(output_file+'.bz2', sep='\t', index=False, compression='bz2')


def print_content(_content, output_file):
    data = {'Word': [], 'Item': [], 'Frequency': []}

    for cascade_id in _content:
        text = _content[cascade_id]
        idxs = nonzero(text)[0]

        for idx in idxs:
            data['Word'].append(idx + 1)
            data['Item'].append(cascade_id + 1)
            data['Frequency'].append(text[idx])

    DataFrame(data).sort_values(by=['Item', 'Word']).to_csv(output_file+'.bz2', sep='\t', index=False, compression='bz2')


if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option("-n", "--networks", dest="network", default="network.dat",
                      help="network file name")

    parser.add_option("-c", "--community", dest="community", default="community.dat",
                      help="community file name")

    parser.add_option("-p", "--prefix", dest="prefix", default="output",
                      help="prefix to use")

    parser.add_option("-s", "--cascadesize", dest="cascadesize", type="int", default=1024,
                      help="Number of cascades to generate")

    parser.add_option("-t", "--tmax", dest="tmax", type="int", default=100,
                      help="Time horizon")


    (options, args) = parser.parse_args()
    dict_opts = vars(options)


    t = time()
    print('Starting the procedure...')

    # Parameters ###################################
    edge_list_file = dict_opts['network'] #'data/s2-1k/network.dat'
    community_file = dict_opts['community']#'data/s2-1k/community.dat'
    prefix = dict_opts['prefix']#'data/s2-1k/s2-1k'

    bias = 0.5
    blshape = 1
    blscale = 1
    tshape = 3
    tscale = 3
    n_words = 1024
    n_cascades = dict_opts['cascadesize']# 2 ** 13
    t_min = 5
    t_max = dict_opts['tmax']  #100
    # ##############################################

    g = Graph()
    g.read_data(edge_list_file, community_file)

    print()

    a, s, p = build_matrices(g, bias, blshape, blscale, tshape, tscale, n_words)

    print()

    _, cascades, content = generate_cascades(n_cascades, g, a, s, p, t_min, t_max, debug=True)

    print_cascades(cascades, prefix + '.cascades')
    print_content(content, prefix + '.content')

    print('... procedure complete! Elapsed time: {0:.3f} seconds'.format(time() - t))
