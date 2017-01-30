"""
uploads a notebook to the gallery
JJGC, February, 2017.
"""

import argparse
import sys
import os
from os import listdir
from os.path import isfile, join
import subprocess as sp
import datetime


def nb_to_html(nb_name):
    """Convert a nb to html
    1. nb_name is the name (full path and extension) of the nb to be converted
    2. the function runs the converter, copies the resulting html file in
    $IC_DATA and deletes the .html from the directory where the .ipynb is.
    """
    print('converting notebook {} to html'.format(nb_name))
    sp.run(["jupyter", "nbconvert", nb_name])

    nbb = nb_name.split('.ipynb')[0] + '.html'
    nbb2 = nbb.split('/')[-1]
    nbb_o = os.environ['IC_DATA'] + '/' + nbb2

    print('copying {} to {}'.format(nbb,nbb2))
    sp.run(["cp", nbb, nbb2])
    print('deleting {}'.format(nbb))
    sp.run(["rm", nbb])

    # print("""in nb_to_html():
    #          nbb = {}
    #          nbb2 = {}
    #          nbb_o = {}""".format(nbb, nbb2, nbb_o) )
    tag_nb = tag_html(nbb2)
    return tag_nb
#
def tag_html(nb_name):
    """Tag the html with the current date"""
    nb = nb_name
    print('tagging nb {}'.format(nb_name))

    now = datetime.datetime.now()
    tag = now.strftime("%Y-%m-%d_%H-%M")
    nbb = nb.split('.')[0] + '_' + tag + '.html'
    # print("""in tag_html():
    #          nb = {}
    #          nbb = {}
    #          """.format(nb, nbb) )
    print('copying {} to {}'.format(nb,nbb))
    sp.run(["cp", nb, nbb])
    print('deleting {}'.format(nb))
    sp.run(["rm", nb])
    return nbb

def upload_notebook(nb_name, path_to_nb, gallery_dir):
    """Upload notebook to the gallery"""

    tag_nb = nb_to_html(path_to_nb + nb_name)  #produces tagged nb in IC_DATA
    tag_nb_full_path = os.environ['IC_DATA'] + '/' + tag_nb
    path_to_gallery = '/data/icgallery/public/'
    nb_in_gallery = path_to_gallery + gallery_dir
    server_name = 'icuser@next.ific.uv.es'
    print("""uploading:
             nb name        = {}
             in path        = {}
             tagged as      = {}
             to gallery dir = {}
             """.format(nb_name, path_to_nb, tag_nb, gallery_dir) )
    print("""command:
             scp {} {}:{}
             """.format(tag_nb_full_path, server_name, nb_in_gallery) )

    sp.run(["scp", tag_nb_full_path, "{}:{}".format(server_name,
                                                    nb_in_gallery)])

if __name__ == "__main__":

    if os.environ['IC_DATA'] == '':
        print("""Environment variable IC_DATA not defined.
                 Please define your IC_DATA pointing to the local
                 directory where you will store data and html files.""")
        sys.exit(1)

    parser = argparse.ArgumentParser(description='gallery parser.')
    parser.add_argument("nb_name",
                        metavar="name of notebook", type=str,
                        help="name of notebook be converted to html")
    parser.add_argument("path_to_nb",
                        metavar="path to notebook", type=str,
                        help="path to notebook be converted to html")
    parser.add_argument("gallery_dir",
                        metavar="gallery directory", type=str,
                        help="name of gallery directory")

    args = parser.parse_args()
    # print(args)
    # print(args.gallery_dir)
    # print(args.path_to_nb)
    # print(args.nb_name)
    upload_notebook(args.nb_name, args.path_to_nb, args.gallery_dir)
