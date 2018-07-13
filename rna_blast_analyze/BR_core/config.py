import configparser
import os
import re
from distutils.util import strtobool
from subprocess import check_output, CalledProcessError, STDOUT
# Beware of lazy or conditional import with config!!!! It would overwrite custom config file.


class tools_paths(object):
    def __init__(self, config_file):
        self.tool_paths = {
            'refold': os.path.abspath(
                os.path.dirname(__file__) + '/../3rd_party_source/refold/'
            ) + os.path.sep,
            'infernal': '',
            'muscle': '',
            'clustal': '',
            'locarna': '',
            'viennarna_bin': '',
            'mfold': '',
            'blast': '',
            'mafft': '',
            't-coffee': '',
            'centroid': '',
            'turbofold': '',
            'rnashapes': '',
            'rapidshapes': '',

        }

        self.data_paths = {
            'rsearch_ribosum': os.path.abspath(
                os.path.dirname(__file__) + '/../3rd_party_source/RSEARCH_matrices/RIBOSUM65.mat'
            ),
            'rfam_dir': os.path.abspath(
                os.path.dirname(__file__) + '/../3rd_party_source/rfamdb'
            ),
            'rfam_url': 'ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz',
            'html_template_dir': os.path.abspath(
                os.path.dirname(__file__) + '/output'
            )
        }

        self.ssh = {
            'ssh_user': '',
            'ssh_pass': '',
            'use_virtual': False,
        }

        if config_file:
            conf_file = config_file
        else:
            conf_file = os.path.sep.join(__file__.split(os.path.sep)[:-1] + ['config.txt'])

        config = configparser.ConfigParser()
        config.read(conf_file)

        if 'TOOL_PATHS' in config:
            for key in config['TOOL_PATHS'].keys():
                if config['TOOL_PATHS'][key].endswith(('\\', '/')):
                    self.tool_paths[key] = config['TOOL_PATHS'][key]
                else:
                    self.tool_paths[key] = config['TOOL_PATHS'][key] + os.sep

        if 'SSH' in config:
            for key in config['SSH'].keys():
                if key == 'use_virtual':
                    self.ssh[key] = strtobool(config['SSH'][key])
                else:
                    self.ssh[key] = config['SSH'][key]

        if 'DATA' in config:
            for key in config['DATA'].keys():
                self.data_paths[key] = config['DATA'][key]

    def override(self, other):
        self.differ(self.tool_paths, other.tool_paths)
        self.differ(self.ssh, other.ssh)
        self.differ(self.data_paths, other.data_paths)

    @staticmethod
    def differ(a, b):
        for key in a.keys():
            if a[key] != b[key]:
                a[key] = b[key]

    @property
    def refold_path(self):
        return self.tool_paths['refold']

    @property
    def infernal_path(self):
        return self.tool_paths['infernal']

    @property
    def muscle_path(self):
        return self.tool_paths['muscle']

    @property
    def clustal_path(self):
        return self.tool_paths['clustal']

    @property
    def locarna_path(self):
        return self.tool_paths['locarna']

    @property
    def viennarna_path(self):
        return self.tool_paths['viennarna_bin']

    @property
    def mfold_path(self):
        return self.tool_paths['mfold']

    @property
    def blast_path(self):
        return self.tool_paths['blast']

    @property
    def mafft_path(self):
        return self.tool_paths['mafft']

    @property
    def tcoffee_path(self):
        return self.tool_paths['t-coffee']

    @property
    def centriod_path(self):
        return self.tool_paths['centroid']

    @property
    def turbofold_path(self):
        return self.tool_paths['turbofold']

    @property
    def rnashapes_path(self):
        return self.tool_paths['rnashapes']

    @property
    def rapidshapes_path(self):
        return self.tool_paths['rapidshapes']

    @property
    def SSH_USER(self):
        return self.ssh['ssh_user']

    @property
    def SSH_PASS(self):
        return self.ssh['ssh_pass']

    @property
    def USE_VIRTUAL(self):
        return self.ssh['use_virtual']

    @property
    def rfam_dir(self):
        return self.data_paths['rfam_dir']

    @property
    def rfam_url(self):
        return self.data_paths['rfam_url']

    @property
    def rsearch_ribosum(self):
        return self.data_paths['rsearch_ribosum']

    @property
    def html_template_dir(self):
        return self.data_paths['html_template_dir']


# load defaults
CONFIG = tools_paths(config_file=None)

if __name__ == '__main__':
    CC = tools_paths(config_file=None)
    for s in CC.tool_paths.values():
        print(s)
    for s in CC.ssh.values():
        print(s)
    for s in CC.rfam.values():
        print(s)
