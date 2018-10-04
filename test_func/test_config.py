# test config load
# test config override
# test if config override is transfered to object loaded in different functions (if not so, then it is unusable)

import os
import unittest
from rna_blast_analyze.BR_core.config import CONFIG, tools_paths

fwd = os.path.dirname(__file__)


class TestConfig(unittest.TestCase):
    def test_config_import(self):
        for key in CONFIG.tool_paths.keys():
            if key == 'refold':
                # do not check
                pass
            else:
                self.assertEqual(CONFIG.tool_paths[key], '', key)

    def test_config_override(self):
        rfam_dir = '/test/test/Documents/rfamdb/'
        self.assertNotEqual(CONFIG.rfam_dir, rfam_dir)
        CONFIG.override(tools_paths(os.path.join(fwd, 'test_data', 'config_test.txt')))
        self.assertEqual(CONFIG.rfam_dir, rfam_dir)


if __name__ == '__main__':
    unittest.main()