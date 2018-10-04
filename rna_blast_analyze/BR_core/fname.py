import sys


def fname():
    return sys._getframe(1).f_code.co_name
