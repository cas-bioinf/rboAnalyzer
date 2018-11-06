import time


def timeit_decorator(func):
    def wrapper(*arg, **kwargs):
        t1 = time.time()
        res = func(*arg, **kwargs)
        t2 = time.time()
        # print '%s took %0.3f ms' % (func.func_name, (t2-t1)*1000.0)
        exec_time = t2 - t1
        return res, exec_time
    return wrapper


def tryit_decorator(func):
    def wrapper(*arg, **kwargs):
        try:
            res = func(*arg, **kwargs)
            return res
        except:
            print('problem in prediction {}'.format(str(func)))
            raise
    return wrapper