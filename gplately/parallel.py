"""Tools to execute routines efficiently by parallelising 
them over several threads. This uses multiple processing units.
"""
from multiprocessing import Pool, Process, Queue, cpu_count 


class Parallel(object):
    """A class that uses multiple processors with `multiprocessing` 
    to execute routines in parallel over several threads.
    
    Parameters
    -----------
    nprocs : int, default=1
        The number of separate executions of a process. By default,
        a single thread is run. 
    """
    def __init__(self, nprocs=1):

        self.nprocs = nprocs


    def parallelise_routine(self, function, *args, **kwargs):
        """Execute a routine over multiple threads on different 
        processors, ultimately reducing computation time.

        `parallelise_routine` permits one item through the process 
        queue when an executed item is extracted with get(). 

        Parameters
        ----------
        self.nprocs : int, default=1
            The number of separate executions of a process. By 
            default, a single thread is run.

        function : method from an instance of an object
            The process to be executed in parallel. Should be 
            supplied as module.class.method (if belonging to a class)
            or module.method. 

        *args : tuple
            Contains all necessary input parameters for the ‘function’.

        **kwargs : dict
            Keyword arguments for the ‘function’.
        """
        if self.nprocs == 1:
            # single thread
            result = function(*args, **kwargs)
            return result

        elif self.nprocs > 1:
            # more than one processor - game on

            results = [[] for i in range(n)]
            processes = []
            q_in = Queue(1)
            q_out = Queue()


            for i in range(self.nprocs):
                pass_args = [function]
                pass_args.extend(args)
                p = Process(target=self._func_queue, args=tuple(pass_args), kwargs=kwargs)
                processes.append(p)

            for p in processes:
                p.daemon = True
                p.start()

            # put items in the queue
            sent = [q_in.put((i,)) for i in range(n)]
            [q_in.put((None,)) for _ in range(nprocs)]

            # get the results
            results = []
            for i in range(len(sent)):
                index, result = q_out.get()
                results[index] = result


            # wait until each processor has finished
            [p.join() for p in processes]

            return results


    def _func_queue(self, function, q_in, q_out, *args, **kwargs):
        while True:
            pos, input_args = q_in.get()
            if pos is None:
                break

            res = function(*input_args, **kwargs)
            q_out.put((pos, res))
        return


