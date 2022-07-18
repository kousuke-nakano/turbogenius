#!usr/bin/env python
# coding: utf-8

"""

TurboGenius command line interface

"""

import sys, os
import click

# Thanks Oto :-))
class OptionEatAll(click.Option):

    def __init__(self, *args, **kwargs):
        self.save_other_options = kwargs.pop('save_other_options', True)
        nargs = kwargs.pop('nargs', -1)
        assert nargs == -1, 'nargs, if set, must be -1 not {}'.format(nargs)
        super(OptionEatAll, self).__init__(*args, **kwargs)
        self._previous_parser_process = None
        self._eat_all_parser = None

    def add_to_parser(self, parser, ctx):
        def parser_process(value, state):
            # method to hook to the parser.process
            done = False
            value = [value]
            if self.save_other_options:
                # grab everything up to the next option
                while state.rargs and not done:
                    for prefix in self._eat_all_parser.prefixes:
                        if state.rargs[0].startswith(prefix):
                            done = True
                    if not done:
                        value.append(state.rargs.pop(0))
            else:
                # grab everything remaining
                value += state.rargs
                state.rargs[:] = []
            value = tuple(value)

            # call the actual process
            self._previous_parser_process(value, state)

        retval = super(OptionEatAll, self).add_to_parser(parser, ctx)
        for name in self.opts:
            our_parser = parser._long_opt.get(name) or parser._short_opt.get(name)
            if our_parser:
                self._eat_all_parser = our_parser
                self._previous_parser_process = our_parser.process
                our_parser.process = parser_process
                break
        return retval

def header(f):
    @click.option("-log","log_level", default="INFO", help= 'logger level, DEBUG, INFO, ERROR')
    def ret(*args, **kwargs):
        logger = getLogger("Turbo-Genius")
        logger.setLevel(kwargs["log_level"])
        stream_handler = StreamHandler()
        #handler_format = Formatter('%(message)s')
        handler_format = Formatter('%(name)s - %(levelname)s - %(lineno)d - %(message)s')
        stream_handler.setFormatter(handler_format)
        logger.addHandler(stream_handler)

        logger_p = getLogger("pyturbo")
        logger_p.setLevel(kwargs["log_level"])
        stream_handler_p = StreamHandler()
        #handler_format_p = Formatter('%(message)s')
        handler_format_p = Formatter('%(name)s - %(levelname)s - %(lineno)d - %(message)s')
        stream_handler_p.setFormatter(handler_format_p)
        logger_p.addHandler(stream_handler_p)

        operation = ",".join( [ val for val, com in zip(["generate", "run", "postprocess"],
                                                        [       "g",   "r",        "post"])
                                if com in kwargs and kwargs[com] ] )
        kwargs["operation"] = operation
        if any([ x in kwargs for x in ["r", "g", "post"]]):
            if operation == "":
                logger.info(f"Jobtype is needed please type --help for more information")
                logger.info(f"")

            else:
                logger.info(f"Selected jobtypes = {kwargs['operation']}")
        try:
            from icecream import ic
            ic(kwargs)
        except:
            pass
        f(*args, **kwargs)
    ret.__name__ = f.__name__
    ret.__doc__ = f.__doc__

    return ret

def decorate_grpost(f):
    f = click.option("-g", "g", is_flag=True, help="Generate an input file")(f)
    f = click.option("-r", "r", is_flag=True, help="Run a program")(f)
    f = click.option("-post", "post", is_flag=True, help="Postprocess")(f)
    return f

@click.group()
def cli():
    pass

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from makefort10_genius import *
from convertfort10mol_genius import *
from prep_genius import *
from vmc_genius import *
from lrdmc_genius import *
from vmc_opt_genius import *
from lrdmc_opt_genius import *
from correlated_sampling_genius import *
from wavefunction import *
from commandline_tools import *

if __name__ == "__main__":
    cli()
