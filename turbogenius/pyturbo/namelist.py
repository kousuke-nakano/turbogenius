#!python -u
# -*- coding: utf-8 -*-

"""

pyturbo: namelist related classes and methods

Todo:
    * docstrings are not completed.
    * refactoring assert sentences. The assert should not be used for any on-the-fly check.
    * implementing __str__ method.
    * implementing sanity_check method.

"""

# python modules
import re
from typing import Optional, Union
from collections import OrderedDict

# turbo-genius modules
from turbogenius.pyturbo.utils.utility import get_str_variable_type_auto

# set logger
from logging import getLogger, StreamHandler, Formatter

logger = getLogger("pyturbo").getChild(__name__)


class Namelist:
    def __init__(self, namelist: Optional[OrderedDict] = None):
        if namelist is None:
            namelist = OrderedDict()
        assert isinstance(
            namelist, OrderedDict
        ), "Please use OrderedDict() for the namelist!"
        self.__namelist = namelist

    @property
    def parameters(self):
        return self.__namelist

    def set_parameter(
        self,
        parameter: str,
        value: Union[int, float, str],
        namelist: Optional[dict] = None,
    ):
        if namelist is None:
            for key, parameters in self.__namelist.items():
                if parameter in parameters.keys():
                    self.__namelist[key][parameter] = value
                    return True
            raise KeyError(
                f"parameter={parameter} is not defined in the defined namelist. Specify a namelist."
            )

        else:  # namelist is not None
            try:
                self.__namelist[namelist][parameter] = value
            except KeyError:
                raise (f"{namelist} is not defined in the namelist.")

    def get_parameter(self, parameter, namelist=None):
        for key, parameters in self.__namelist.items():
            if namelist is not None and key != namelist:
                continue
            if parameter in parameters.keys():
                return parameters[parameter]

        return None

    def get_parameters(self):
        return self.__namelist

    def comment_out(self, parameter):
        for key, parameters in self.__namelist.items():
            if parameter in parameters.keys():
                value = self.__namelist[key].pop(parameter)
                self.__namelist[key]["!" + parameter] = value

    def write(self, file_name):
        output = []

        for key, parameters in self.__namelist.items():
            output.append(f"{key}\n")
            namelist_values = parameters
            for key, value in namelist_values.items():
                if (
                    type(value) == str
                    and not re.match(".*true.*", value)
                    and not re.match(".*false.*", value)
                ):
                    output.append(f"    {key}='{value}'\n")
                else:
                    output.append(f"    {key}={value}\n")
            output.append("/\n\n")

        with open(file_name, "w") as f:
            f.writelines(output)

    @staticmethod
    def read_parameters_from_file(file_name):
        with open(file_name, "r") as f:
            input_lines = f.readlines()

        namelist = []
        namelist_l = []
        namelist_d = {}
        namelist_d_ordered = OrderedDict()
        read_flag = False

        for line in input_lines:
            if re.match("^[\s]*!.*", line):
                continue
            if re.match(".*ATOMIC_POSITIONS.*", line):
                break
            if re.match(".*&.*", line):
                namelist.append(line.replace("\n", ""))
                namelist_d = {}
                read_flag = True
                continue
            if re.match(".*/.*", line):
                namelist_l.append(namelist_d)
                read_flag = False
                continue
            if read_flag:
                # print(line)
                # key,value,_= line.replace("\n","").replace(" ","").split("=")
                key, value, *_ = re.split(
                    "[=,!]", line.replace("\n", "").replace(" ", "")
                )
                converted_value = get_str_variable_type_auto(value)
                namelist_d[key] = converted_value

        for name, parameters in zip(namelist, namelist_l):
            namelist_d_ordered[name] = parameters

        return namelist_d_ordered

    @classmethod
    def parse_namelist_from_file(cls, file_name):
        namelist_d_ordered = cls.read_parameters_from_file(file_name)
        return cls(namelist=namelist_d_ordered)


if __name__ == "__main__":
    logger = getLogger("pyturbo")
    logger.setLevel("INFO")
    stream_handler = StreamHandler()
    stream_handler.setLevel("DEBUG")
    handler_format = Formatter(
        "%(name)s - %(levelname)s - %(lineno)d - %(message)s"
    )
    stream_handler.setFormatter(handler_format)
    logger.addHandler(stream_handler)

    # moved to examples
