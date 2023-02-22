#!python -u
# -*- coding: utf-8 -*-

from __future__ import print_function

# set logger
from logging import getLogger

logger = getLogger("pyturbo").getChild(__name__)


def prompt(text, checker):
    """Loop input() *FOREVER* while input is invalid."""
    while True:
        output = input(text)
        if checker(output):
            return output


def get_nonlocalmoves_setting(nonlocalmoves: str):
    if nonlocalmoves == "tmove":
        typereg = 0
        npow = 0.0
    elif nonlocalmoves == "dla":
        typereg = 6
        npow = 1.0
    elif nonlocalmoves == "dlatm":
        typereg = 6
        npow = 0.0
    elif nonlocalmoves == "la":
        typereg = 0
        npow = 1.0
    else:
        logger.error(f"nonlocalmoves={nonlocalmoves} is not implemented.")
        raise NotImplementedError

    return typereg, npow


def get_optimizer_flags(
    optimizer: str = "sr",
    opt_onebody: bool = True,
    opt_twobody: bool = True,
    opt_det_mat: bool = True,
    opt_jas_mat: bool = True,
    opt_det_basis_exp: bool = True,
    opt_jas_basis_exp: bool = True,
    opt_det_basis_coeff: bool = True,
    opt_jas_basis_coeff: bool = True,
    qmc_type: str = "vmc",
):
    logger.debug(optimizer)

    # check if either onebody or twobody jastrow is optimized
    if opt_onebody or opt_twobody:
        iesd = 1
        if opt_onebody:
            iesdonebodyoff = ".false."
        else:
            iesdonebodyoff = ".true."
        if opt_twobody:
            iesdtwobodyoff = ".false."
        else:
            iesdtwobodyoff = ".true."
    else:
        iesd = 0
        iesdonebodyoff = ".false."
        iesdtwobodyoff = ".false."

    # jastrow matrix element optimization
    if opt_jas_mat:
        iesfree = 1
        twobodyoff = ".false."
    else:
        if opt_onebody:
            iesfree = 1
            twobodyoff = ".true."
        else:
            iesfree = 0
            twobodyoff = ".false."

    # determinant matrix element optimization
    if opt_det_mat:
        iessw = 1
    else:
        iessw = 0

    # det. basis set optimization
    if opt_det_basis_coeff or opt_det_basis_exp:
        iesup = 1
        if opt_det_basis_coeff and not opt_det_basis_exp:  # opt only coeff.
            if optimizer == "sr":  # stochastic reconfiguration
                optimizer_number = return_optimizer_number(
                    optimizer=optimizer, qmc_type=qmc_type, opt_exponent=False
                )
            elif optimizer == "lr":  # linear method
                optimizer_number = return_optimizer_number(
                    optimizer=optimizer, qmc_type=qmc_type, opt_exponent=False
                )
            else:
                raise NotImplementedError
        else:  # opt also exponents
            if optimizer == "sr":  # stochastic reconfiguration
                optimizer_number = return_optimizer_number(
                    optimizer=optimizer, qmc_type=qmc_type, opt_exponent=True
                )
            elif optimizer == "lr":  # linear method
                optimizer_number = return_optimizer_number(
                    optimizer=optimizer, qmc_type=qmc_type, opt_exponent=True
                )
            else:
                raise NotImplementedError
    else:
        iesup = 0
        if optimizer == "sr":  # stochastic reconfiguration
            optimizer_number = return_optimizer_number(
                optimizer=optimizer, qmc_type=qmc_type, opt_exponent=False
            )
        elif optimizer == "lr":  # linear method
            optimizer_number = return_optimizer_number(
                optimizer=optimizer, qmc_type=qmc_type, opt_exponent=False
            )
        else:
            raise NotImplementedError

    # jas. basis set optimization
    if opt_jas_basis_coeff or opt_jas_basis_exp:
        iesm = 1
        if opt_jas_basis_coeff and not opt_jas_basis_exp:  # opt only coeff.
            if optimizer == "sr":  # stochastic reconfiguration
                assert optimizer_number == return_optimizer_number(
                    optimizer=optimizer, qmc_type=qmc_type, opt_exponent=False
                ), "There is a conflict!! Not implemented."
            elif optimizer == "lr":  # linear method
                assert optimizer_number == return_optimizer_number(
                    optimizer=optimizer, qmc_type=qmc_type, opt_exponent=False
                ), "There is a conflict!! Not implemented."
            else:
                raise NotImplementedError
        else:
            if optimizer == "sr":  # stochastic reconfiguration
                assert optimizer_number == return_optimizer_number(
                    optimizer=optimizer, qmc_type=qmc_type, opt_exponent=True
                ), "There is a conflict!! Not implemented."
            elif optimizer == "lr":  # linear method
                assert optimizer_number == return_optimizer_number(
                    optimizer=optimizer, qmc_type=qmc_type, opt_exponent=True
                ), "There is a conflict!! Not implemented."
            else:
                raise NotImplementedError
    else:
        iesm = 0

        if optimizer == "sr":  # stochastic reconfiguration
            assert optimizer_number == return_optimizer_number(
                optimizer=optimizer, qmc_type=qmc_type, opt_exponent=False
            ), "There is a conflict!! Not implemented."
        elif optimizer == "lr":  # linear method
            assert optimizer_number == return_optimizer_number(
                optimizer=optimizer, qmc_type=qmc_type, opt_exponent=False
            ), "There is a conflict!! Not implemented."
        else:
            raise NotImplementedError

    return (
        optimizer_number,
        iesdonebodyoff,
        iesdtwobodyoff,
        twobodyoff,
        iesd,
        iesfree,
        iessw,
        iesup,
        iesm,
    )


def return_optimizer_number(
    optimizer: str = "sr", qmc_type: str = "vmc", opt_exponent: bool = False
):
    if optimizer == "lr":
        if qmc_type == "vmc":
            if opt_exponent:
                optimizer_number = -8
            else:
                optimizer_number = -4
        elif qmc_type == "lrdmc":
            if opt_exponent:
                optimizer_number = -28
            else:
                optimizer_number = -24
        else:
            raise NotImplementedError

    elif optimizer == "sr":
        if qmc_type == "vmc":
            if opt_exponent:
                optimizer_number = -5
            else:
                optimizer_number = -9
        elif qmc_type == "lrdmc":
            if opt_exponent:
                optimizer_number = -25
            else:
                optimizer_number = -29
        else:
            raise NotImplementedError
    else:
        raise NotImplementedError

    return optimizer_number
