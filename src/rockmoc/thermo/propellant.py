"""Propellant → (γ, R, Tc) for a perfect-gas nozzle design (RocketCEA if installed)."""

from __future__ import annotations

import contextlib
from functools import lru_cache

PSI = 6894.757293168


@lru_cache(maxsize=256)
def of_props(ox: str, fuel: str, of: float, pc: float, eps: float = 40.0) -> dict:
    """Chamber/throat properties at mixture ratio ``of`` and chamber pressure ``pc`` [Pa].

    Returns ``{gamma, R, Tc, cea}``. Without RocketCEA (``pip install rocketcea``) an
    engineering stand-in is used and ``cea`` is False — fine for trends, not for design.
    """
    try:
        from rocketcea.cea_obj import CEA_Obj

        c = CEA_Obj(oxName=ox, fuelName=fuel)
        p = pc / PSI
        mw, gamma = c.get_Chamber_MolWt_gamma(p, of, eps)
        with contextlib.suppress(Exception):  # ponytail: CEA throat call is optional
            gamma = c.get_Throat_MolWt_gamma(p, of, eps)[1]
        return {"gamma": float(gamma), "R": 8314.462618 / float(mw),
                "Tc": float(c.get_Tcomb(p, of)) / 1.8, "cea": True}  # CEA returns °R
    except Exception:  # noqa: BLE001  # ponytail: any CEA failure → stand-in
        from rockmoc.plot.sweeps import _synthetic_of_props

        g, R, Tc = _synthetic_of_props(of)
        return {"gamma": g, "R": R, "Tc": Tc, "cea": False}
