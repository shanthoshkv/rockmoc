"""Built-in example gases for students (engineering frozen estimates)."""

from __future__ import annotations

from rockmoc.thermo.perfect_gas import PerfectGas

# name -> (gamma, R [J/kg/K], short note)
_PRESETS: dict[str, tuple[float, float, str]] = {
    "air": (1.40, 287.05, "Dry air, standard teaching default"),
    "air_hot": (1.33, 287.05, "Hot air / combustion-product-ish teaching value"),
    "lox_rp1": (1.22, 355.0, "LOX/RP-1 frozen engineering estimate near stoich"),
    "lox_rp1_fuel_rich": (1.25, 370.0, "LOX/RP-1 fuel-rich frozen estimate"),
    "lox_lh2": (1.20, 550.0, "LOX/LH2 frozen engineering estimate"),
    "lox_ch4": (1.20, 380.0, "LOX/CH4 (methalox) frozen engineering estimate"),
    "n2o_htpb": (1.25, 320.0, "N2O/HTPB hybrid-ish teaching estimate"),
    "monoprop_h2o2": (1.30, 380.0, "High-test peroxide decomposition estimate"),
}


def list_example_gases() -> list[str]:
    """Names accepted by :func:`example_gas`."""
    return sorted(_PRESETS)


def example_gas(name: str = "air") -> PerfectGas:
    """Return a named :class:`~rockmoc.thermo.perfect_gas.PerfectGas` preset.

    These are **engineering frozen estimates** for homework and demos — not
    flight thermo. Prefer RocketCEA / Cantera for propellant-accurate γ and R.

    Parameters
    ----------
    name:
        One of :func:`list_example_gases` (case-insensitive). Aliases:
        ``rp1`` → ``lox_rp1``, ``lh2`` → ``lox_lh2``, ``methalox`` → ``lox_ch4``.
    """
    key = name.strip().lower().replace("-", "_").replace(" ", "_")
    aliases = {
        "rp1": "lox_rp1",
        "loxrp1": "lox_rp1",
        "lh2": "lox_lh2",
        "hydrolox": "lox_lh2",
        "methalox": "lox_ch4",
        "methane": "lox_ch4",
        "default": "air",
    }
    key = aliases.get(key, key)
    if key not in _PRESETS:
        known = ", ".join(list_example_gases())
        raise ValueError(
            f"Unknown example gas {name!r}.\n"
            f"  Fix: choose one of: {known}\n"
            f"  Example: example_gas('lox_rp1')"
        )
    g, R, _note = _PRESETS[key]
    return PerfectGas(gamma=g, R=R)


def example_gas_info(name: str = "air") -> dict[str, float | str]:
    """Metadata for a preset (γ, R, note)."""
    gas = example_gas(name)
    key = name.strip().lower().replace("-", "_").replace(" ", "_")
    aliases = {"rp1": "lox_rp1", "lh2": "lox_lh2", "methalox": "lox_ch4", "default": "air"}
    key = aliases.get(key, key)
    _g, _R, note = _PRESETS[key]
    return {
        "name": key,
        "gamma": gas.gamma,
        "R": gas.R,
        "note": note,
        "disclaimer": "Frozen engineering estimate for teaching — not flight thermo.",
    }


def describe_example_gases() -> str:
    """Pretty table of all presets for terminals / notebooks."""
    lines = [
        f"{'name':18} {'γ':>6} {'R':>8}  note",
        "-" * 72,
    ]
    for name in list_example_gases():
        g, R, note = _PRESETS[name]
        lines.append(f"{name:18} {g:6.2f} {R:8.1f}  {note}")
    lines.append("")
    lines.append("Disclaimer: frozen teaching estimates only.")
    return "\n".join(lines)


# Convenience callables matching common student code style
def air() -> PerfectGas:
    return example_gas("air")


def lox_rp1() -> PerfectGas:
    return example_gas("lox_rp1")


def lox_lh2() -> PerfectGas:
    return example_gas("lox_lh2")


def lox_ch4() -> PerfectGas:
    return example_gas("lox_ch4")
