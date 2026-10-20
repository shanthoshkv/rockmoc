"""Production accuracy pipeline: fine MOC → strict wall → optional BL → field Cf.

``accuracy="production"`` enables dense mesh, Sauer throat metadata, expansion
table (multi-γ if provided), BL correction, and field-integrated performance.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal

from rockmoc.analysis.boundary_layer import BLCorrectionResult, correct_wall_bl
from rockmoc.analysis.throat import sauer_throat_ivl
from rockmoc.analysis.uq import mesh_certificate
from rockmoc.designers.base import DesignResult, DesignSpec
from rockmoc.designers.ideal import design_ideal
from rockmoc.designers.mln import design_mln
from rockmoc.designers.rao import design_rao
from rockmoc.designers.tic import design_tic
from rockmoc.exceptions import DesignError
from rockmoc.performance.field_cf import FieldCfResult, field_integrated_cf
from rockmoc.performance.integrated import off_design_altitude_sweep
from rockmoc.thermo.expansion_table import (
    ExpansionTable,
    build_multi_gamma_table,
    build_perfect_gas_table,
)
from rockmoc.thermo.perfect_gas import PerfectGas
from rockmoc.thermo.variable_gamma import VariableGammaGas


@dataclass
class ProductionResult:
    """Full production design package."""

    design: DesignResult
    field: FieldCfResult
    bl: BLCorrectionResult | None
    table: ExpansionTable | None
    throat_ivl_meta: dict
    certificate: dict
    altitude_sweep: list[dict]
    meta: dict = field(default_factory=dict)

    def summary(self) -> dict[str, Any]:
        s = self.design.summary()
        s["production"] = {
            "field_Cf": self.field.summary(),
            "bl": self.bl.summary() if self.bl else None,
            "certificate": self.certificate,
            "throat_ivl": self.throat_ivl_meta,
            "table_mode": self.table.mode if self.table else None,
            **self.meta,
        }
        return s

    def plot(self, **kwargs: Any) -> list[Any]:
        """Dark dashboard: contour, MOC, BL, altitude performance."""
        from rockmoc.plot.dashboard import display_production_result

        show = kwargs.pop("show", True)
        save_dir = kwargs.pop("save_dir", None)
        prefix = kwargs.pop("prefix", "production")
        return display_production_result(
            self, show=show, save_dir=save_dir, prefix=prefix
        )


def _resolve_gas_and_table(
    spec: DesignSpec,
    *,
    multi_gamma: tuple[float, float, float] | None = None,
    use_table: bool = True,
    M_max: float = 8.0,
) -> tuple[PerfectGas | VariableGammaGas, ExpansionTable | None]:
    gas = spec.gas
    if not use_table:
        return gas, None
    if multi_gamma is not None:
        gc, gt, ge = multi_gamma
        tab = build_multi_gamma_table(
            gamma_c=gc,
            gamma_t=gt,
            gamma_e=ge,
            R=gas.R,
            M_max=M_max,
            n=300,
        )
        return VariableGammaGas(tab), tab
    tab = build_perfect_gas_table(gas, M_max=M_max, n=300)
    return VariableGammaGas(tab), tab


def design_production(
    spec: DesignSpec,
    *,
    kind: Literal["mln", "ideal", "tic", "rao"] = "mln",
    apply_bl: bool = True,
    multi_gamma: tuple[float, float, float] | None = None,
    pc: float = 5.0e6,
    Tc: float = 3200.0,
    pa: float = 0.0,
    n_chars_min: int = 120,
    certificate: bool = True,
    altitude_points: list[float] | None = None,
    rao_length_fraction: float = 0.80,
    store_net: bool = True,
) -> ProductionResult:
    """Run the production accuracy pipeline.

    Parameters
    ----------
    multi_gamma:
        Optional (γ_c, γ_t, γ_e) for variable-γ table. If None, constant-γ
        table from ``spec.gas``.
    apply_bl:
        Apply integral BL displacement correction to the wall.
    """
    spec.validate()
    # Fine mesh
    n = max(int(spec.n_characteristics), n_chars_min)
    # Build thermo
    gas_m, table = _resolve_gas_and_table(
        spec, multi_gamma=multi_gamma, use_table=True
    )
    # Pass VariableGammaGas through MOC so γ varies with M along the nozzle
    Me_tmp = spec.resolved_exit_mach()
    trunc = spec.truncate_fraction
    if kind == "tic" and trunc is None:
        trunc = 0.80
    design_spec = DesignSpec(
        gas=gas_m,  # type: ignore[arg-type]  # VariableGammaGas or PerfectGas
        throat_radius=spec.throat_radius,
        exit_mach=Me_tmp,
        n_characteristics=n,
        throat_mach=max(spec.throat_mach, 1.05),
        throat_rc_over_rt=spec.throat_rc_over_rt,
        sharp_corner=spec.sharp_corner,
        delta=spec.delta,
        truncate_fraction=trunc,
    )

    Me = design_spec.resolved_exit_mach()

    # Throat IVL metadata (analytic Sauer — no CFD)
    ivl = sauer_throat_ivl(
        design_spec.gas,
        throat_radius=design_spec.throat_radius,
        rc_over_rt=design_spec.throat_rc_over_rt,
        n_points=41,
    )

    if kind == "mln":
        design = design_mln(design_spec, store_net=store_net)
    elif kind == "ideal":
        design = design_ideal(design_spec, store_net=store_net)
    elif kind == "tic":
        if design_spec.truncate_fraction is None:
            design_spec = DesignSpec(
                gas=design_spec.gas,
                throat_radius=design_spec.throat_radius,
                exit_mach=Me,
                n_characteristics=n,
                throat_mach=design_spec.throat_mach,
                throat_rc_over_rt=design_spec.throat_rc_over_rt,
                truncate_fraction=0.80,
                delta=design_spec.delta,
            )
        design = design_tic(design_spec, store_net=store_net)
    elif kind == "rao":
        design = design_rao(design_spec, length_fraction=rao_length_fraction)
    else:
        raise DesignError(f"unknown kind {kind}")

    contour = design.contour
    bl_res = None
    if apply_bl:
        bl_res = correct_wall_bl(
            contour,
            design_spec.gas,
            throat_radius=design_spec.throat_radius,
            pc=pc,
            Tc=Tc,
            turbulent=True,
            axisymmetric=(design_spec.delta == 1.0),
        )
        contour_for_perf = bl_res.corrected
    else:
        contour_for_perf = contour

    field = field_integrated_cf(
        contour_for_perf,
        design_spec.gas,
        throat_radius=design_spec.throat_radius,
        exit_mach_design=Me,
        pc=pc,
        Tc=Tc,
        pa=pa,
        net=design.net,
    )

    cert: dict = {"skipped": True}
    if certificate and kind in ("mln", "ideal", "tic"):
        try:
            cert = mesh_certificate(
                throat_radius=design_spec.throat_radius,
                exit_mach=Me,
                gamma=design_spec.gas.gamma,
                n_coarse=max(16, n // 3),
                n_fine=n,
                tol_rel=0.05,
            )
        except Exception as exc:
            cert = {"passed": False, "error": str(exc)}

    if altitude_points is None:
        altitude_points = [0.0, 10e3, 20e3, 40e3, 60e3]
    sweep = off_design_altitude_sweep(
        contour_for_perf,
        design_spec.gas,
        throat_radius=design_spec.throat_radius,
        exit_mach=Me,
        pc=pc,
        Tc=Tc,
        altitudes_m=altitude_points,
    )

    # Attach field metrics into design meta
    design.meta = {
        **design.meta,
        "production_pipeline": True,
        "field_cf": field.summary(),
        "bl_applied": apply_bl,
    }

    return ProductionResult(
        design=design,
        field=field,
        bl=bl_res,
        table=table,
        throat_ivl_meta=ivl.summary(),
        certificate=cert,
        altitude_sweep=sweep,
        meta={
            "kind": kind,
            "n_characteristics": n,
            "pc_Pa": pc,
            "Tc_K": Tc,
            "pa_Pa": pa,
            "multi_gamma": multi_gamma,
        },
    )
