"""Final thesis figure: strategy trade-off and geometric sensitivity.

The top row establishes why parallel sectioning concentrates more target in
fewer sections. The bottom row maps which surface-first serial section is the
earliest section with the maximum absolute target volume after small changes
in orientation and surface morphology.

All serial profiles are indexed in collection order. Section 1 is the first
0.1 mm slab whose outer boundary touches the surface-facing support plane;
subsequent sections move into the tissue.

Outputs are written to D:\\Vibratome_sections.
"""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import MultipleLocator

from app import (
    CUBE_COLOR,
    OVERLAP_COLOR,
    TARGET_COLOR,
    build_figure,
)
from geometry import (
    TARGET_CUBOID,
    TARGET_FOLLOW,
    TOP_FLAT,
    TOP_SPHERE_CONCAVE,
    TissueBlock,
    compute_volumes,
    consecutive_cuts,
    knife_normal,
    make_columns,
    target_z_bounds,
    tissue_box,
)


OUT = Path(r"D:\Vibratome_sections")
STEM = "parameter_sweep_geometric_model"
FIGURE_STEM = f"{STEM}_V6"

CUBE = (3.0, 3.0, 3.0)
TARGET = (3.0, 3.0, 0.2)
DEPTH = 0.1
KNIFE = 0.1
SECTION_SIZE = (6.0, 6.0, KNIFE)
GRID_N = 128

SECTION_BLUE = "#0072B2"
TARGET_ORANGE = "#E69F00"
FACE_COLOR = "#0072B2"
DIAGONAL_COLOR = "#7A5195"
FLAT_50_COLOR = "#CC3311"
FLAT_100_COLOR = "#0072B2"
BOWL_50_COLOR = "#D97706"
BOWL_100_COLOR = "#006D5B"
PARALLEL_BAND = "#009E73"
PERPENDICULAR_BAND = "#CC79A7"

ANGLE_VALUES = np.arange(0.0, 16.0, 1.0)
APEX_VALUES = np.arange(0.0, 1.0 + 1e-9, 0.1)


@dataclass
class SerialProfile:
    rx: float
    ry: float
    rz: float
    top_shape: str
    curve_offset: float
    target_mode: str
    section_numbers: np.ndarray
    midplanes: np.ndarray
    section_mm3: np.ndarray
    target_mm3: np.ndarray
    peak_sections: np.ndarray
    first_peak_section: int
    peak_target_mm3: float


@dataclass
class HeatmapData:
    name: str
    x_values: np.ndarray
    y_values: np.ndarray
    section_number: np.ndarray
    peak_target_mm3: np.ndarray
    tie_count: np.ndarray
    peak_sections: list[list[list[int]]]


def surface_support(block: TissueBlock, normal: np.ndarray) -> float:
    """Maximum projection of the actual tissue along the positive normal."""
    n = np.asarray(normal, dtype=float)
    if block.is_flat:
        return float(np.max(tissue_box(block.size).vertices() @ n))

    hx, hy = 0.5 * block.size[0], 0.5 * block.size[1]
    coords_x = np.linspace(-hx, hx, 201)
    coords_y = np.linspace(-hy, hy, 201)
    X, Y = np.meshgrid(coords_x, coords_y, indexing="xy")
    Z = block.z_surface(X, Y)
    top_projection = n[0] * X + n[1] * Y + n[2] * Z
    bottom_projection = n[0] * X + n[1] * Y + n[2] * block.z_bottom
    return float(max(np.max(top_projection), np.max(bottom_projection)))


def surface_first_profile(
    rx: float,
    ry: float,
    rz: float = 0.0,
    top_shape: str = TOP_FLAT,
    curve_offset: float = 0.0,
    target_mode: str = TARGET_CUBOID,
    columns=None,
) -> SerialProfile:
    """Return consecutive cuts ordered from first tissue contact inward."""
    block = TissueBlock(CUBE, top_shape, curve_offset, "y")
    normal = knife_normal(rx, ry, rz)
    support = surface_support(block, normal)
    pivot = (support - KNIFE / 2.0) * normal
    cols = columns
    if not block.is_flat and cols is None:
        cols = make_columns(block, TARGET, DEPTH, target_mode, n=GRID_N)

    cuts = consecutive_cuts(
        cube_size=CUBE,
        target_size=TARGET,
        depth_from_top=DEPTH,
        thickness=KNIFE,
        rx=rx,
        ry=ry,
        rz=rz,
        pivot=pivot,
        top_shape=top_shape,
        curve_offset=curve_offset,
        cylinder_axis="y",
        target_mode=target_mode,
        grid_n=GRID_N,
        _columns=cols,
    )
    cuts = sorted(
        (cut for cut in cuts if cut.tissue_volume > 1e-10),
        key=lambda cut: cut.midplane,
        reverse=True,
    )
    section_volume = np.array([cut.tissue_volume for cut in cuts], dtype=float)
    target_volume = np.array([cut.target_volume for cut in cuts], dtype=float)
    midplanes = np.array([cut.midplane for cut in cuts], dtype=float)
    peak = float(np.max(target_volume))
    peak_mask = np.isclose(target_volume, peak, rtol=1e-6, atol=1e-10)
    peak_sections = np.flatnonzero(peak_mask) + 1
    return SerialProfile(
        rx=float(rx),
        ry=float(ry),
        rz=float(rz),
        top_shape=top_shape,
        curve_offset=float(curve_offset),
        target_mode=target_mode,
        section_numbers=np.arange(1, len(cuts) + 1, dtype=int),
        midplanes=midplanes,
        section_mm3=section_volume,
        target_mm3=target_volume,
        peak_sections=peak_sections,
        first_peak_section=int(peak_sections[0]),
        peak_target_mm3=peak,
    )


def best_cut_curve_for_thickness(
    thickness: float,
    csv_path: Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Load or compute phase-optimised face and diagonal tilt curves."""
    tilts = np.arange(0.0, 91.0, 1.0)
    if csv_path.exists():
        data = np.loadtxt(csv_path, delimiter=",", skiprows=1)
        if data.shape == (len(tilts), 3) and np.allclose(data[:, 0], tilts):
            print(f"Loaded {thickness * 1000:.0f} µm best-cut curve:", csv_path, flush=True)
            return data[:, 0], data[:, 1], data[:, 2]

    print(
        f"Computing phase-optimised {thickness * 1000:.0f} µm best-cut curves...",
        flush=True,
    )

    def best_for(rx: float, rz: float) -> float:
        normal = knife_normal(rx, 0.0, rz)
        best = 0.0
        for phase in np.linspace(0.0, thickness, 8, endpoint=False):
            cuts = consecutive_cuts(
                cube_size=CUBE,
                target_size=TARGET,
                depth_from_top=DEPTH,
                thickness=thickness,
                rx=rx,
                ry=0.0,
                rz=rz,
                pivot=phase * normal,
            )
            best = max(best, max(cut.target_volume for cut in cuts))
        return best

    face = np.array([best_for(float(angle), 0.0) for angle in tilts])
    diagonal = np.array([best_for(float(angle), 45.0) for angle in tilts])
    np.savetxt(
        csv_path,
        np.column_stack([tilts, face, diagonal]),
        delimiter=",",
        header="tilt_deg,face_aligned_mm3,diagonal_mm3",
        comments="",
    )
    return tilts, face, diagonal


def compute_best_cut_curves() -> tuple[
    np.ndarray,
    dict[float, tuple[np.ndarray, np.ndarray]],
]:
    tilts_100, face_100, diagonal_100 = best_cut_curve_for_thickness(
        0.1,
        OUT / "max_target_vs_tilt.csv",
    )
    tilts_50, face_50, diagonal_50 = best_cut_curve_for_thickness(
        0.05,
        OUT / "max_target_vs_tilt_0p05.csv",
    )
    if not np.array_equal(tilts_50, tilts_100):
        raise AssertionError("Best-cut angle grids differ between thicknesses")
    return tilts_100, {
        0.05: (face_50, diagonal_50),
        0.10: (face_100, diagonal_100),
    }


def bowl_best_cut_curve_for_thickness(
    thickness: float,
    csv_path: Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Best-cut curve for the 0.7 mm spherical bowl and following target."""
    tilts = np.arange(0.0, 91.0, 1.0)
    if csv_path.exists():
        data = np.loadtxt(csv_path, delimiter=",", skiprows=1)
        if data.shape == (len(tilts), 3) and np.allclose(data[:, 0], tilts):
            print(
                f"Loaded {thickness * 1000:.0f} µm bowl best-cut curve:",
                csv_path,
                flush=True,
            )
            return data[:, 0], data[:, 1], data[:, 2]

    print(
        f"Computing {thickness * 1000:.0f} µm bowl best-cut curve...",
        flush=True,
    )
    block = TissueBlock(CUBE, TOP_SPHERE_CONCAVE, 0.7, "y")
    columns = make_columns(block, TARGET, DEPTH, TARGET_FOLLOW, n=64)

    def best_for(rx: float, rz: float) -> float:
        normal = knife_normal(rx, 0.0, rz)
        best = 0.0
        for phase in np.linspace(0.0, thickness, 6, endpoint=False):
            cuts = consecutive_cuts(
                cube_size=CUBE,
                target_size=TARGET,
                depth_from_top=DEPTH,
                thickness=thickness,
                rx=rx,
                ry=0.0,
                rz=rz,
                pivot=phase * normal,
                top_shape=TOP_SPHERE_CONCAVE,
                curve_offset=0.7,
                target_mode=TARGET_FOLLOW,
                grid_n=64,
                _columns=columns,
            )
            best = max(best, max(cut.target_volume for cut in cuts))
        return best

    face = np.array([best_for(float(angle), 0.0) for angle in tilts])
    diagonal = np.array([best_for(float(angle), 45.0) for angle in tilts])
    np.savetxt(
        csv_path,
        np.column_stack([tilts, face, diagonal]),
        delimiter=",",
        header="tilt_deg,face_aligned_mm3,diagonal_mm3",
        comments="",
    )
    return tilts, face, diagonal


def compute_bowl_best_cut_curves() -> dict[float, tuple[np.ndarray, np.ndarray]]:
    tilts_50, face_50, diagonal_50 = bowl_best_cut_curve_for_thickness(
        0.05,
        OUT / "max_target_vs_tilt_sphere_concave_follow_h0.7_t0.05.csv",
    )
    tilts_100, face_100, diagonal_100 = bowl_best_cut_curve_for_thickness(
        0.10,
        OUT / "max_target_vs_tilt_sphere_concave_follow_h0.7.csv",
    )
    if not np.array_equal(tilts_50, tilts_100):
        raise AssertionError("Bowl best-cut angle grids differ between thicknesses")
    return {
        0.05: (face_50, diagonal_50),
        0.10: (face_100, diagonal_100),
    }


def compute_flat_xy() -> HeatmapData:
    print("Computing flat cube–cuboid X/Y rotation map...", flush=True)
    n = len(ANGLE_VALUES)
    sections = np.empty((n, n), dtype=int)
    peaks = np.empty((n, n), dtype=float)
    ties = np.empty((n, n), dtype=int)
    all_peaks: list[list[list[int]]] = [[[] for _ in range(n)] for _ in range(n)]
    for iy, ry in enumerate(ANGLE_VALUES):
        for ix, rx in enumerate(ANGLE_VALUES):
            profile = surface_first_profile(
                rx,
                ry,
                top_shape=TOP_FLAT,
                target_mode=TARGET_CUBOID,
            )
            sections[iy, ix] = profile.first_peak_section
            peaks[iy, ix] = profile.peak_target_mm3
            ties[iy, ix] = len(profile.peak_sections)
            all_peaks[iy][ix] = profile.peak_sections.tolist()
        print(f"  Y rotation {ry:.0f}° complete", flush=True)
    return HeatmapData(
        "flat_xy", ANGLE_VALUES.copy(), ANGLE_VALUES.copy(), sections, peaks, ties, all_peaks
    )


def compute_concavity_x() -> HeatmapData:
    print("Computing bowl-depth/X-rotation map...", flush=True)
    nx = len(ANGLE_VALUES)
    ny = len(APEX_VALUES)
    sections = np.empty((ny, nx), dtype=int)
    peaks = np.empty((ny, nx), dtype=float)
    ties = np.empty((ny, nx), dtype=int)
    all_peaks: list[list[list[int]]] = [[[] for _ in range(nx)] for _ in range(ny)]
    for iy, apex in enumerate(APEX_VALUES):
        shape = TOP_FLAT if np.isclose(apex, 0.0) else TOP_SPHERE_CONCAVE
        block = TissueBlock(CUBE, shape, float(apex), "y")
        cols = (
            None
            if block.is_flat
            else make_columns(block, TARGET, DEPTH, TARGET_FOLLOW, n=GRID_N)
        )
        for ix, rx in enumerate(ANGLE_VALUES):
            profile = surface_first_profile(
                rx,
                0.0,
                top_shape=shape,
                curve_offset=float(apex),
                target_mode=TARGET_FOLLOW,
                columns=cols,
            )
            sections[iy, ix] = profile.first_peak_section
            peaks[iy, ix] = profile.peak_target_mm3
            ties[iy, ix] = len(profile.peak_sections)
            all_peaks[iy][ix] = profile.peak_sections.tolist()
        print(f"  Bowl depth {apex:.1f} mm complete", flush=True)
    return HeatmapData(
        "concavity_x", ANGLE_VALUES.copy(), APEX_VALUES.copy(), sections, peaks, ties, all_peaks
    )


def compute_concave_xy(apex: float = 0.7) -> HeatmapData:
    print(f"Computing X/Y rotation map for {apex:.1f} mm bowl...", flush=True)
    n = len(ANGLE_VALUES)
    sections = np.empty((n, n), dtype=int)
    peaks = np.empty((n, n), dtype=float)
    ties = np.empty((n, n), dtype=int)
    all_peaks: list[list[list[int]]] = [[[] for _ in range(n)] for _ in range(n)]
    block = TissueBlock(CUBE, TOP_SPHERE_CONCAVE, apex, "y")
    cols = make_columns(block, TARGET, DEPTH, TARGET_FOLLOW, n=GRID_N)
    for iy, ry in enumerate(ANGLE_VALUES):
        for ix, rx in enumerate(ANGLE_VALUES):
            profile = surface_first_profile(
                rx,
                ry,
                top_shape=TOP_SPHERE_CONCAVE,
                curve_offset=apex,
                target_mode=TARGET_FOLLOW,
                columns=cols,
            )
            sections[iy, ix] = profile.first_peak_section
            peaks[iy, ix] = profile.peak_target_mm3
            ties[iy, ix] = len(profile.peak_sections)
            all_peaks[iy][ix] = profile.peak_sections.tolist()
        print(f"  Y rotation {ry:.0f}° complete", flush=True)
    return HeatmapData(
        "concave07_xy",
        ANGLE_VALUES.copy(),
        ANGLE_VALUES.copy(),
        sections,
        peaks,
        ties,
        all_peaks,
    )


def write_curve_data(
    tilts: np.ndarray,
    curves: dict[float, tuple[np.ndarray, np.ndarray]],
    bowl_curves: dict[float, tuple[np.ndarray, np.ndarray]],
) -> None:
    path = OUT / f"{STEM}_best_cut.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "tilt_deg",
                "face_50um_mm3",
                "diagonal_50um_mm3",
                "face_100um_mm3",
                "diagonal_100um_mm3",
                "bowl_face_50um_mm3",
                "bowl_diagonal_50um_mm3",
                "bowl_face_100um_mm3",
                "bowl_diagonal_100um_mm3",
            ]
        )
        face_50, diagonal_50 = curves[0.05]
        face_100, diagonal_100 = curves[0.10]
        bowl_face_50, bowl_diagonal_50 = bowl_curves[0.05]
        bowl_face_100, bowl_diagonal_100 = bowl_curves[0.10]
        for row in zip(
            tilts,
            face_50,
            diagonal_50,
            face_100,
            diagonal_100,
            bowl_face_50,
            bowl_diagonal_50,
            bowl_face_100,
            bowl_diagonal_100,
        ):
            writer.writerow(row)
    print("wrote", path, flush=True)


def write_profile_data(profiles: dict[str, SerialProfile]) -> None:
    path = OUT / f"{STEM}_serial_profiles.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "orientation",
                "section_number_surface_inward",
                "absolute_section_volume_mm3",
                "absolute_target_volume_mm3",
                "is_highest_yielding",
            ]
        )
        for name, profile in profiles.items():
            peaks = set(profile.peak_sections.tolist())
            for section, section_volume, target_volume in zip(
                profile.section_numbers,
                profile.section_mm3,
                profile.target_mm3,
            ):
                writer.writerow(
                    [
                        name,
                        int(section),
                        f"{section_volume:.10g}",
                        f"{target_volume:.10g}",
                        int(section) in peaks,
                    ]
                )
    print("wrote", path, flush=True)


def write_heatmap_data(data: HeatmapData) -> None:
    path = OUT / f"{STEM}_{data.name}.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "x_value",
                "y_value",
                "earliest_highest_yielding_section",
                "all_highest_yielding_sections",
                "number_of_tied_maxima",
                "maximum_absolute_target_volume_mm3",
            ]
        )
        for iy, y_value in enumerate(data.y_values):
            for ix, x_value in enumerate(data.x_values):
                writer.writerow(
                    [
                        f"{x_value:.6g}",
                        f"{y_value:.6g}",
                        int(data.section_number[iy, ix]),
                        ";".join(str(x) for x in data.peak_sections[iy][ix]),
                        int(data.tie_count[iy, ix]),
                        f"{data.peak_target_mm3[iy, ix]:.10g}",
                    ]
                )
    print("wrote", path, flush=True)


def write_metadata() -> None:
    path = OUT / f"{STEM}_metadata.json"
    metadata = {
        "cube_mm": CUBE,
        "target_mm": TARGET,
        "target_depth_mm": DEPTH,
        "knife_thickness_mm": KNIFE,
        "best_cut_curve_section_thicknesses_mm": [0.05, 0.10],
        "best_cut_curve_orientations": [
            "face-aligned",
            "diagonal (45 degrees in Z)",
        ],
        "best_cut_curve_models": [
            "flat cube-cuboid",
            "spherical concave 0.7 mm bowl with surface-following target",
        ],
        "serial_profile_scenarios": [
            "flat 0 degrees",
            "flat 15 degrees in X",
            "flat 90 degrees in X",
            "spherical concave 0.7 mm bowl, 0 degrees",
            "spherical concave 0.7 mm bowl, 15 degrees in X",
            "spherical concave 0.7 mm bowl, 90 degrees in X",
        ],
        "serial_index_definition": (
            "Section 1 is the first knife slab whose outer boundary touches "
            "the tissue support plane on the surface-facing side. Numbering "
            "then proceeds inward along the negative knife normal."
        ),
        "heatmap_value": (
            "Earliest section among all sections tied for maximum absolute "
            "target volume. Cells with more than one tied maximum are marked."
        ),
        "flat_xy": {
            "target_mode": TARGET_CUBOID,
            "rx_deg": ANGLE_VALUES.tolist(),
            "ry_deg": ANGLE_VALUES.tolist(),
        },
        "concavity_x": {
            "target_mode": TARGET_FOLLOW,
            "surface": "spherical concave",
            "bowl_depth_mm": APEX_VALUES.tolist(),
            "rx_deg": ANGLE_VALUES.tolist(),
            "ry_deg": 0.0,
        },
        "concave07_xy": {
            "target_mode": TARGET_FOLLOW,
            "surface": "spherical concave",
            "bowl_depth_mm": 0.7,
            "rx_deg": ANGLE_VALUES.tolist(),
            "ry_deg": ANGLE_VALUES.tolist(),
        },
        "curved_xy_grid_n": GRID_N,
    }
    path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print("wrote", path, flush=True)


def choose_scene_section(profile: SerialProfile, middle_target: bool = False) -> int:
    if middle_target:
        positive = np.flatnonzero(profile.target_mm3 > 1e-10)
        return int(positive[len(positive) // 2])
    return int(profile.first_peak_section - 1)


def style_scene(fig) -> None:
    hidden = {
        "Display cuboid (not used in volumes)",
        "section-edges",
        "original-cube",
        "Rotation pivot",
        "Tissue bottom",
    }
    for trace in fig.data:
        name = trace.name or ""
        if name in hidden:
            trace.visible = False
            continue
        trace.showlegend = False
        if name == "Tissue in this cut":
            trace.color = SECTION_BLUE
            trace.opacity = 0.62
        elif name == "Target in this cut":
            trace.color = TARGET_ORANGE
            trace.opacity = 0.96
        elif name == "Target region":
            trace.color = TARGET_COLOR
        elif name == "Tissue cube":
            trace.opacity = 0.09
    fig.update_layout(
        title=None,
        showlegend=False,
        margin=dict(l=0, r=0, t=0, b=0),
        paper_bgcolor="white",
        plot_bgcolor="white",
    )
    axis_style = dict(
        range=[-1.85, 1.85],
        showticklabels=False,
        ticks="",
        showgrid=False,
        zeroline=False,
        showline=True,
        linecolor="#334155",
        linewidth=4,
        showbackground=False,
        showspikes=False,
    )
    fig.update_scenes(
        xaxis=dict(
            **axis_style,
            title=dict(text=""),
        ),
        yaxis=dict(
            **axis_style,
            title=dict(text=""),
        ),
        zaxis=dict(
            **axis_style,
            title=dict(text=""),
        ),
        camera=dict(eye=dict(x=1.7, y=1.45, z=1.18)),
        aspectmode="cube",
    )
    for text, x_pos, y_pos in (
        ("Z", 0.15, 0.51),
        ("Y", 0.31, 0.17),
        ("X", 0.74, 0.17),
    ):
        fig.add_annotation(
            x=x_pos,
            y=y_pos,
            xref="paper",
            yref="paper",
            text=f"<b>{text}</b>",
            showarrow=False,
            font=dict(size=34, color="#111827"),
        )


def export_scene(
    profile: SerialProfile,
    path: Path,
    middle_target: bool = False,
) -> None:
    index = choose_scene_section(profile, middle_target=middle_target)
    normal = knife_normal(profile.rx, profile.ry, profile.rz)
    pivot = profile.midplanes[index] * normal
    (
        _volumes,
        tissue,
        target,
        section,
        overlap_points,
        tissue_cut_points,
        _stack,
        draw,
    ) = compute_volumes(
        cube_size=CUBE,
        target_size=TARGET,
        depth_from_top=DEPTH,
        section_size=SECTION_SIZE,
        rx=profile.rx,
        ry=profile.ry,
        rz=profile.rz,
        pivot=pivot,
        top_shape=profile.top_shape,
        curve_offset=profile.curve_offset,
        cylinder_axis="y",
        target_mode=profile.target_mode,
        grid_n=GRID_N,
        vis_n=44,
    )
    fig = build_figure(
        tissue,
        target,
        section,
        overlap_points,
        tissue_cut_points,
        pivot,
        draw,
        CUBE,
    )
    style_scene(fig)
    fig.write_image(str(path), width=760, height=650, scale=2)
    print("wrote", path, flush=True)


def _ink_bbox(image: np.ndarray, pad: int = 8) -> tuple[int, int, int, int]:
    rgb = image[:, :, :3]
    threshold = 0.97 if np.issubdtype(image.dtype, np.floating) else 248
    ink = np.any(rgb < threshold, axis=2)
    rows = np.where(ink.any(axis=1))[0]
    cols = np.where(ink.any(axis=0))[0]
    if rows.size == 0 or cols.size == 0:
        return 0, image.shape[0], 0, image.shape[1]
    return (
        max(0, int(rows[0]) - pad),
        min(image.shape[0], int(rows[-1]) + pad + 1),
        max(0, int(cols[0]) - pad),
        min(image.shape[1], int(cols[-1]) + pad + 1),
    )


def crop_shared(images: list[np.ndarray]) -> list[np.ndarray]:
    boxes = [_ink_bbox(image) for image in images]
    r0 = min(box[0] for box in boxes)
    r1 = max(box[1] for box in boxes)
    c0 = min(box[2] for box in boxes)
    c1 = max(box[3] for box in boxes)
    return [image[r0:r1, c0:c1] for image in images]


def draw_morphology_cross_section(
    ax: plt.Axes,
    apex: float,
    profile: SerialProfile,
    tilted_profile: SerialProfile,
    show_axis_reference: bool = False,
) -> None:
    shape = TOP_FLAT if np.isclose(apex, 0.0) else TOP_SPHERE_CONCAVE
    block = TissueBlock(CUBE, shape, apex, "y")
    qmax = CUBE[1] / 2.0
    q = np.linspace(-qmax, qmax, 500)
    x = np.zeros_like(q)
    y = q
    top = block.z_surface(x, y)
    target_lo, target_hi = target_z_bounds(
        x, y, block, TARGET, DEPTH, TARGET_FOLLOW
    )
    cut_mid = profile.midplanes[profile.first_peak_section - 1]
    cut_lo = np.full_like(q, cut_mid - KNIFE / 2.0)
    cut_hi = np.full_like(q, cut_mid + KNIFE / 2.0)
    tissue_cut_lo = np.maximum(block.z_bottom, cut_lo)
    tissue_cut_hi = np.minimum(top, cut_hi)
    overlap_lo = np.maximum(target_lo, cut_lo)
    overlap_hi = np.minimum(target_hi, cut_hi)

    ax.fill_between(q, block.z_bottom, top, color=CUBE_COLOR, alpha=0.16, lw=0)
    ax.fill_between(
        q,
        target_lo,
        target_hi,
        where=np.isfinite(target_lo) & np.isfinite(target_hi),
        color=TARGET_COLOR,
        alpha=0.42,
        lw=0,
    )
    ax.fill_between(
        q,
        tissue_cut_lo,
        tissue_cut_hi,
        where=tissue_cut_hi > tissue_cut_lo,
        color=SECTION_BLUE,
        alpha=0.62,
        lw=0,
    )
    ax.fill_between(
        q,
        overlap_lo,
        overlap_hi,
        where=np.isfinite(overlap_lo)
        & np.isfinite(overlap_hi)
        & (overlap_hi > overlap_lo),
        color=TARGET_ORANGE,
        alpha=0.96,
        lw=0,
    )
    tilted_normal = knife_normal(
        tilted_profile.rx,
        tilted_profile.ry,
        tilted_profile.rz,
    )
    tilted_mid = tilted_profile.midplanes[tilted_profile.first_peak_section - 1]
    in_plane_projection = tilted_normal[0] * x + tilted_normal[1] * y
    tilted_lo = (
        tilted_mid - KNIFE / 2.0 - in_plane_projection
    ) / tilted_normal[2]
    tilted_hi = (
        tilted_mid + KNIFE / 2.0 - in_plane_projection
    ) / tilted_normal[2]
    ax.plot(q, tilted_lo, color="#1D4ED8", lw=1.0, ls="--")
    ax.plot(q, tilted_hi, color="#1D4ED8", lw=1.0, ls="--")
    ax.plot(q, top, color="#64748b", lw=0.9)
    ax.set_xlim(-qmax, qmax)
    ax.set_ylim(0.15, 1.65)
    ax.set_box_aspect(0.50)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(f"{apex:.1f} mm", fontsize=8, pad=1)
    if show_axis_reference:
        origin = (0.10, 0.16)
        ax.annotate(
            "",
            xy=(0.31, origin[1]),
            xytext=origin,
            xycoords="axes fraction",
            arrowprops=dict(arrowstyle="-|>", color="#111827", lw=0.9),
        )
        ax.annotate(
            "",
            xy=(origin[0], 0.42),
            xytext=origin,
            xycoords="axes fraction",
            arrowprops=dict(arrowstyle="-|>", color="#111827", lw=0.9),
        )
        ax.text(
            0.33,
            0.13,
            r"$y$",
            transform=ax.transAxes,
            fontsize=8,
            ha="left",
            va="center",
        )
        ax.text(
            0.08,
            0.44,
            r"$z$",
            transform=ax.transAxes,
            fontsize=8,
            ha="center",
            va="bottom",
        )
    for spine in ax.spines.values():
        spine.set_visible(False)


def generate_visuals(
    flat_parallel: SerialProfile,
    flat_perpendicular: SerialProfile,
    bowl_perpendicular: SerialProfile,
    flat_xy: HeatmapData,
    concavity_x: HeatmapData,
    concave_xy: HeatmapData,
) -> dict[str, Path]:
    del flat_xy, concavity_x, concave_xy  # data objects document matching scenarios
    print("Exporting model visualizations...", flush=True)
    profiles = {
        "parallel": flat_parallel,
        "perpendicular": flat_perpendicular,
        "bowl_perpendicular": bowl_perpendicular,
        "flat_00": surface_first_profile(0.0, 0.0),
        "flat_1010": surface_first_profile(10.0, 10.0),
        "bowl_00": surface_first_profile(
            0.0,
            0.0,
            top_shape=TOP_SPHERE_CONCAVE,
            curve_offset=0.7,
            target_mode=TARGET_FOLLOW,
        ),
        "bowl_1010": surface_first_profile(
            10.0,
            10.0,
            top_shape=TOP_SPHERE_CONCAVE,
            curve_offset=0.7,
            target_mode=TARGET_FOLLOW,
        ),
    }
    paths: dict[str, Path] = {}
    for name, profile in profiles.items():
        path = OUT / f"_{STEM}_{name}.png"
        export_scene(
            profile,
            path,
            middle_target=(name == "perpendicular"),
        )
        paths[name] = path
    return paths


def configure_matplotlib() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.size": 8.5,
            "axes.linewidth": 0.8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def clean_axis(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def add_panel_letter(ax: plt.Axes, letter: str, x: float = -0.10, y: float = 1.09) -> None:
    ax.text(
        x,
        y,
        letter,
        transform=ax.transAxes,
        fontsize=13,
        fontweight="bold",
        va="top",
        ha="left",
    )


def show_images(
    axes: list[plt.Axes],
    paths: list[Path],
    titles: list[str],
    title_fontsize: float = 8.5,
) -> None:
    images = crop_shared([plt.imread(path) for path in paths])
    for ax, image, title in zip(axes, images, titles):
        ax.imshow(image)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(title, fontsize=title_fontsize, pad=2)
        for spine in ax.spines.values():
            spine.set_visible(False)


def heatmap_palette(
    heatmaps: list[HeatmapData],
) -> tuple[ListedColormap, BoundaryNorm, int, int]:
    minimum = min(int(np.min(data.section_number)) for data in heatmaps)
    maximum = max(int(np.max(data.section_number)) for data in heatmaps)
    colors = plt.get_cmap("viridis")(np.linspace(0.06, 0.94, maximum - minimum + 1))
    cmap = ListedColormap(colors)
    boundaries = np.arange(minimum - 0.5, maximum + 1.5, 1.0)
    norm = BoundaryNorm(boundaries, cmap.N)
    return cmap, norm, minimum, maximum


def plot_heatmap(
    ax: plt.Axes,
    data: HeatmapData,
    cmap: ListedColormap,
    norm: BoundaryNorm,
    x_label: str,
    y_label: str,
    equal: bool,
) -> object:
    x_edges = np.concatenate(
        ([data.x_values[0] - 0.5], 0.5 * (data.x_values[:-1] + data.x_values[1:]), [data.x_values[-1] + 0.5])
    )
    y_step = float(data.y_values[1] - data.y_values[0])
    y_edges = np.concatenate(
        (
            [data.y_values[0] - y_step / 2.0],
            0.5 * (data.y_values[:-1] + data.y_values[1:]),
            [data.y_values[-1] + y_step / 2.0],
        )
    )
    mesh = ax.pcolormesh(
        x_edges,
        y_edges,
        data.section_number,
        cmap=cmap,
        norm=norm,
        shading="flat",
        rasterized=True,
    )

    tied_y, tied_x = np.where(data.tie_count > 1)
    if tied_x.size:
        xs = data.x_values[tied_x]
        ys = data.y_values[tied_y]
        ax.scatter(xs, ys, marker="x", s=14, color="white", linewidth=1.1, zorder=4)
        ax.scatter(xs, ys, marker="x", s=14, color="#111827", linewidth=0.35, zorder=5)

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    if np.allclose(np.diff(data.x_values), 1.0):
        ax.set_xticks(data.x_values)
        ax.set_xticklabels([f"{int(value)}" for value in data.x_values])
    if np.allclose(np.diff(data.y_values), 1.0):
        ax.set_yticks(data.y_values)
        ax.set_yticklabels([f"{int(value)}" for value in data.y_values])
    else:
        ax.set_yticks(data.y_values)
        ax.set_yticklabels([f"{value:.1f}" for value in data.y_values])
    ax.set_xticks(x_edges[1:-1], minor=True)
    ax.set_yticks(y_edges[1:-1], minor=True)
    ax.grid(which="minor", color="white", lw=0.25, alpha=0.22)
    ax.tick_params(which="minor", length=0)
    ax.tick_params(which="major", labelsize=6.0, length=2.0, pad=1.5)
    ax.set_aspect("equal" if equal else "auto")
    return mesh


def plot_serial_bars(
    ax: plt.Axes,
    profile: SerialProfile,
    title: str,
    x_max: int,
    show_xlabel: bool,
) -> None:
    x = profile.section_numbers
    ax.bar(
        x,
        profile.section_mm3,
        width=0.84,
        color=SECTION_BLUE,
        alpha=0.72,
        edgecolor="white",
        linewidth=0.25,
        label="Absolute section volume",
        zorder=1,
    )
    ax.bar(
        x,
        profile.target_mm3,
        width=0.48,
        color=TARGET_ORANGE,
        edgecolor="white",
        linewidth=0.25,
        label="Absolute target volume",
        zorder=2,
    )
    target_positive = int(np.sum(profile.target_mm3 > 1e-9))
    ax.set_title(
        f"{title}\n{target_positive} target-positive sections",
        fontsize=7.0,
        pad=2,
    )
    ax.set_xlim(0.4, x_max + 0.6)
    ax.set_ylim(0, 1.02)
    if show_xlabel:
        ax.set_xlabel("Collected section number", fontsize=7)
    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.tick_params(labelsize=6.5)
    clean_axis(ax)


def make_figure(
    tilts: np.ndarray,
    face: np.ndarray,
    diagonal: np.ndarray,
    strategy_profiles: dict[str, SerialProfile],
    flat_xy: HeatmapData,
    concavity_x: HeatmapData,
    concave_xy: HeatmapData,
    visual_paths: dict[str, Path],
) -> tuple[Path, Path]:
    configure_matplotlib()
    heatmaps = [flat_xy, concavity_x, concave_xy]
    cmap, norm, minimum_section, maximum_section = heatmap_palette(heatmaps)

    fig = plt.figure(figsize=(16.2, 10.8))
    outer = fig.add_gridspec(
        2,
        1,
        height_ratios=[0.92, 1.08],
        hspace=0.30,
        left=0.055,
        right=0.955,
        top=0.900,
        bottom=0.075,
    )
    top = outer[0].subgridspec(
        1, 3, width_ratios=[2.2, 3.5, 6.3], wspace=0.27
    )
    bottom = outer[1].subgridspec(1, 3, wspace=0.25)

    # a) Geometry for parallel and perpendicular strategies.
    a_grid = top[0, 0].subgridspec(1, 2, wspace=0.03)
    ax_a = [fig.add_subplot(a_grid[0, i]) for i in range(2)]
    show_images(
        ax_a,
        [visual_paths["parallel"], visual_paths["perpendicular"]],
        ["Parallel (0°)", "Perpendicular (90°)"],
    )
    add_panel_letter(ax_a[0], "a", x=-0.18, y=1.28)
    ax_a[0].text(
        1.0,
        1.28,
        "Sectioning geometry",
        transform=ax_a[0].transAxes,
        ha="center",
        va="top",
        fontsize=10,
        fontweight="semibold",
    )

    # b) Original best-cut-vs-tilt graph.
    ax_b = fig.add_subplot(top[0, 1])
    ax_b.axvspan(0, 15, color=PARALLEL_BAND, alpha=0.12, lw=0)
    ax_b.axvspan(75, 90, color=PERPENDICULAR_BAND, alpha=0.12, lw=0)
    ax_b.plot(tilts, face, color=FACE_COLOR, lw=2.1, label="Face-aligned")
    ax_b.plot(
        tilts,
        diagonal,
        color=DIAGONAL_COLOR,
        lw=2.1,
        label="Diagonal (45° in Z)",
    )
    ax_b.set_xlim(0, 90)
    ax_b.set_ylim(0, 1.0)
    ax_b.set_xlabel("Tilt from the top face (°)")
    ax_b.set_ylabel("Target in best serial 0.1 mm cut (mm³)")
    ax_b.legend(frameon=False, fontsize=8, loc="upper right")
    ax_b.set_title("Maximum target volume per section", fontsize=10, fontweight="semibold")
    clean_axis(ax_b)
    add_panel_letter(ax_b, "b", x=-0.15, y=1.10)

    # c) Surface-indexed absolute section and target volumes.
    c_grid = top[0, 2].subgridspec(1, 2, wspace=0.16)
    ax_c1 = fig.add_subplot(c_grid[0, 0])
    ax_c2 = fig.add_subplot(c_grid[0, 1], sharey=ax_c1)
    plot_serial_bars(ax_c1, strategy_profiles["parallel"], "Parallel")
    plot_serial_bars(ax_c2, strategy_profiles["perpendicular"], "Perpendicular")
    ax_c1.set_ylabel("Absolute volume (mm³)")
    ax_c2.tick_params(axis="y", labelleft=False, left=False)
    ax_c1.legend(
        handles=[
            Patch(facecolor=SECTION_BLUE, alpha=0.72, label="Total section"),
            Patch(facecolor=TARGET_ORANGE, label="Target in section"),
        ],
        frameon=False,
        fontsize=8,
        ncol=2,
        loc="upper center",
        bbox_to_anchor=(1.03, 1.10),
    )
    add_panel_letter(ax_c1, "c", x=-0.13, y=1.20)
    ax_c1.text(
        1.03,
        1.20,
        "Absolute volumes across consecutive sections",
        transform=ax_c1.transAxes,
        ha="center",
        va="top",
        fontsize=10,
        fontweight="semibold",
    )

    # d) Flat cube-cuboid X/Y rotations.
    d_grid = bottom[0, 0].subgridspec(2, 1, height_ratios=[0.75, 2.35], hspace=0.13)
    d_visual = d_grid[0].subgridspec(1, 2, wspace=0.03)
    ax_dv = [fig.add_subplot(d_visual[0, i]) for i in range(2)]
    show_images(
        ax_dv,
        [visual_paths["flat_00"], visual_paths["flat_1010"]],
        ["X=0°, Y=0°", "X=10°, Y=10°"],
    )
    add_panel_letter(ax_dv[0], "d", x=-0.14, y=1.35)
    ax_dv[0].text(
        1.0,
        1.35,
        "Flat cube–cuboid",
        transform=ax_dv[0].transAxes,
        ha="center",
        va="top",
        fontsize=10,
        fontweight="semibold",
    )
    ax_d = fig.add_subplot(d_grid[1])
    mesh = plot_heatmap(
        ax_d,
        flat_xy,
        cmap,
        norm,
        "X rotation (°)",
        "Y rotation (°)",
        equal=True,
    )

    # e) Bowl depth and X rotation.
    e_grid = bottom[0, 1].subgridspec(2, 1, height_ratios=[0.75, 2.35], hspace=0.13)
    e_visual = e_grid[0].subgridspec(1, 3, wspace=0.04)
    ax_ev = [fig.add_subplot(e_visual[0, i]) for i in range(3)]
    for ax, apex in zip(ax_ev, (0.0, 0.5, 1.0)):
        iy = int(round(apex * 10))
        profile = surface_first_profile(
            0.0,
            0.0,
            top_shape=TOP_FLAT if np.isclose(apex, 0.0) else TOP_SPHERE_CONCAVE,
            curve_offset=apex,
            target_mode=TARGET_FOLLOW,
        )
        del iy
        draw_morphology_cross_section(ax, apex, profile)
    add_panel_letter(ax_ev[0], "e", x=-0.18, y=1.42)
    ax_ev[1].text(
        0.5,
        1.42,
        "Spherical concavity × X rotation",
        transform=ax_ev[1].transAxes,
        ha="center",
        va="top",
        fontsize=10,
        fontweight="semibold",
    )
    ax_e = fig.add_subplot(e_grid[1])
    plot_heatmap(
        ax_e,
        concavity_x,
        cmap,
        norm,
        "X rotation (°)",
        "Bowl depth (mm)",
        equal=False,
    )

    # f) X/Y rotations at 0.7 mm bowl depth.
    f_grid = bottom[0, 2].subgridspec(2, 1, height_ratios=[0.75, 2.35], hspace=0.13)
    f_visual = f_grid[0].subgridspec(1, 2, wspace=0.03)
    ax_fv = [fig.add_subplot(f_visual[0, i]) for i in range(2)]
    show_images(
        ax_fv,
        [visual_paths["bowl_00"], visual_paths["bowl_1010"]],
        ["X=0°, Y=0°", "X=10°, Y=10°"],
    )
    ax_fv[0].text(
        1.0,
        1.35,
        "0.7 mm bowl · surface-following target",
        transform=ax_fv[0].transAxes,
        ha="center",
        va="top",
        fontsize=10,
        fontweight="semibold",
    )
    ax_f = fig.add_subplot(f_grid[1])
    ax_fv[0].text(
        0.02,
        0.95,
        "f",
        transform=ax_fv[0].transAxes,
        fontsize=13,
        fontweight="bold",
        va="top",
        ha="left",
    )
    plot_heatmap(
        ax_f,
        concave_xy,
        cmap,
        norm,
        "X rotation (°)",
        "Y rotation (°)",
        equal=True,
    )

    cax = fig.add_axes([0.968, 0.105, 0.011, 0.34])
    cbar = fig.colorbar(mesh, cax=cax)
    cbar.set_ticks(np.arange(minimum_section, maximum_section + 1))
    cbar.set_label("Earliest highest-yielding section number", fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    fig.legend(
        handles=[
            Patch(facecolor=CUBE_COLOR, alpha=0.20, edgecolor="none", label="Tissue"),
            Patch(facecolor=SECTION_BLUE, alpha=0.70, edgecolor="none", label="Section"),
            Patch(facecolor=TARGET_COLOR, alpha=0.42, edgecolor="none", label="Target"),
            Patch(facecolor=TARGET_ORANGE, edgecolor="none", label="Target in section"),
            Line2D(
                [0],
                [0],
                marker="x",
                color="#111827",
                lw=0,
                markersize=5,
                label="Multiple sections tied for maximum",
            ),
        ],
        loc="upper center",
        bbox_to_anchor=(0.50, 0.993),
        ncol=5,
        frameon=False,
        fontsize=8,
        columnspacing=1.25,
        handlelength=1.2,
    )

    png = OUT / f"{STEM}.png"
    pdf = OUT / f"{STEM}.pdf"
    temporary_png = OUT / f"_{STEM}_temporary.png"
    temporary_pdf = OUT / f"_{STEM}_temporary.pdf"
    fig.savefig(temporary_png, dpi=300, facecolor="white")
    fig.savefig(temporary_pdf, facecolor="white")
    plt.close(fig)
    if png.exists():
        png.unlink()
    if pdf.exists():
        pdf.unlink()
    temporary_png.replace(png)
    temporary_pdf.replace(pdf)
    print("wrote", png, flush=True)
    print("wrote", pdf, flush=True)
    return png, pdf


def make_figure_v2(
    tilts: np.ndarray,
    curves: dict[float, tuple[np.ndarray, np.ndarray]],
    bowl_curves: dict[float, tuple[np.ndarray, np.ndarray]],
    strategy_profiles: dict[str, SerialProfile],
    flat_xy: HeatmapData,
    concavity_x: HeatmapData,
    concave_xy: HeatmapData,
    visual_paths: dict[str, Path],
) -> tuple[Path, Path]:
    """Compose the revised six-panel figure with aligned panel headers."""
    configure_matplotlib()
    heatmaps = [flat_xy, concavity_x, concave_xy]
    cmap, norm, minimum_section, maximum_section = heatmap_palette(heatmaps)

    fig = plt.figure(figsize=(17.2, 11.4))
    outer = fig.add_gridspec(
        2,
        1,
        height_ratios=[1.04, 1.00],
        hspace=0.22,
        left=0.050,
        right=0.920,
        top=0.940,
        bottom=0.070,
    )
    top = outer[0].subgridspec(
        1,
        3,
        width_ratios=[2.1, 3.4, 7.0],
        wspace=0.25,
    )
    bottom = outer[1].subgridspec(1, 3, wspace=0.24)

    def header(cell, letter: str, title: str):
        grid = cell.subgridspec(
            2,
            1,
            height_ratios=[0.12, 0.88],
            hspace=0.02,
        )
        ax_header = fig.add_subplot(grid[0])
        ax_header.axis("off")
        ax_header.text(
            0.00,
            0.50,
            letter,
            fontsize=13,
            fontweight="bold",
            va="center",
            ha="left",
        )
        ax_header.text(
            0.50,
            0.50,
            title,
            fontsize=10,
            fontweight="semibold",
            va="center",
            ha="center",
        )
        return grid[1], ax_header

    def bottom_header(cell, letter: str, title: str):
        grid = cell.subgridspec(
            3,
            1,
            height_ratios=[0.14, 0.66, 2.20],
            hspace=0.08,
        )
        ax_header = fig.add_subplot(grid[0])
        ax_header.axis("off")
        ax_header.text(
            0.00,
            0.50,
            letter,
            fontsize=13,
            fontweight="bold",
            va="center",
            ha="left",
        )
        ax_header.text(
            0.50,
            0.50,
            title,
            fontsize=10,
            fontweight="semibold",
            va="center",
            ha="center",
        )
        return grid[1], grid[2], ax_header

    # a) Flat and concave parallel/perpendicular geometries.
    a_content, _a_header = header(top[0, 0], "a", "Sectioning geometry")
    a_grid = a_content.subgridspec(2, 2, wspace=0.03, hspace=0.24)
    ax_a = [fig.add_subplot(a_grid[row, col]) for row in range(2) for col in range(2)]
    show_images(
        ax_a,
        [
            visual_paths["parallel"],
            visual_paths["perpendicular"],
            visual_paths["bowl_00"],
            visual_paths["bowl_perpendicular"],
        ],
        [
            "Flat · parallel\nX=0°",
            "Flat · perpendicular\nX=90°",
            "0.7 mm bowl · parallel\nX=0°",
            "0.7 mm bowl · perpendicular\nX=90°",
        ],
        title_fontsize=5.7,
    )

    # b) Face-aligned best-cut curves for two models and two section thicknesses.
    b_content, _b_header = header(
        top[0, 1],
        "b",
        "Maximum target volume per section",
    )
    ax_b = fig.add_subplot(b_content)
    ax_b.axvspan(0, 15, color=PARALLEL_BAND, alpha=0.12, lw=0)
    ax_b.axvspan(75, 90, color=PERPENDICULAR_BAND, alpha=0.12, lw=0)
    face_50, _diagonal_50 = curves[0.05]
    face_100, _diagonal_100 = curves[0.10]
    bowl_face_50, _bowl_diagonal_50 = bowl_curves[0.05]
    bowl_face_100, _bowl_diagonal_100 = bowl_curves[0.10]
    ax_b.plot(
        tilts,
        face_50,
        color=FLAT_50_COLOR,
        lw=2.0,
        ls="-",
    )
    ax_b.plot(
        tilts,
        face_100,
        color=FLAT_100_COLOR,
        lw=2.0,
        ls="-",
    )
    ax_b.plot(
        tilts,
        bowl_face_50,
        color=BOWL_50_COLOR,
        lw=2.0,
        ls="-",
    )
    ax_b.plot(
        tilts,
        bowl_face_100,
        color=BOWL_100_COLOR,
        lw=2.0,
        ls="-",
    )
    ax_b.set_xlim(0, 90)
    ax_b.set_ylim(0, 1.0)
    ax_b.set_xlabel("Tilt from the top face (°)")
    ax_b.set_ylabel("Target in best serial cut (mm³)")
    ax_b.legend(
        handles=[
            Line2D(
                [0],
                [0],
                color=FLAT_50_COLOR,
                lw=2.0,
                label="Cube–cuboid · 50 µm",
            ),
            Line2D(
                [0],
                [0],
                color=FLAT_100_COLOR,
                lw=2.0,
                label="Cube–cuboid · 100 µm",
            ),
            Line2D(
                [0],
                [0],
                color=BOWL_50_COLOR,
                lw=2.0,
                label="0.7 mm bowl · 50 µm",
            ),
            Line2D(
                [0],
                [0],
                color=BOWL_100_COLOR,
                lw=2.0,
                label="0.7 mm bowl · 100 µm",
            ),
        ],
        frameon=False,
        fontsize=6.1,
        loc="upper right",
        ncol=2,
        handlelength=2.6,
        columnspacing=0.9,
    )
    clean_axis(ax_b)

    # c) Six surface-indexed absolute-volume profiles.
    c_content, c_header = header(
        top[0, 2],
        "c",
        "Absolute volumes across consecutive sections",
    )
    c_header.texts[1].set_x(0.34)
    c_header.legend(
        handles=[
            Patch(facecolor=SECTION_BLUE, alpha=0.72, label="Total section"),
            Patch(facecolor=TARGET_ORANGE, label="Target in section"),
        ],
        frameon=False,
        fontsize=6.8,
        ncol=2,
        loc="center right",
    )
    c_grid = c_content.subgridspec(2, 3, hspace=0.46, wspace=0.13)
    profile_specs = [
        ("flat_parallel", "Flat · 0°"),
        ("flat_15", "Flat · X=15°"),
        ("flat_perpendicular", "Flat · 90°"),
        ("bowl_parallel", "0.7 mm bowl · 0°"),
        ("bowl_15", "0.7 mm bowl · X=15°"),
        ("bowl_perpendicular", "0.7 mm bowl · 90°"),
    ]
    x_max = max(len(strategy_profiles[key].section_numbers) for key, _ in profile_specs)
    c_axes: list[plt.Axes] = []
    for i, (key, title) in enumerate(profile_specs):
        row, col = divmod(i, 3)
        ax = fig.add_subplot(c_grid[row, col])
        plot_serial_bars(
            ax,
            strategy_profiles[key],
            title,
            x_max=x_max,
            show_xlabel=(row == 1),
        )
        if col == 0:
            ax.set_ylabel("Absolute volume (mm³)", fontsize=7)
        else:
            ax.tick_params(axis="y", labelleft=False, left=False)
        c_axes.append(ax)

    # d) Flat cube-cuboid X/Y rotations.
    d_visual_spec, d_heat_spec, _d_header = bottom_header(
        bottom[0, 0],
        "d",
        "Flat cube–cuboid: X–Y rotation",
    )
    d_visual = d_visual_spec.subgridspec(1, 2, wspace=0.03)
    ax_dv = [fig.add_subplot(d_visual[0, i]) for i in range(2)]
    show_images(
        ax_dv,
        [visual_paths["flat_00"], visual_paths["flat_1010"]],
        ["X=0°, Y=0°", "X=10°, Y=10°"],
    )
    ax_d = fig.add_subplot(d_heat_spec)
    mesh = plot_heatmap(
        ax_d,
        flat_xy,
        cmap,
        norm,
        "X rotation (°)",
        "Y rotation (°)",
        equal=True,
    )

    # e) Concavity and X rotation, including 10° dashed section boundaries.
    e_visual_spec, e_heat_spec, e_header = bottom_header(
        bottom[0, 1],
        "e",
        "Spherical concavity × X rotation",
    )
    e_header.texts[1].set_x(0.34)
    e_header.legend(
        handles=[
            Line2D(
                [0],
                [0],
                color="#1D4ED8",
                lw=1.1,
                ls="--",
                label="10° section boundaries",
            )
        ],
        loc="center right",
        frameon=False,
        fontsize=6.3,
        handlelength=2.2,
    )
    e_visual = e_visual_spec.subgridspec(1, 3, wspace=0.04)
    ax_ev = [fig.add_subplot(e_visual[0, i]) for i in range(3)]
    for index, (ax, apex) in enumerate(zip(ax_ev, (0.0, 0.5, 1.0))):
        shape = TOP_FLAT if np.isclose(apex, 0.0) else TOP_SPHERE_CONCAVE
        profile_0 = surface_first_profile(
            0.0,
            0.0,
            top_shape=shape,
            curve_offset=apex,
            target_mode=TARGET_FOLLOW,
        )
        profile_10 = surface_first_profile(
            10.0,
            0.0,
            top_shape=shape,
            curve_offset=apex,
            target_mode=TARGET_FOLLOW,
        )
        draw_morphology_cross_section(
            ax,
            apex,
            profile_0,
            profile_10,
            show_axis_reference=(index == 0),
        )
    ax_e = fig.add_subplot(e_heat_spec)
    plot_heatmap(
        ax_e,
        concavity_x,
        cmap,
        norm,
        "X rotation (°)",
        "Bowl depth (mm)",
        equal=False,
    )

    # f) X/Y rotations at 0.7 mm bowl depth.
    f_visual_spec, f_heat_spec, _f_header = bottom_header(
        bottom[0, 2],
        "f",
        "0.7 mm bowl: X–Y rotation",
    )
    f_visual = f_visual_spec.subgridspec(1, 2, wspace=0.03)
    ax_fv = [fig.add_subplot(f_visual[0, i]) for i in range(2)]
    show_images(
        ax_fv,
        [visual_paths["bowl_00"], visual_paths["bowl_1010"]],
        ["X=0°, Y=0°", "X=10°, Y=10°"],
    )
    ax_f = fig.add_subplot(f_heat_spec)
    plot_heatmap(
        ax_f,
        concave_xy,
        cmap,
        norm,
        "X rotation (°)",
        "Y rotation (°)",
        equal=True,
    )

    cax = fig.add_axes([0.940, 0.095, 0.011, 0.34])
    cbar = fig.colorbar(mesh, cax=cax)
    cbar.set_ticks(np.arange(minimum_section, maximum_section + 1))
    cbar.set_label("Earliest highest-yielding section number", fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    fig.legend(
        handles=[
            Patch(facecolor=CUBE_COLOR, alpha=0.20, edgecolor="none", label="Tissue"),
            Patch(facecolor=SECTION_BLUE, alpha=0.70, edgecolor="none", label="Section"),
            Patch(facecolor=TARGET_COLOR, alpha=0.42, edgecolor="none", label="Target"),
            Patch(facecolor=TARGET_ORANGE, edgecolor="none", label="Target in section"),
            Line2D(
                [0],
                [0],
                marker="x",
                color="#111827",
                lw=0,
                markersize=6,
                label="× = multiple sections tied for highest yield",
            ),
        ],
        loc="upper center",
        bbox_to_anchor=(0.50, 0.992),
        ncol=5,
        frameon=False,
        fontsize=8,
        columnspacing=1.25,
        handlelength=1.2,
    )

    png = OUT / f"{FIGURE_STEM}.png"
    pdf = OUT / f"{FIGURE_STEM}.pdf"
    temporary_png = OUT / f"_{FIGURE_STEM}_temporary.png"
    temporary_pdf = OUT / f"_{FIGURE_STEM}_temporary.pdf"
    fig.savefig(
        temporary_png,
        dpi=300,
        facecolor="white",
        bbox_inches="tight",
        pad_inches=0.08,
    )
    fig.savefig(
        temporary_pdf,
        facecolor="white",
        bbox_inches="tight",
        pad_inches=0.08,
    )
    plt.close(fig)
    if png.exists():
        png.unlink()
    if pdf.exists():
        pdf.unlink()
    temporary_png.replace(png)
    temporary_pdf.replace(pdf)
    print("wrote", png, flush=True)
    print("wrote", pdf, flush=True)
    return png, pdf


def validate(
    profiles: dict[str, SerialProfile],
    heatmaps: list[HeatmapData],
) -> None:
    parallel = profiles["flat_parallel"]
    perpendicular = profiles["flat_perpendicular"]
    if parallel.first_peak_section != 2:
        raise AssertionError("Reference parallel peak should begin in section 2")
    if not np.isclose(parallel.peak_target_mm3, 0.9, atol=1e-8):
        raise AssertionError("Reference parallel peak should be 0.9 mm3")
    if not np.isclose(perpendicular.peak_target_mm3, 0.06, atol=1e-8):
        raise AssertionError("Reference perpendicular peak should be 0.06 mm3")
    if int(np.sum(parallel.target_mm3 > 1e-9)) != 2:
        raise AssertionError("Parallel target should occupy two surface-aligned sections")
    if int(np.sum(perpendicular.target_mm3 > 1e-9)) != 30:
        raise AssertionError("Perpendicular target should occupy thirty sections")
    for data in heatmaps:
        if np.any(data.section_number < 1):
            raise AssertionError(f"{data.name}: invalid section number")
        if np.any(data.tie_count < 1):
            raise AssertionError(f"{data.name}: invalid tie count")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    print("Stage 1/5: strategy data", flush=True)
    tilts, curves = compute_best_cut_curves()
    bowl_curves = compute_bowl_best_cut_curves()
    bowl_block = TissueBlock(CUBE, TOP_SPHERE_CONCAVE, 0.7, "y")
    bowl_columns = make_columns(
        bowl_block,
        TARGET,
        DEPTH,
        TARGET_FOLLOW,
        n=GRID_N,
    )
    strategy_profiles = {
        "flat_parallel": surface_first_profile(0.0, 0.0),
        "flat_15": surface_first_profile(15.0, 0.0),
        "flat_perpendicular": surface_first_profile(90.0, 0.0),
        "bowl_parallel": surface_first_profile(
            0.0,
            0.0,
            top_shape=TOP_SPHERE_CONCAVE,
            curve_offset=0.7,
            target_mode=TARGET_FOLLOW,
            columns=bowl_columns,
        ),
        "bowl_15": surface_first_profile(
            15.0,
            0.0,
            top_shape=TOP_SPHERE_CONCAVE,
            curve_offset=0.7,
            target_mode=TARGET_FOLLOW,
            columns=bowl_columns,
        ),
        "bowl_perpendicular": surface_first_profile(
            90.0,
            0.0,
            top_shape=TOP_SPHERE_CONCAVE,
            curve_offset=0.7,
            target_mode=TARGET_FOLLOW,
            columns=bowl_columns,
        ),
    }
    write_curve_data(tilts, curves, bowl_curves)
    write_profile_data(strategy_profiles)

    print("Stage 2/5: parameter-sweep data", flush=True)
    flat_xy = compute_flat_xy()
    concavity_x = compute_concavity_x()
    concave_xy = compute_concave_xy(0.7)
    for data in (flat_xy, concavity_x, concave_xy):
        write_heatmap_data(data)
    write_metadata()
    validate(strategy_profiles, [flat_xy, concavity_x, concave_xy])

    print("Stage 3/5: model visualizations", flush=True)
    visual_paths = generate_visuals(
        strategy_profiles["flat_parallel"],
        strategy_profiles["flat_perpendicular"],
        strategy_profiles["bowl_perpendicular"],
        flat_xy,
        concavity_x,
        concave_xy,
    )

    print("Stage 4/5: final figure", flush=True)
    make_figure_v2(
        tilts,
        curves,
        bowl_curves,
        strategy_profiles,
        flat_xy,
        concavity_x,
        concave_xy,
        visual_paths,
    )
    print("Stage 5/5: complete", flush=True)


if __name__ == "__main__":
    main()
