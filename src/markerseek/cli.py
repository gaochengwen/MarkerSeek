"""CLI wiring for MarkerSeek."""

from __future__ import annotations

import argparse

from .analysis import MarkerSeekError, run_analysis
from .output import write_analysis_outputs


def positive_int(raw_value: str) -> int:
    try:
        value = int(raw_value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("Value must be an integer.") from exc
    if value <= 0:
        raise argparse.ArgumentTypeError("Value must be a positive integer.")
    return value


def nonnegative_int(raw_value: str) -> int:
    try:
        value = int(raw_value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("Value must be an integer.") from exc
    if value < 0:
        raise argparse.ArgumentTypeError("Value must be a non-negative integer.")
    return value


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="markerseek", description="Pi analysis for chloroplast GenBank files.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    analyze = subparsers.add_parser("analyze", help="Run the complete MarkerSeek workflow.")
    analyze.add_argument("inputs", nargs="+", help="GenBank files or directories.")
    analyze.add_argument("--outdir", default="output", help="Directory for PDF, PNG, and TSV outputs.")
    analyze.add_argument("--reference", help="Reference GenBank file. Defaults to the first input.")
    analyze.add_argument("--window", type=positive_int, default=600, help="Sliding window size in bp.")
    analyze.add_argument("--step", type=positive_int, default=200, help="Sliding step size in bp.")
    analyze.add_argument(
        "--hotspot-mode",
        choices=("top-percent", "top-n", "threshold"),
        default="top-percent",
        help="How hotspots are chosen.",
    )
    analyze.add_argument(
        "--hotspot-value",
        type=float,
        default=3.0,
        help="Top percentage, top count, or minimum Pi threshold depending on --hotspot-mode.",
    )
    analyze.add_argument("--mafft-bin", default="mafft", help="MAFFT executable or absolute path.")
    analyze.add_argument(
        "--mafft-threads",
        type=positive_int,
        default=None,
        help="Number of MAFFT worker threads. Defaults to MAFFT's own thread setting.",
    )
    analyze.add_argument(
        "--label-mode",
        choices=("peak-only", "all", "none"),
        default="peak-only",
        help="How hotspot labels are drawn on the figure.",
    )
    analyze.add_argument(
        "--label-max",
        type=positive_int,
        default=None,
        help="Maximum number of hotspot labels to draw (default: no limit).",
    )
    analyze.add_argument(
        "--label-min-distance",
        type=nonnegative_int,
        default=0,
        help="Minimum midpoint spacing in bp between labeled hotspots (default: 0 = label every peak that exceeds the threshold).",
    )
    analyze.add_argument(
        "--similarity-window",
        type=positive_int,
        default=200,
        help="Sliding window size in bp for the mVISTA-style similarity figure.",
    )
    analyze.add_argument(
        "--similarity-step",
        type=positive_int,
        default=60,
        help="Sliding step size in bp for the similarity figure.",
    )
    analyze.add_argument(
        "--similarity-floor",
        type=float,
        default=0.5,
        help="Lower bound (fraction in 0–1) of the similarity y-axis.",
    )
    analyze.add_argument(
        "--no-similarity-plot",
        action="store_true",
        help="Skip generating the similarity figure (similarity_plot.pdf/png).",
    )
    analyze.add_argument("--lsc", help="Manual region override in 1-based inclusive form start:end.")
    analyze.add_argument("--irb", help="Manual region override in 1-based inclusive form start:end.")
    analyze.add_argument("--ssc", help="Manual region override in 1-based inclusive form start:end.")
    analyze.add_argument("--ira", help="Manual region override in 1-based inclusive form start:end.")

    return parser


def parse_manual_region(raw_value: str) -> tuple[int, int]:
    if ":" not in raw_value:
        raise MarkerSeekError("Manual regions must use start:end syntax.")
    start_text, end_text = raw_value.split(":", maxsplit=1)
    try:
        start = int(start_text)
        end = int(end_text)
    except ValueError as exc:
        raise MarkerSeekError("Manual region coordinates must be integers in start:end format.") from exc
    if start < 1 or end < 1:
        raise MarkerSeekError("Manual region coordinates must be 1-based positive integers.")
    return start - 1, end


def manual_regions_from_args(args: argparse.Namespace) -> dict[str, tuple[int, int]] | None:
    region_text = {"LSC": args.lsc, "IRb": args.irb, "SSC": args.ssc, "IRa": args.ira}
    provided = {name: value for name, value in region_text.items() if value is not None}
    if not provided:
        return None
    if len(provided) != 4:
        raise MarkerSeekError("Manual region overrides require all of --lsc, --irb, --ssc, and --ira.")
    return {name: parse_manual_region(text) for name, text in region_text.items()}


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command != "analyze":
        parser.error("Unsupported command.")

    try:
        result = run_analysis(
            args.inputs,
            outdir=args.outdir,
            reference=args.reference,
            window=args.window,
            step=args.step,
            hotspot_mode=args.hotspot_mode,
            hotspot_value=args.hotspot_value,
            mafft_bin=args.mafft_bin,
            mafft_threads=args.mafft_threads,
            manual_regions=manual_regions_from_args(args),
        )
    except MarkerSeekError as exc:
        parser.exit(2, f"Error: {exc}\n")

    write_analysis_outputs(
        result,
        args.outdir,
        hotspot_mode=args.hotspot_mode,
        hotspot_value=args.hotspot_value,
        label_mode=args.label_mode,
        label_max=args.label_max,
        label_min_distance_bp=args.label_min_distance,
        similarity_window=args.similarity_window,
        similarity_step=args.similarity_step,
        similarity_floor=args.similarity_floor,
        include_similarity_plot=not args.no_similarity_plot,
    )
    return 0
