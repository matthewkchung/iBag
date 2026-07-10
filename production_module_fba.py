"""
Screen heterologous production modules in the Bacillus coagulans iBag597 model.

The module CSVs in ../ecoli_20_modules_sg use BiGG-style metabolite ids such as
atp_c and accoa_c. iBag597 uses ids such as ATP[c] and Acetyl_CoA[c], so this
script maps native metabolites explicitly and creates only non-native pathway
intermediates from the metabolite table.
"""

from __future__ import annotations

import argparse
import ast
import csv
import re
from zipfile import ZipFile
from dataclasses import dataclass
from pathlib import Path
from xml.etree import ElementTree as ET

import cobra
import pandas as pd
from cobra import Metabolite, Reaction


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MODEL = ROOT / "iBag" / "ModelFiles" / "iBag597.xml"
DEFAULT_METABOLITES = ROOT / "ecoli_20_modules_sg" / "metabolite_table.csv"
DEFAULT_REACTIONS = ROOT / "ecoli_20_modules_sg" / "reaction_table.csv"
DEFAULT_PATHWAYS = ROOT / "ecoli_20_modules_sg" / "pathway_table.csv"
DEFAULT_OUTPUT = ROOT / "iBag" / "production_module_screen.csv"
DEFAULT_IBAG_BIGG_MAP = ROOT / "iBag" / "iBag597BiGG.xlsx"


NATIVE_METABOLITE_MAP = {
    "ac_c": "Acetate[c]",
    "acald_c": "Acetaldehyde[c]",
    "accoa_c": "Acetyl_CoA[c]",
    "actp_c": "Acetylphosphate[c]",
    "adp_c": "ADP[c]",
    "akg_c": "2_Oxoglutarate[c]",
    "asp__L_c": "L_Aspartate[c]",
    "atp_c": "ATP[c]",
    "btcoa_c": "Butanoyl_CoA[c]",
    "co2_c": "CO2[c]",
    "coa_c": "CoA[c]",
    "etoh_c": "Ethanol[c]",
    "glu__L_c": "L_Glutamate[c]",
    "h2o_c": "H2O[c]",
    "h_c": "H+[c]",
    "hom__L_c": "L_Homoserine[c]",
    "lac__D_c": "D_Lactate[c]",
    "nad_c": "NAD[c]",
    "nadh_c": "NADH[c]",
    "nadp_c": "NADP[c]",
    "nadph_c": "NADPH[c]",
    "nh4_c": "NH3[c]",
    "oaa_c": "Oxaloacetate[c]",
    "phom_c": "O_Phospho_L_homoserine[c]",
    "pi_c": "Phosphate[c]",
    "ppcoa_c": "Propanoyl_CoA[c]",
    "pyr_c": "Pyruvate[c]",
    "succ_c": "Succinate[c]",
    "succoa_c": "Succinyl_CoA[c]",
    "thr__L_c": "L_Threonine[c]",
}


REACTION_ARROW_RE = re.compile(r"\s*(<=>|=>|-->)\s*")
TERM_RE = re.compile(r"^\s*(?:(-?\d+(?:\.\d+)?)\s+)?([A-Za-z0-9_]+_[a-z])\s*$")
REACTION_SPAN_RE = re.compile(
    r"[A-Za-z0-9_]+_[a-z](?:\s*\+\s*(?:-?\d+(?:\.\d+)?\s+)?[A-Za-z0-9_]+_[a-z])*"
    r"\s*(?:<=>|=>|-->)\s*"
    r"(?:-?\d+(?:\.\d+)?\s+)?[A-Za-z0-9_]+_[a-z](?:\s*\+\s*(?:-?\d+(?:\.\d+)?\s+)?[A-Za-z0-9_]+_[a-z])*"
)

XLSX_NS = {"a": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}


@dataclass
class ScreenResult:
    module_id: str
    name: str
    product_metabolite: str
    demand_reaction: str
    added_reactions: int
    added_metabolites: int
    imbalanced_reactions: int
    growth: float | None
    max_product: float | None
    product_at_min_growth: float | None
    status: str
    notes: str


def parse_list(value: object) -> list[str]:
    if value is None or pd.isna(value):
        return []
    try:
        parsed = ast.literal_eval(str(value))
    except (SyntaxError, ValueError):
        return [str(value)]
    if isinstance(parsed, (list, tuple)):
        return [str(item) for item in parsed]
    return [str(parsed)]


def load_metabolite_table(path: Path) -> dict[str, dict[str, object]]:
    df = pd.read_csv(path)
    return {str(row["id"]): row.to_dict() for _, row in df.iterrows()}


def _xlsx_shared_strings(archive: ZipFile) -> list[str]:
    if "xl/sharedStrings.xml" not in archive.namelist():
        return []
    root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
    return [
        "".join(text.text or "" for text in item.findall(".//a:t", XLSX_NS))
        for item in root.findall("a:si", XLSX_NS)
    ]


def _xlsx_cell_value(cell: ET.Element, shared_strings: list[str]) -> str:
    value = cell.find("a:v", XLSX_NS)
    if value is None or value.text is None:
        return ""
    if cell.attrib.get("t") == "s":
        return shared_strings[int(value.text)]
    return value.text


def _column_index(cell_ref: str) -> int:
    letters = re.match(r"[A-Z]+", cell_ref)
    if not letters:
        return 0
    index = 0
    for letter in letters.group(0):
        index = index * 26 + ord(letter) - ord("A") + 1
    return index - 1


def read_first_sheet_xlsx(path: Path) -> list[dict[str, str]]:
    """Read the first sheet from a small xlsx file without requiring openpyxl."""
    with ZipFile(path) as archive:
        shared_strings = _xlsx_shared_strings(archive)
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        first_sheet = workbook.find("a:sheets/a:sheet", XLSX_NS)
        if first_sheet is None:
            return []

        rels_root = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        rels = {rel.attrib["Id"]: rel.attrib["Target"] for rel in rels_root}
        rel_id = first_sheet.attrib["{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"]
        sheet_path = "xl/" + rels[rel_id].lstrip("/")
        sheet = ET.fromstring(archive.read(sheet_path))

        rows = []
        for row in sheet.findall(".//a:sheetData/a:row", XLSX_NS):
            values: list[str] = []
            for cell in row.findall("a:c", XLSX_NS):
                index = _column_index(cell.attrib.get("r", "A1"))
                while len(values) <= index:
                    values.append("")
                values[index] = _xlsx_cell_value(cell, shared_strings).strip()
            rows.append(values)

    if not rows:
        return []
    headers = rows[0]
    return [dict(zip(headers, row + [""] * (len(headers) - len(row)))) for row in rows[1:]]


def sanitize_reaction_string(value: object) -> str:
    text = str(value or "").strip().rstrip(",")
    matches = REACTION_SPAN_RE.findall(text)
    if matches:
        return matches[-1].strip().rstrip(",")
    return text


def sanitize_reaction_id(value: object) -> str:
    return str(value or "").strip().split(",", 1)[0].strip()


def clean_annotation_value(value: object) -> str:
    text = str(value or "").strip()
    return "" if text.upper() == "NONE" else text


def load_heterologous_reaction_table(path: Path) -> pd.DataFrame:
    """
    Extract reactions whose iBag abbreviation is NONE from iBag597BiGG.xlsx.

    The workbook header spells this column "Abbrevation (iBag)", so column
    matching is normalized to avoid depending on the typo.
    """
    rows = read_first_sheet_xlsx(path)
    if not rows:
        return pd.DataFrame(columns=["id", "name", "rxn_str", "kegg_id", "bigg_id"])

    def normalized_key(row: dict[str, str], wanted: str) -> str:
        wanted_norm = re.sub(r"[^a-z0-9]", "", wanted.lower())
        for key in row:
            if re.sub(r"[^a-z0-9]", "", key.lower()) == wanted_norm:
                return key
        raise KeyError(f"Column not found in {path}: {wanted}")

    sample = rows[0]
    id_col = normalized_key(sample, "Abbreviation csv")
    name_col = normalized_key(sample, "name")
    rxn_col = normalized_key(sample, "rxn csv")
    ibag_col = normalized_key(sample, "Abbrevation iBag")
    kegg_col = normalized_key(sample, "KEGG id")
    bigg_col = normalized_key(sample, "BiGG id")

    reaction_rows = []
    for row in rows:
        if row.get(ibag_col, "").strip().upper() != "NONE":
            continue
        reaction_rows.append(
            {
                "id": sanitize_reaction_id(row.get(id_col, "")),
                "name": row.get(name_col, "").strip(),
                "rxn_str": sanitize_reaction_string(row.get(rxn_col, "")),
                "kegg_id": clean_annotation_value(row.get(kegg_col, "")),
                "bigg_id": clean_annotation_value(row.get(bigg_col, "")),
            }
        )
    return pd.DataFrame(reaction_rows)


def model_metabolite_id(bigg_id: str) -> str:
    if bigg_id in NATIVE_METABOLITE_MAP:
        return NATIVE_METABOLITE_MAP[bigg_id]
    if bigg_id.endswith("_c"):
        return f"{bigg_id[:-2]}[c]"
    if bigg_id.endswith("_e"):
        return f"{bigg_id[:-2]}[e]"
    return bigg_id


def ensure_metabolite(model: cobra.Model, bigg_id: str, met_table: dict[str, dict[str, object]]) -> Metabolite:
    mapped_id = model_metabolite_id(bigg_id)
    if mapped_id in model.metabolites:
        return model.metabolites.get_by_id(mapped_id)

    bare_id = bigg_id.rsplit("_", 1)[0]
    row = met_table.get(bare_id, {})
    charge = row.get("charge")
    if charge is not None and not pd.isna(charge):
        charge = int(float(charge))
    else:
        charge = None

    met = Metabolite(
        mapped_id,
        name=str(row.get("name") or bare_id),
        formula=None if pd.isna(row.get("formula")) else row.get("formula"),
        charge=charge,
        compartment="c" if bigg_id.endswith("_c") else "e",
    )
    kegg_id = row.get("kegg_id")
    if kegg_id is not None and not pd.isna(kegg_id) and str(kegg_id).strip():
        met.annotation["kegg.compound"] = str(kegg_id).strip()
    model.add_metabolites([met])
    return met


def parse_reaction_stoichiometry(
    rxn_str: str,
    model: cobra.Model,
    met_table: dict[str, dict[str, object]],
) -> tuple[dict[Metabolite, float], bool]:
    parts = REACTION_ARROW_RE.split(rxn_str, maxsplit=1)
    if len(parts) != 3:
        raise ValueError(f"Cannot parse reaction string: {rxn_str}")
    left, arrow, right = parts
    reversible = arrow == "<=>"
    stoich: dict[Metabolite, float] = {}

    def add_side(side: str, sign: float) -> None:
        for term in side.split("+"):
            term = term.strip()
            if not term:
                continue
            match = TERM_RE.match(term)
            if not match:
                raise ValueError(f"Cannot parse reaction term '{term}' in '{rxn_str}'")
            coeff = float(match.group(1) or 1.0)
            met = ensure_metabolite(model, match.group(2), met_table)
            stoich[met] = stoich.get(met, 0.0) + sign * coeff

    add_side(left, -1.0)
    add_side(right, 1.0)
    return stoich, reversible


def add_module_reactions(
    model: cobra.Model,
    module_rxns: list[str],
    rxn_df: pd.DataFrame,
    met_table: dict[str, dict[str, object]],
) -> tuple[int, int]:
    before_mets = len(model.metabolites)
    added_rxns = 0
    rxn_rows = rxn_df.set_index("id")
    for reaction_id in module_rxns:
        if reaction_id in model.reactions:
            continue
        if reaction_id not in rxn_rows.index:
            raise KeyError(f"Reaction {reaction_id} is not present in the reaction table")
        row = rxn_rows.loc[reaction_id]
        stoich, reversible = parse_reaction_stoichiometry(row["rxn_str"], model, met_table)
        reaction = Reaction(reaction_id)
        reaction.name = str(row.get("name") or reaction_id)
        reaction.lower_bound = -1000.0 if reversible else 0.0
        reaction.upper_bound = 1000.0
        reaction.add_metabolites(stoich)
        ec_number = row.get("ec_number")
        if ec_number is not None and not pd.isna(ec_number) and str(ec_number).strip():
            reaction.annotation["ec-code"] = str(ec_number).strip()
        kegg_id = None if row.get("kegg_id") is None or pd.isna(row.get("kegg_id")) else clean_annotation_value(row.get("kegg_id"))
        if kegg_id:
            reaction.annotation["kegg.reaction"] = kegg_id
        bigg_id = None if row.get("bigg_id") is None or pd.isna(row.get("bigg_id")) else clean_annotation_value(row.get("bigg_id"))
        if bigg_id:
            reaction.annotation["bigg.reaction"] = bigg_id
        model.add_reactions([reaction])
        added_rxns += 1
    return added_rxns, len(model.metabolites) - before_mets


def add_all_heterologous_reactions(
    model: cobra.Model,
    workbook: Path,
    met_table: dict[str, dict[str, object]],
) -> tuple[pd.DataFrame, int, int]:
    rxn_df = load_heterologous_reaction_table(workbook)
    added_rxns, added_mets = add_module_reactions(model, rxn_df["id"].tolist(), rxn_df, met_table)
    return rxn_df, added_rxns, added_mets


def add_product_demand(model: cobra.Model, product_id: str, met_table: dict[str, dict[str, object]]) -> str:
    product_met = ensure_metabolite(model, f"{product_id}_c", met_table)
    demand_id = f"DM_{product_met.id.replace('[', '_').replace(']', '')}"
    if demand_id in model.reactions:
        return demand_id
    demand = model.add_boundary(product_met, type="demand", reaction_id=demand_id)
    demand.lower_bound = 0.0
    demand.upper_bound = 1000.0
    return demand.id


def maximize_reaction(model: cobra.Model, reaction_id: str) -> float | None:
    with model:
        model.objective = reaction_id
        solution = model.optimize()
        if solution.status != "optimal":
            return None
        return float(solution.objective_value)


def count_imbalanced_reactions(model: cobra.Model, reaction_ids: list[str]) -> int:
    imbalanced = 0
    for reaction_id in reaction_ids:
        if reaction_id not in model.reactions:
            continue
        if model.reactions.get_by_id(reaction_id).check_mass_balance():
            imbalanced += 1
    return imbalanced


def screen_module(
    base_model: cobra.Model,
    base_growth: float,
    pathway_row: pd.Series,
    rxn_df: pd.DataFrame,
    met_table: dict[str, dict[str, object]],
    min_growth_fraction: float,
) -> ScreenResult:
    model = base_model.copy()
    module_id = str(pathway_row["id"])
    module_rxns = parse_list(pathway_row["rxns"])
    product_id = str(pathway_row["product_id"])
    notes = []

    try:
        added_rxns, added_mets = add_module_reactions(model, module_rxns, rxn_df, met_table)
        demand_id = add_product_demand(model, product_id, met_table)
    except Exception as exc:
        return ScreenResult(module_id, str(pathway_row["name"]), product_id, "", 0, 0, 0, None, None, None, "error", str(exc))

    imbalanced_reactions = count_imbalanced_reactions(model, module_rxns)
    if imbalanced_reactions:
        notes.append(f"{imbalanced_reactions} module reaction(s) fail mass/charge balance check")

    growth_solution = model.optimize()
    growth = None
    if growth_solution.status == "optimal":
        growth = float(growth_solution.objective_value)
    else:
        notes.append(f"growth status={growth_solution.status}")

    max_product = maximize_reaction(model, demand_id)

    product_at_min_growth = None
    if base_growth > 0:
        with model:
            model.reactions.get_by_id("EXBiomass").lower_bound = base_growth * min_growth_fraction
            model.objective = demand_id
            solution = model.optimize()
            if solution.status == "optimal":
                product_at_min_growth = float(solution.objective_value)
            else:
                notes.append(f"min-growth product status={solution.status}")

    return ScreenResult(
        module_id=module_id,
        name=str(pathway_row["name"]),
        product_metabolite=model_metabolite_id(f"{product_id}_c"),
        demand_reaction=demand_id,
        added_reactions=added_rxns,
        added_metabolites=added_mets,
        imbalanced_reactions=imbalanced_reactions,
        growth=growth,
        max_product=max_product,
        product_at_min_growth=product_at_min_growth,
        status="ok",
        notes="; ".join(notes),
    )


def write_results(results: list[ScreenResult], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(ScreenResult.__dataclass_fields__))
        writer.writeheader()
        for result in results:
            writer.writerow(result.__dict__)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--metabolites", type=Path, default=DEFAULT_METABOLITES)
    parser.add_argument("--reactions", type=Path, default=DEFAULT_REACTIONS)
    parser.add_argument("--pathways", type=Path, default=DEFAULT_PATHWAYS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--ibag-bigg-map", type=Path, default=DEFAULT_IBAG_BIGG_MAP)
    parser.add_argument(
        "--load-all-heterologous",
        action="store_true",
        help="Preload every reaction in iBag597BiGG.xlsx whose iBag abbreviation is NONE.",
    )
    parser.add_argument("--min-growth-fraction", type=float, default=0.9)
    args = parser.parse_args()

    base_model = cobra.io.read_sbml_model(str(args.model))
    base_solution = base_model.optimize()
    if base_solution.status != "optimal":
        raise RuntimeError(f"Base model growth optimization failed: {base_solution.status}")
    base_growth = float(base_solution.objective_value)
    met_table = load_metabolite_table(args.metabolites)
    rxn_df = pd.read_csv(args.reactions)
    pathway_df = pd.read_csv(args.pathways)

    if args.load_all_heterologous:
        heterologous_rxn_df, added_rxns, added_mets = add_all_heterologous_reactions(
            base_model,
            args.ibag_bigg_map,
            met_table,
        )
        heterologous_rxns = heterologous_rxn_df["id"].tolist()
        imbalanced = count_imbalanced_reactions(base_model, heterologous_rxns)
        print(
            f"Loaded {added_rxns}/{len(heterologous_rxns)} heterologous reactions "
            f"from {args.ibag_bigg_map} ({added_mets} new metabolites, "
            f"{imbalanced} imbalanced reactions)"
        )

    results = [
        screen_module(base_model, base_growth, row, rxn_df, met_table, args.min_growth_fraction)
        for _, row in pathway_df.iterrows()
    ]
    results.sort(
        key=lambda result: (
            result.product_at_min_growth is not None,
            result.product_at_min_growth or -1.0,
            result.max_product or -1.0,
        ),
        reverse=True,
    )
    write_results(results, args.output)

    print(f"Wrote {len(results)} module results to {args.output}")
    print(f"Wild-type growth objective: {base_growth:.6g}")
    print("Top modules by product flux at minimum wild-type growth:")
    for result in results[:10]:
        print(
            f"{result.module_id:14s} growth={result.growth!s:>10s} "
            f"max_product={result.max_product!s:>10s} "
            f"product_at_{args.min_growth_fraction:.0%}_WT_growth={result.product_at_min_growth!s:>10s} "
            f"imbalanced_rxns={result.imbalanced_reactions}"
        )


if __name__ == "__main__":
    main()
